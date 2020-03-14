# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 10:43:00 2019

@author: Jonathan White
"""

import os
from fitDots import dotFinder
import skimage
import networkx as nx
import numpy as np
from operator import itemgetter
import json
import pickle
from matplotlib import pyplot as plt
import pandas as pd
#import seaborn as sns
import copy
from scipy.spatial import cKDTree as KDTree
from bisect import bisect_left, bisect
from multiprocessing import Pool

'''
def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    # if i != len(a) and a[i] == x:
    return i
    # raise ValueError
'''

class dotAligner:
    '''
    Class for aligning images using dots that occur in each image.
    '''

    def __init__(self, fname_base='MMStack_Pos%d.ome.tif', fname_nums=0, directory=None, z=0, channel=1):
        '''
        :param fname_base: String text base name of files to load. May include the formating operator, '%', to insert numbers.
        :param fname_nums: tuple of number to insert in fname_base using the formating operator, '%'
        :param directory: the directory in which folders containing files for each hyb and background are stored.
        :param z: the z index of stack to analyze
        :param channel: the channel of the stack to analyze
        '''
        # self.hybs = []
        self.graphs = []
        self.PSFs = []
        self.psfSorted = False
        self.directory = None
        self.z = z
        self.channel = channel
        self.imfname = fname_base % fname_nums
        self.alignedPSFs = None

        if directory:
            self.directory = directory

    def findPSFs(self, nhybs=80, bnd=2048, ref=False, threshold=None):
        '''
        Finds psfs in all hybs
        :param nhybs: number of hybs in which to find PSFs; during testing use lower number to speed computation
        :param bnd: for testing: choose a subset of the image up to bnd pixels to reduce runtime
        :param ref: If true, find bead PSFs
        :param threshold: If provided reject local maxima with counts fewer than this
        :return:
        '''

        os.chdir(self.directory)

        # get background data
        os.chdir('final_background')
        bg = skimage.io.imread(self.imfname)
        # self.bg = bg[0][:bnd,:bnd,1]
        self.bg = bg[self.z][:bnd, :bnd, self.channel]
        os.chdir('..')
        if ref:
            os.chdir('beads')
            im = skimage.io.imread(self.imfname)
            im = im[self.z][:bnd, :bnd, self.channel]
            mod = dotFinder(im, self.bg, 2)
            model_im, singlePSFs, overlapping_peaks = mod.fitSinglePSFs(threshold=threshold)
            self.refPSFs = singlePSFs
            return self.refPSFs, mod

        # read in hyb data, find dots
        mods = []
        for hyb in range(nhybs):
            os.chdir('HybCycle_%d' % hyb)
            im = skimage.io.imread(self.imfname)
            # im = im[0][:,:,1]
            # im = im[0][:bnd,:bnd,1]
            im = im[self.z][:bnd, :bnd, self.channel]
            # im[im < threshold] = 0
            mod = dotFinder(im, self.bg, 2)
            mods.append(mod)
            model_im, singlePSFs, overlapping_peaks = mod.fitSinglePSFs(threshold=threshold)
            self.PSFs.append(singlePSFs)
            # self.hybs.append(mod)
            os.chdir('..')
        os.chdir('..')
        return self.PSFs, mods

    def savePSFs(self, filename):
        '''
        Saves PSFs found using findPSFs
        :param filename: file name to which to save PSFss
        :return:
        '''
        os.chdir(self.directory)
        with open(filename, 'wb') as f:
            pickle.dump(self.PSFs, f)

    def savePSFsDF(self, filename):
        self.psfsDF.to_csv(filename)

    def saveRefPSFs(self, filename):
        os.chdir(self.directory)
        with open(filename, 'wb') as f:
            pickle.dump(self.refPSFs, f)

    def loadPSFs(self, filename):
        '''
        Loads saved PSFs computed earlier by findPSFs
        :param filename: filename of saved PSFs
        :return:
        '''

        # os.chdir(self.directory)
        if filename[-3:] == 'pkl':
            with open(filename, 'rb') as f:
                self.PSFs = pickle.load(f)
        else:
            self.psfsDF = pd.read_csv(filename, sep=',', header=0, index_col='Hyb')


    def find_matching_dots2(self, maxError=0.5, nHybs=40):
        '''
        Find dots that match each other between hybs in
        :param maxError:
        :return:
        '''
        self.spatMatchedDots = spatMatDotHolder()
        spatMatDot.nHybs = nHybs

        # for each hyb
        for hybNum, hybData in self.alignedPSFs.groupby(level=0):
            print('Hyb:', hybNum)
            self.spatMatchedDots.set_current_hyb(hybNum)

            # find a home for each PSF in the hyb
            for i, psf in hybData.iterrows():
                if hybNum == 0:
                    self.spatMatchedDots.insert_new_spatDot(spatMatDot(hybNum, psf))
                else:
                    found = self.spatMatchedDots.matchPSF(psf)
                    if not found:
                        self.spatMatchedDots.insert_new_spatDot(spatMatDot(hybNum, psf))

    def find_matching_dots(self, refHyb=0, search_dist=0.2):
        '''
        matches the dots in each aligned hybridization
        :refHyb: the integer number of the hyb to use as a reference and to look
                for matches to
        :dist: the maximum distance within which to match dots.
        :return:
        self.matchedDots contains sublists for each dot found in hyb0 with its
            matches in other hybs. Entries the sublists for each dot are either two tuples containing
            the rank of a the dot (by CompCoord) in its hyb in the 0th index, and the fitPSF object in the 1st index.
            None if no dots were found, or list of candidate dot two-tuples if ambiguous
            The first entry in each sublist is the dot found in in hyb0
            Each subsequent entry in each sublist after
            the original dot in hyb 0 is either a single dot that was matched; a list
            of candidate dots, if there there was not an unambiguous match; or None
            if no matches were found.
        '''
        self.matchedDots = []

        if not self.psfSorted:
            # self.sortedPSFs = self.alignedPSFs.sort_values(by='row')
            self.sortPSFCoord()
            self.psfSorted = True

        # Now match dots.

        # For each dot in the refHyb
        for i, psfi in enumerate(self.sortedPSFs[refHyb]):
            # for i, psfi in self.sortedPSFs.loc[0].iterrows():
            # Initialize sum of coordintates with which the average
            # of the coordinates will be calculated on each round
            # row_ave = psfi.row
            # col_ave = psfi.col

            imatches = [(i, psfi)]

            # for all of the other hybs:
            for j, hyb in enumerate(self.sortedPSFs):  # [1:]):
                # for j, hyb in self.sortedPSFs.groupby(level=0):
                if j == refHyb:
                    # if j == 0:
                    continue
                # find the highest and lowest index of dots within the search distnace of psfi
                #lo = index(self.sortedPSFKeys[j], (psfi.row - search_dist, psfi.col))
                lo = bisect_left(self.sortedPSFKeys[j], (psfi.row - search_dist, psfi.col))
                #hi = index(self.sortedPSFKeys[j], (psfi.row + search_dist, psfi.col))
                hi = bisect_left(self.sortedPSFKeys[j], (psfi.row + search_dist, psfi.col))
                # lo = index(self.sortedPSFKeys[j], psfi['row']-search_dist)
                # hi = index(self.sortedPSFKeys[j], psfi['row']+search_dist)
                # lo = index(np.array(hyb['row']), float(psfi['row']) - search_dist)
                # hi = index(np.array(hyb['row']), float(psfi['row']) + search_dist)

                # initialize list of dots within 1 pixel width of psfi
                candidates = []

                # for each dot within 1 row of psf 0
                for k, dot in enumerate(hyb[lo:hi]):
                    # for ind, dot in hyb.iloc[lo:(hi+1)].iterrows():
                    # if the dot is within the search distance of psfi, add it to the candidates list
                    if (dot.row - psfi.row) ** 2 + (dot.col - psfi.col) ** 2 < search_dist:
                        candidates.append((k + lo, dot))

                # now decide what to do with candidates

                if len(candidates) == 1:
                    # Hooray! We have a match
                    if candidates[0][1].matched:
                        print("Warning: PSF at %d, %d already matched." % (candidates[0][1].row, candidates[0][1].col))
                    candidates[0][1].matched = True
                    imatches.append(candidates[0])

                elif len(candidates) > 1:

                    # remove any candidates have already been matched
                    ic = 0
                    removed = []
                    while ic < len(candidates):
                        if candidates[ic][1].matched:
                            removed.append(candidates.pop(ic))
                        else:
                            ic += 1

                    # if this shortens the list of candidates to 1, then add it, mark as matched
                    if len(candidates) == 1:
                        if candidates[0][1].matched:
                            print("Warning: PSF at %d, %d already matched." % (
                            candidates[0][1].row, candidates[0][1].col))
                        imatches.append(candidates[0])
                        candidates[0][1].matched = True

                    # if their are still multiple candidates, add the list
                    elif len(candidates) > 1:
                        imatches.append(candidates)

                    # if all were removed, then just add originial candidates
                    else:
                        imatches.append(removed)
                        print("Warning: multiple already matched PSFs matched to psf at %d, %d." % (psfi.row, psfi.col))

                # Otherwise, no matches :(
                else:
                    imatches.append(None)
            # add matches to list of all matched dots
            self.matchedDots.append(imatches)

        # if there are multiple dots within 1 pixel of each other, try clustering by coordinates and brightness

        return


    def sortPSFCoord(self):
        '''
        sorts PSFs, makes shift corrections
        :return:
        self.sortedPSFs- list that contains list of PSFs in each hyb. The PSF lists for each hyb are sorted by their row
        then their column, that is by row than column coordinate. The coordinates are are aligned using the graph based
            alignnment.
        self.sortedPSFKeys - tuples of row, column for each PSF in sorted PSFs.
        '''

        # first, sort dots in each hyb by coordinate
        self.sortedPSFs = []
        self.sortedPSFKeys = []
        # self.graphAlignedDots = []
        for hybi, hybPSFs in enumerate(self.PSFs):
            # now, sort
            shybPSFs = sorted(hybPSFs, key=lambda psf: (psf.row, psf.col))

            # make shift correction
            # for i, hyb in enumerate(shybPSFs):
            for dot in shybPSFs:
                # dot.row -= self.mean_row_shifts[hybi]
                # dot.col -= self.mean_col_shifts[hybi]
                dot.row -= self.offsets['row_offset'].iloc(0)[hybi]
                dot.col -= self.offsets['col_offset'].iloc(0)[hybi]

            # build list of keys and initialize matched status of each psf
            shybPSFsKeys = [(psf.row, psf.col) for psf in shybPSFs]
            for psf in shybPSFs: psf.matched = False


            self.sortedPSFs.append(shybPSFs)
            self.sortedPSFKeys.append(shybPSFsKeys)

    def plotHybPSFs(self, hyb, log=True):
        '''
        Plots Fit PSFs on an image
        :param hyb: hyb number of the image
        :return:
        '''
        hybPSFS = self.PSFs[hyb]
        os.chdir(self.directory)
        os.chdir('HybCycle_%d' % hyb)
        im = skimage.io.imread(self.imfname)
        im = im[self.z][:, :, self.channel]
        if log:
            plt.imshow(np.log2(im), cmap='Greys')
            scale = 'Log2 Intensity (A.U.)'
        else:
            plt.imshow(im, cmap='Greys')
            scale = 'Intensity (A.U)'
        xpsfs = []
        ypsfs = []
        for psf in hybPSFS:
            xpsfs.append(psf.col)
            ypsfs.append(psf.row)
        plt.plot(xpsfs, ypsfs, '.r')
        cbar = plt.colorbar()
        cbar.set_label(scale)
        plt.legend()
        plt.show()
        os.chdir('..')

    def loadRefPSFs(self, refPSFs, df=True):
        '''
        Load File of saved bead psfs
        :param refPSFs: string filename of saved PSFs
        :return:
        '''
        #os.chdir(self.directory)
        if df:
            self.refPSFs = pd.read_csv(refPSFs)
        else:
            with open(refPSFs, 'rb') as f:
                refPSFs = pickle.load(f)
                self.refPSFs = self.convtPSFList2DF(refPSFs)

    def find_offsets_from_ref(self, maxXYSearchError=1, maxZSearchError=3, maxXYMatchError=0.1, maxZMatchError=2,
                              trav=True, minMatches=3, refFirst=False, brightThresh=None, stopRefProp=0.5, nLongestCheck=1000, hybsToAlign = None):
        '''
        Align images using bead
        :param maxXYMatchError: maximum allowed difference between matched displacement vectors between beads in reference
        image and in hyb image
        :return:
        data frame of offsets and offset standard errors
        '''
        # self.convtPSFList2DF()
        #print('Initiaizing Reference Graph')
        if refFirst:
            self.refGraph = refDotGraph(self.psfsDF.loc[0], brightThresh=brightThresh, maxError=maxXYMatchError)
            self.matchedRefDots = self.psfsDF.loc[0].copy
        else:
            self.refGraph = refDotGraph(self.refPSFs, brightThresh=brightThresh, maxError=maxXYMatchError)
            self.matchedRefDots = self.refPSFs.copy()
        # self.refGraph = refDotGraph(self.refPSFs)
        if not trav:
            self.refGraph.sortEdgesRows()

        hybs = []
        row_offset = []
        col_offset = []
        z_offset = []
        row_offset_SE = []
        col_offset_SE = []
        z_offset_SE = []
        n_offsets_calculated = []
        dataArrays = [hybs, row_offset, col_offset, z_offset, row_offset_SE, col_offset_SE, z_offset_SE, n_offsets_calculated]

        self.matchesDF = pd.DataFrame(columns=('ref_row', 'ref_col', 'ref_z', 'comp_row', 'comp_col', 'comp_z',
                                               'aligned_row', 'aligned_col', 'aligned_z'))

        if hybsToAlign:
            for hyb in hybsToAlign:
                self.alignHyb(dataArrays, hyb, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError,
                              trav, minMatches, brightThresh, stopRefProp, nLongestCheck)
        else:
            # for i, hybPSFs in enumerate(self.PSFs):
            for hyb, hybPSFs in self.psfsDF.groupby(level=0):
                self.alignHyb(dataArrays, hyb, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError,
                              trav, minMatches,brightThresh, stopRefProp, nLongestCheck)

        self.offsets = pd.DataFrame()
        self.offsets['hyb'] = hybs
        self.offsets['row_offset'] = row_offset
        self.offsets['col_offset'] = col_offset
        self.offsets['z_offset'] = z_offset
        self.offsets['row_offset_SE'] = row_offset_SE
        self.offsets['col_offset_SE'] = col_offset_SE
        self.offsets['z_offset_SE'] = z_offset_SE
        self.offsets['n_offsets_calculated'] = n_offsets_calculated
        return self.offsets, self.matchesDF

    def alignHyb(self, dataArrays, hyb, maxXYSearchError=1, maxZSearchError=3, maxXYMatchError=0.1, maxZMatchError=2,
                              trav=True, minMatches=3, brightThresh=None, stopRefProp=0.5, nLongestCheck=1000):
        # print('Initializing Hyb %d Graph' % i)
        # hybGraph = dotGraph(self.PSFs[i])
        hybs, row_offset, col_offset, z_offset, row_offset_SE, col_offset_SE, z_offset_SE, n_offsets_calculated = dataArrays
        print('Aligning hyb %d' % hyb)
        hybPSFs = self.psfsDF.loc[hyb]

        if trav:
            offsets, matchesDF = self.refGraph.alignWith2(hybPSFs, maxXYSearchError, maxZSearchError, maxXYMatchError,
                                                          maxZMatchError, minMatches=minMatches,
                                                          stopRefProp=stopRefProp, hyb_num=hyb,
                                                          nLongestCheck=nLongestCheck)
            r_mean, c_mean, z_mean, r_mean_se, c_mean_se, z_mean_se, r_weight_mean, c_weight_mean, n_matched = offsets
        else:
            hybGraph = dotGraph(hybPSFs, brightThresh)  # self.PSFs[i])
            hybGraph.sortEdgesRows()
            offsets, matchesDF = self.refGraph.alignWith(hybGraph, maxXYMatchError, maxZMatchError,
                                                         minMatches=minMatches)
            r_mean, c_mean, z_mean, r_mean_se, c_mean_se, z_mean_se, r_weight_mean, c_weight_mean, n_matched = offsets
        print('Offsets and SEs:', r_mean, c_mean, z_mean, r_mean_se, c_mean_se, z_mean_se)
        print('Calculated from %d matches' % n_matched)

        hybs.append(hyb)
        row_offset.append(r_mean)
        col_offset.append(c_mean)
        z_offset.append(z_mean)
        row_offset_SE.append(r_mean_se)
        col_offset_SE.append(c_mean_se)
        z_offset_SE.append(z_mean_se)
        n_offsets_calculated.append(n_matched)
        self.matchesDF = self.matchesDF.append(matchesDF)
        return

    def find_offsets_from_ref_par(self, maxXYSearchError=1, maxZSearchError=3, maxXYMatchError=0.1, maxZMatchError=2,
                                  minMatches=3, refFirst=False, brightThresh=None, stopRefProp=0.5, nLongestCheck=1000):
        '''
        Align images using bead
        :param maxXYError: maximum allowed difference between matched displacement vectors between beads in reference
        image and in hyb image
        :return:
        data frame of offsets and offset standard errors
        '''
        # self.convtPSFList2DF()
        # print('Initiaizing Reference Graph')
        if refFirst:
            self.refGraph = refDotGraph(self.psfsDF.loc[0], brightThresh=brightThresh, maxError=maxXYMatchError)
            self.matchedRefDots = self.psfsDF.loc[0].copy
        else:
            self.refGraph = refDotGraph(self.refPSFs, brightThresh=brightThresh, maxError=maxXYMatchError)
            self.matchedRefDots = self.refPSFs.copy()
        # self.refGraph = refDotGraph(self.refPSFs)

        self.maxXYSearchError = maxXYSearchError
        self.maxZSearchError = maxZSearchError
        self.maxXYMatchError = maxXYMatchError
        self.maxZMatchError = maxZMatchError
        self.minMatches = minMatches
        self.stopRefProp = stopRefProp
        self.nLongestCheck = nLongestCheck
        '''
        def align_wrap(hyb):
            hybPSFs = self.psfsDF.loc[hyb]
            return self.refGraph.alignWith2(hybPSFs, maxXYError, maxZError, minMatches=minMatches,
                                            stopRefProp=stopRefProp, hyb_num=hyb, nLongestCheck=nLongestCheck)
        '''
        maxHyb = max(self.psfsDF.index)
        print('check if main')
        if __name__ == '__main__':
            pool = Pool(processes=5)
            print('Starting parallel Alignments:')
            alignment = pool.map(self.align_wrap, range(maxHyb))
            pool.close()
            pool.join()
        else:
            print('name:', __name__)

        offset_list = []
        self.matchesDF = pd.DataFrame(columns=('ref_row', 'ref_col', 'ref_z', 'comp_row', 'comp_col', 'comp_z',
                                               'aligned_row', 'aligned_col', 'aligned_z'))
        for hyb_results in alignment:
            offset_list.append(hyb_results[0])
            self.matchesDF.append(hyb_results[1])

        self.offsets = pd.DataFrame(data=offset_list, columns=('hyb', 'row_offset', 'col_offset', 'z_offset', 'row_offset_SE', 'col_offset_SE', 'z_offset_SE', 'n_offsets_calculated'))

        return self.offsets, self.matchesDF
        #else:
            #sprint('not main')

    def align_wrap(self, hyb):
        hybPSFs = self.psfsDF.loc[hyb]
        return self.refGraph.alignWith2(hybPSFs, self.maxXYSearchError, self.maxZSearchError, self.maxXYMatchError, self.maxZMatchError,
        minMatches=self.minMatches, stopRefProp=self.stopRefProp, hyb_num=hyb, nLongestCheck=self.nLongestCheck)


    def saveOffsets(self, filename):
        '''
        Save offsets to csv file
        :param filename: string name of file to save to
        :return:
        '''
        self.offsets.to_csv(filename)

    def loadOffsets(self, filename):
        '''
        Load saved offsets
        :param filename: string name of saved offsets file
        :return:
        '''
        self.offsets = pd.read_csv(filename, header=0)

    def saveMatchesDF(self, filename):
        self.matchesDF.to_csv(filename)

    def loadMatchesDF(self, filename):
        self.matchesDF = pd.read_csv(filename) #, index_col=('ref_row', 'ref_col', 'ref_z'))

    def compareMatchesToRef(self):

        plt.plot(self.matchesDF['comp_row'], self.matchesDF['comp_col'], 'xb', label='Pre-Aligned Dots')
        plt.plot(self.matchesDF['aligned_row'], self.matchesDF['aligned_col'], '+g', label='Alinged Dots')
        plt.plot(self.matchesDF['ref_row'], self.matchesDF['ref_col'], '.r', label='Reference Dots')
        plt.legend()
        plt.xlabel('Column (pixels)')
        plt.ylabel('Row (pixels)')
        plt.show()


        plt.plot(self.matchesDF['comp_row'], self.matchesDF['comp_z'], 'xb', label='Pre-Aligned Dots')
        plt.plot(self.matchesDF['aligned_row'], self.matchesDF['aligned_z'], '+g', label='Alinged Dots')
        plt.plot(self.matchesDF['ref_row'], self.matchesDF['ref_z'], '.r', label='Reference Dots')
        plt.legend()
        plt.xlabel('Row (pixels)')
        plt.ylabel('Z (slice)')
        plt.show()

        row_stdv = np.std(self.matchesDF['aligned_row'] - self.matchesDF['ref_row'])
        col_stdv = np.std(self.matchesDF['aligned_col'] - self.matchesDF['ref_col'])
        z_stdv = np.std(self.matchesDF['aligned_z'] - self.matchesDF['ref_z'])

        pa_row_stdv = np.std(self.matchesDF['comp_row'] - self.matchesDF['ref_row'])
        pa_col_stdv = np.std(self.matchesDF['comp_col'] - self.matchesDF['ref_col'])
        pa_z_stdv = np.std(self.matchesDF['comp_z'] - self.matchesDF['ref_z'])

        print("Pre-aligned Row Stdv:", pa_row_stdv, 'Pre-aligned Column Stdv:', pa_col_stdv, 'Pre-aligned z Stdv:', pa_z_stdv )

        print("Aligned Row Stdv:", row_stdv, 'Aligned Column Stdv:', col_stdv, 'Aligned z stdv:', z_stdv)



    def alignPSFs(self):
        '''
        Align PSFs in different hybs using the dot reference alignment.
        :return:
        alignedPSFs, data frame with a multiIndex: [hyb][psf_index]. Also saves to data member self.alignedPSFs
        '''
        row = []
        col = []
        row_offset_SE = []
        col_offset_SE = []
        hybn = []
        hyb_ind = []
        amp = []
        RMSE = []
        sigma = []
        for hyb, hybPSFs in enumerate(self.PSFs):
            for i, psf in enumerate(hybPSFs):
                amp.append(psf.amp)
                RMSE.append(psf.RMSE)
                sigma.append(psf.width)
                row.append(psf.row - self.offsets['row_offset'][hyb])
                # row.append(psf.row + self.offsets['row_offset'][hyb])
                col.append(psf.col - self.offsets['col_offset'][hyb])
                row_offset_SE.append(self.offsets['row_offset_SE'][hyb])
                col_offset_SE.append(self.offsets['col_offset_SE'][hyb])
                hybn.append(hyb)
                hyb_ind.append(i)
        index_arrays = [hybn, hyb_ind]
        index_tuples = list(zip(*index_arrays))
        index = pd.MultiIndex.from_tuples(index_tuples, names=['Hyb', 'Dot'])
        data = np.array([amp, RMSE, sigma, row, col, row_offset_SE, col_offset_SE]).transpose()
        column_names = ['amp', 'RMSE', 'sigma', 'row', 'col', 'row_offset_SE', 'col_offset_SE']
        self.alignedPSFs = pd.DataFrame(data=data, index=index, columns=column_names)
        return self.alignedPSFs

    def alignPSFsDF(self):
        self.alignedPSFs = pd.DataFrame()
        for hyb, hybPSFs in self.psfsDF.groupby(level=0):
            hybPSFs['row'] -= self.offsets['row_offset'][hyb]
            hybPSFs['col'] -= self.offsets['col_offset'][hyb]
            self.alignedPSFs = pd.concat((self.alignedPSFs, hybPSFs), axis=0)

    def convtPSFList2DF(self, PSFs=None):
        row = []
        col = []
        hybn = []
        hyb_ind = []
        amp = []
        RMSE = []
        sigma = []
        if PSFs is None:
            saveMember = True
            PSFs = self.PSFs
        else:
            saveMember = False

        for hyb, hybPSFs in enumerate(PSFs):
            for i, psf in enumerate(hybPSFs):
                amp.append(psf.amp)
                RMSE.append(psf.RMSE)
                sigma.append(psf.width)
                row.append(psf.row)
                col.append(psf.col)
                hybn.append(hyb)
                hyb_ind.append(i)
        index_arrays = [hybn, hyb_ind]
        index_tuples = list(zip(*index_arrays))
        index = pd.MultiIndex.from_tuples(index_tuples, names=['Hyb', 'Dot'])
        data = np.array([amp, RMSE, sigma, row, col]).transpose()
        column_names = ['amp', 'RMSE', 'sigma', 'row', 'col']
        psfsDF = pd.DataFrame(data=data, index=index, columns=column_names)
        if saveMember:
            self.psfsDF = psfsDF
        return psfsDF

    def saveAlignedPSFs(self, filename):
        '''
        Saves
        :param filename: String ame of save file.
        :return:
        '''
        self.alignedPSFs.to_csv(filename)

    def loadAlignedPSFs(self, filename):
        '''
        Load saved alined PSFs into self.alignedPSFs data member
        :param filename: string filename of saved PSFs file
        :return:
        '''
        self.alignedPSFs = pd.read_csv(filename, sep=',', header=0, index_col='Hyb')

    def compareHybAlignmentToRef(self, hyb, ref0=False):
        '''
        Plots the bead image with the bead peaks with aligned peaks from the given hyb number
        :param hyb: int number of the hyb to compare alignment to beads
        :return:
        '''
        # os.chdir('beads')
        if ref0:
            self.refPSFs = self.alignedPSFs.loc[0]
        beadim = skimage.io.imread(self.imfname)
        beadim = beadim[self.z][:, :, self.channel]

        plt.imshow(np.log2(beadim), cmap='Greys')

        # plot aligned
        plt.plot(self.alignedPSFs.loc[hyb]['col'], self.alignedPSFs.loc[hyb]['row'],
                 '+g', label='Aligned Fits from Hyb %d' % hyb)

        # plot fit peaks of bead image
        ref_rows = []
        ref_cols = []
        for i, psf in self.refPSFs.iterrows():
            ref_rows.append(psf['row'])
            ref_cols.append(psf['col'])
        # ref_rows = 2048 - np.array(ref_rows)
        plt.plot(ref_cols, ref_rows, 'xb', label='Reference Image Fit Centers')

        # self.alignedPSFs.loc[hyb].plot(x='row', y='col', )
        plt.legend()
        cbar = plt.colorbar()
        cbar.set_label('Log2 Intensity')
        # plt.plot()
        plt.show()

def align_par(hyb):
    return

class dotGraph(nx.Graph):
    '''
    Inherits from nx.Graph. ALso contains data and methods for sorting and matching.
    '''

    def __init__(self, PSFs=None, brightThresh=None, nBrightest=None):
        super(dotGraph, self).__init__()
        #self.edgesXSort = []
        if type(PSFs) != None:
            self.buildGraphPSFs(PSFs, brightThresh=brightThresh)

    def buildGraphPSFs(self, PSFs, brightThresh=None):
        '''
        Makes graph representing an image were PSFs are nodes nd edges containing associated displacement vectors
        connect all nodes/dots
        :param PSFs: dataframe of PSFs as returned by convtPSFLISt2DF. #list fitPSF objects
        :param brightThresh: If given, only consider PSFs brighter than this
        :return:
        '''

        # sPSFs = self.sortDotsBrightness(PSFs, brightThresh)
        # self.sPSFs = PSFs.sort_values(by=['amp'])
        
        if brightThresh:
            sPSFs = PSFs.sort_values(by='amp', ascending=False)
            sPSFs = sPSFs[sPSFs['amp'] > brightThresh]
        else:
            sPSFs = PSFs

        self.graph = nx.Graph()

        # print('Sorted PSFs, going on to iterative adding...')

        # Add add the dots to graph as a node
        # including its brightness and coordinates
        # all nodes are connected by edges holding the displacement vector between them.
        self.disp_vecs = []
        self.dispVecEdges = []
        self.dispVecRs = []
        # for psf in sPSFs:

        # for i, psf in self.sPSFs.iterrows():
        for i, psf in sPSFs.iterrows():
            # self.add_node((psf.row, psf.col), amp=psf.amp)
            if 'amp' in psf:
                self.add_node((psf['row'], psf['col'], psf['z']), amp=psf['amp'])
            else:
                self.add_node((psf['row'], psf['col'], psf['z']))
            # self.graph.add_node((psf['row'], psf['col']), amp=psf['amp'])
            # add edges between the new node and the already added nodes
            # for each node in the graph
            for node in self.__iter__():
                # if it is not the node that was just added
                # if node != (psf.row, psf.col):
                if node != (psf['row'], psf['col'], psf['z']):
                    # calculate the length of the node. First the euclidean length
                    # length = np.sqrt((node[0] - psf.row) ** 2 + (node[1] - psf.col) ** 2)
                    # then the row and column distances between the node in a smaller row and the node in a larger row
                    # to ensure the row component of the displacement vector is always positive
                    row_dist, col_dist, z_dist = self.findDispVec(node[0], node[1], node[2], psf['row'], psf['col'], psf['z'])

                    r = np.sqrt(row_dist**2 + col_dist**2)

                    # add an edge between the PSF added and the the node that was already in the graph
                    self.add_edge(node, (psf['row'], psf['col'], psf['z']), rd=row_dist, cd=col_dist, zd=z_dist, r=r, uamp=self.nodes[node]['amp'], vamp=psf['amp'])
                    self.disp_vecs.append((row_dist, col_dist, z_dist))
                    # add
                    if node[0] < psf['row']:
                        self.dispVecEdges.append(((node[0], node[1], node[2]), (psf['row'], psf['col'], psf['z'])))
                    else:
                        self.dispVecEdges.append(((psf['row'], psf['col'], psf['z']), (node[0], node[1], node[2])))
                    self.dispVecRs.append(r)


        self.disp_vecs = pd.DataFrame(self.disp_vecs, columns=('row', 'col','z'))#np.array(self.disp_vecs)

    def buildKDTree(self):
        #print('Building KDTree')
        self.edgeKDTree = latKDTree(self.disp_vecs)

    def sortEdgesRows(self):
        # make list of edges sorted by row distance
        edgeList = []
        for (u, v, eattrs) in self.edges.data():
            edgeList.append((u, v, eattrs))

        self.rsortedEdges = sorted(edgeList, key=lambda edge: edge[2]['rd'])

    def findDispVec(self, r1, c1, z1, r2, c2, z2):
        '''
        Find dispacement vector such that the row component is always positive.
        :param r1: row coordinate of first psf
        :param c1: column coordinate of first psf
        :param z1: z coordinate of first psf
        :param r2: row coordinate of second psf
        :param c2: column coordinate of second psf
        :param z2: z coordinate of second psf
        :return:
        row - the row component of the displacement vector between two psfs (always positive)
        col - the column component of the displacement vector between two psfs
        '''
        # if  < psf.row:
        if r1 < r2:
            row_dist = r2 - r1  # psf.row - node[0]
            col_dist = c2 - c1  # psf.col - node[1]
            z_dist = z2 - z1
        else:
            row_dist = r1 - r2  # node[0] - psf.row
            col_dist = c1 - c2  # node[1] - psf.col
            z_dist = z1 - z2
        return row_dist, col_dist, z_dist

    def sortDotsBrightness(self, PSFs, brightThresh=None):
        # sort the PSFs by brightness
        sPSFs = sorted(PSFs, key=lambda psf: psf.amp, reverse=True)

        if brightThresh:  # use dots brighter than threshold for alignment
            # binary search for thresh
            # initialize low index and high index. note that low indices hold
            # brighter dots and high indices hold dimmer dots
            lo = 0
            hi = len(sPSFs)
            while lo < hi:
                mid = (lo + hi) // 2
                # if a[mid] < x: lo = mid+1
                if sPSFs[mid].amp > brightThresh:
                    lo = mid + 1
                else:
                    hi = mid
            sPSFs = sPSFs[:lo]

        return sPSFs

    def sortNodesByEdges(self):
        '''
        Returns a list of Nodes sorted by their edges
        :return:
        sortedNodes - a list of nodes sorted by their edge displacement vector. Each element of sortedNodes
            is a 2 element list. The first element of the 2 element list is the node object.
            The second element is a tuple of tuples containing the node's edge displacement vectors. This tuple
            of displacement vectors is used for sorting and matching.
        '''
        self.sortedNodes = []
        toRemove = []
        for node, nbrs in self.adj.items():
            # make a sorted list of the edges from remaining nodes

            edgeList = []

            # if a node has fewer than three neighbors, mark for removal and continue to next node.
            if len(nbrs) < 3:
                toRemove.append(node)
            else:
                for nbr, eattrs in nbrs.items():
                    edgeList.append((eattrs['rd'], eattrs['cd']))
                edgeList = tuple(sorted(edgeList))  # , key=edgeValue))
                self.sortedNodes.append((node, edgeList))
        self.sortedNodes = sorted(self.sortedNodes, key=itemgetter(1))
        # self.graphs_Sorted_Nodes.append(nodeList)
        for node in toRemove:
            self.remove_node(node)
        return self.sortedNodes


class refDotGraph(dotGraph):
    '''
    Inherits from dotGraph. Holds the dot graph for the reference image for alignment.
    Contains methods for aligning with other images.
    '''

    def __init__(self, PSFs, brightThresh=None, nBrightest=None, maxError=0.1):
        print('Building Reference Graph')
        super(refDotGraph, self).__init__(PSFs, brightThresh, nBrightest)
        self.maxError = maxError

    def lookUpPosVec(self, rdist, cdist, zdist, maxXYError, maxZError):
        '''
        Look to see is a position vector separating dots in the reference image matches a given vector
        :param rdist: the row component of position vector under consideration
        :param cdist: the column component of position vector under consideration
        :param zdist: the z component of position vector under consideration
        :return:
        edges - list of matching edges denoted as 2 element tuple of tuples of coordinates of nodes at end of edge
        disp_vec_matches - list of position vectors corresponding to each matching edge
        '''
        matches, matchInds = self.edgeKDTree.query_cylinder_point((rdist, cdist, zdist), maxXYError, maxZError)
        disp_vec_matches = matches#self.disp_vecs[matches]
        edges = []
        for matchInd in matchInds:
            edges.append(self.dispVecEdges[matchInd])
        return edges, disp_vec_matches

    def alignWith(self, dGraph, maxXYError, maxZError = 2, minMatches=3, hyb_num=0):
        '''
        Finds alignment with another dot Graph
        :param dGraph:
        :return:
        '''
        print('Finding intersection')
        matchDict = self.findEdgeIntersection(dGraph, maxXYError=maxXYError, maxZError=maxZError)
        offset = self.estOffset(matchDict)
        matchesDF = self.getMatchDF(matchDict, offset, minMatches, hyb_num)
        # refMatchGraph, compMatchGraph = self.buildMatchedGraphs(matches)
        # offset = self.estimateOffset(refMatchGraph, compMatchGraph, maxError)
        # return offset, matches, refMatchGraph, compMatchGraph, matchDict
        return offset, matchesDF

    def alignWith2(self, compPSFs, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError, minMatches=3,
                   brightThresh=0, stopRefProp=0.5, hyb_num=0, nLongestCheck=1000):
        '''
        Reference dots with dots from a comparison image
        :param compPSFs: DataFrame of dots for a hybridization formatatted like output of convtPSFList2DF
        :param brightThresh: Ignore dots dimmer than this
        :return:
        Offset of hyb from reference
        '''
        self.buildKDTree()
        print('Starting Alignment: Searching for matching edge in reference graph')
        if brightThresh:
            compPSFs = compPSFs[compPSFs['amp'] > brightThresh]
        #matchDict = self.findMatches(compPSFs, maxXYError, maxZError, minMatches, stopRefProp=stopRefProp, nBrightestCheck=nBrightestCheck)
        matchDict = self.findMatchesLongestEdges(compPSFs, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError, minMatches, stopRefProp=stopRefProp,
                                                 nLongestCheck=nLongestCheck)
        offset = self.estOffset(matchDict, minMatches=minMatches)
        matchesDF = self.getMatchDF(matchDict, offset, minMatches, hyb_num)
        return offset, matchesDF

    def findEdgeIntersection(self, dGraph, maxXYError, maxZError, n=None):
        '''
        Looks for matching edges in reference and comparison graphs. Keeps track of nodes that are matched in matching
        edges in the matchDict dictionary. See return for description of matchDict.
        :param dGraph: another dot graph with which to find the intersection
        :param maxXYError: the maximum allowed lateral error between vectors
        :param maxZError: the maximum allowed z error between vectors
        :param n: the number of the dot graph comparing to
        :return:
        dictionary of matched nodes. Keys are tuples containing first a tuple of the coordinates of the dot in the
        reference graph, then a tuple of the coordinents in the comparison dot. The value is the number of edges
        each pair was matched in.
        '''

        ref_rsortedEdges = copy.copy(self.rsortedEdges)
        comp_rsortedEdges = copy.copy(dGraph.rsortedEdges)

        maxXYErrorSq = maxXYError ** 2

        matches = []
        matchDict = {}
        print('n_refEdges:', len(ref_rsortedEdges), 'n_compEdges:', len(comp_rsortedEdges))

        # for each edge in dGraph
        for edge in comp_rsortedEdges:
            # look for matches within error
            redge = 0
            edge_ref_matches = []
            # print('Matching Edge:', edge)
            while redge < len(ref_rsortedEdges):
                # if the reference edge is mre than the maxError rows lower than the edge under consideration
                if ref_rsortedEdges[redge][2]['rd'] < edge[2]['rd'] - maxXYError:
                    # remove the reference edge from consideration
                    ref_rsortedEdges.pop(0)

                # if the reference edge is more than the maxError rows higher than the edge under consideration
                elif ref_rsortedEdges[redge][2]['rd'] > edge[2]['rd'] + maxXYError:
                    # no more matches will be found, move onto the next edge
                    break

                # if the reference edge is within error of the edge under consideration
                elif (ref_rsortedEdges[redge][2]['rd'] - edge[2]['rd']) ** 2 + (ref_rsortedEdges[redge][2]['cd'] - edge[2]['cd']) ** 2 < maxXYErrorSq\
                        and ref_rsortedEdges[redge][2]['zd'] - edge[2]['zd'] < maxZError:
                    # add to matches
                    edge_ref_matches.append(ref_rsortedEdges[redge])
                redge += 1

            if len(edge_ref_matches) > 1:
                # There is more than one match! Don't bother considering it.
                continue

            elif len(edge_ref_matches) == 1:
                matches.append([edge_ref_matches[0], edge])
                if (edge_ref_matches[0][0], edge[0]) in matchDict:
                    matchDict[(edge_ref_matches[0][0], edge[0])] += 1
                else:
                    matchDict[(edge_ref_matches[0][0], edge[0])] = 1

            # otherwise, no matches

        return matchDict

    def findMatches(self, compPSFs, maxXYError, maxZError, minMatches, stopRefProp = 0.5, nBrightestCheck=1000):
        '''
        Find dots from an image that match those in the reference image. To do this, first find a position
        vector between dots in the image that matches a position vector in the reference image to match the dots
        on either end of the the position vector. Then, look for matches to the known neighbors of the dots in the reference
        image by following the displacement vectors stored in its edges.
        :param compPSFs:
        :return:
        '''
        matchDict = {}
        edgeMatchDict = {}
        # sPSFs = self.sortDotsBrightness(PSFs, brightThresh)
        sCompPSFs = compPSFs.sort_values(by='amp', ascending=False)
        compPSFsTree = latKDTree(sCompPSFs.loc[:, 'row':'z'])
        #self.compPSFsTree = KDTree(sCompPSFs.loc[:, 'row':'col'])#:'z'])
        nPSFs = len(sCompPSFs)
        found_enough = False
        # First, look for a position vector that matches a position vector in the reference image
        print('Searching for matching edge vectors')
        for i in range(nPSFs):
            for j in range(i):
                iPSF = sCompPSFs.iloc[i]
                jPSF = sCompPSFs.iloc[j]
                rdist, cdist, zdist= self.findDispVec(iPSF['row'], iPSF['col'], iPSF['z'], jPSF['row'], jPSF['col'], jPSF['z'])
                # check whether there is a matching position in the reference image
                # match = self.edgeKDTree.query_ball_point((rdist, cdist), self.maxError)
                edgeMatch, dispVecMatches = self.lookUpPosVec(rdist, cdist, zdist, maxXYError, maxZError)
                if len(edgeMatch) == 1 and edgeMatch[0] not in edgeMatchDict:
                    # Start graph traversal which adds matches to matchDict
                    edgeMatchDict[edgeMatch[0]] = 1
                    self.travRefFindMatches(compPSFsTree, edgeMatch, iPSF, jPSF, matchDict, edgeMatchDict, maxXYError, maxZError, minMatches, stopRefProp)

                    if self._find_n_well_matched(matchDict, minMatches) > self.number_of_nodes() * stopRefProp:
                        found_enough = True
                        break

                    if i > nBrightestCheck:
                        # time to give up
                        print('giving up search')
                        found_enough = True
                        break
            if found_enough:
                break
        return matchDict

    def findMatchesLongestEdges(self, compPSFs, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError,
                                minMatches, nLongestCheck=100, stopRefProp = 0.5):
        '''
        Search dots in the readout graph for matches in the reference graph by trying to find dots connected
        by the longest edges in the reference graph. Once an edge in the reference graph is found, traverses the
        reference graph to match the rest of the nodes.
        :param compPSFs: DataFrame of dots in the readout image. Must have columns
        :param maxXYMatchError: Maximum allowed lateral error between matched position vectors.
        :param maxZMatchError: Maximum allowed z error between matched position vectors
        :param minMatches: Minimum number of edges in which a reference readout dot must be matched for the pair
        to be included in the offset calculation
        :param nLongestCheck: Number of longest edges in the reference  graph to search before giving up.
        :param stopRefProp: the proportion of nodes in the reference graph, that when matched with at least minMatches
        edges to dots in the readout image, to consider enough to stop searching for more.
        :return:
        '''
        self.matchDict = {}
        self.edgeMatchDict = {}
        self.compPSFsTree = latKDTree(compPSFs.loc[:, 'row':'z'])

        edgeList = []
        for (u, v, eattrs) in self.edges.data():
            edge = []
            for obj in (u, v, eattrs):
                for entry in obj:
                    if type(obj) == tuple:
                        edge.append(entry)
                    else:
                        edge.append(obj[entry])

            edgeList.append(edge)
        edgeDF = pd.DataFrame(edgeList, columns=('u_row', 'u_col', 'u_z', 'v_row', 'v_col', 'v_z', 'row_dist', 'col_dist', 'z_dist', 'r', 'uamp', 'vamp'))
        edgeDF.sort_values(by='r', ascending=False, inplace=True)

        j = 1
        for i, edge in edgeDF.iterrows():
            # search for corresponding edge in readout image
            self._find_ref_edge(edge, compPSFs, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError, minMatches, stopRefProp)
            if self._find_n_well_matched(self.matchDict, minMatches) > self.number_of_nodes() * stopRefProp:
                break

            if j > nLongestCheck:
                # time to give up
                print('giving up search')
                break
            else:
                j += 1

        return self.matchDict

    def _find_ref_edge(self, edge, compPSFs, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError,
                       minMatches, stopRefProp):
        '''
        Given an edge in the reference graph, look for the corresponding matching edge in the dataframe
        of dots in the readout image.
        :param edge: Pandas series with data representing an edge: 'u_row', 'u_col', 'u_z', 'v_row', 'v_col', 'v_z', 'row_dist',
        'col_dist', 'z_dist', 'r', 'uamp', 'vamp'
        :param compPSFs: DataFrame of dots in the readout image. Must have columns
        :param maxXYMatchError: Maximum allowed lateral error between matched position vectors.
        :param maxZMatchError: Maximum allowed z error between matched position vectors
        :param minMatches: Minimum number of edges in which a reference readout dot must be matched for the pair
        :param stopRefProp: the proportion of nodes in the reference graph, that when matched with at least minMatches 
        edges to dots in the readout image, to consider enough to stop searching for more.                             
        :return:
        '''
        bright_thresh_factor = 2

        maxXYErrorsq = maxXYSearchError ** 2
        
        #subset points in read out image to search through
        # u is the dot with lower row, v is the dot with higher row
        urmax = 2048 - edge['row_dist'] + maxXYSearchError
        vrmin = edge['row_dist'] - maxXYSearchError

        u_cands = compPSFs.iloc[list(compPSFs['row'] < urmax)]
        v_cands = compPSFs.iloc[list(compPSFs['row'] > vrmin)]

        if edge['col_dist'] > 0:
            ucmax = 2048 - edge['col_dist'] + maxXYSearchError
            vcmin = edge['col_dist'] - maxXYSearchError

            u_cands = u_cands.iloc[list(u_cands['col'] < ucmax)]
            v_cands = v_cands.iloc[list(v_cands['col'] > vcmin)]
        elif edge['col_dist'] < 0:
            ucmin = - edge['col_dist'] - maxXYSearchError
            vcmax = 2048 + edge['col_dist'] + maxXYSearchError

            u_cands = u_cands.iloc[list(u_cands['col'] > ucmin)]
            v_cands = v_cands.iloc[list(v_cands['col'] < vcmax)]

        # the u dot from the edge list is the note that was added first to the reference graph, here u is the dot of less row and v is the dot of higher row

        if edge['u_row'] < edge['v_row']:
            # if the edge list u dot is the same as the u dot here
            ucands_amp_gt_lbnd = np.greater(np.array(u_cands['amp']), float(edge['uamp'] / bright_thresh_factor))
            ucands_amp_lt_ubnd = np.less(np.array(u_cands['amp']), float(edge['uamp']) * bright_thresh_factor)
            u_cands = u_cands.loc[list(np.logical_and(ucands_amp_gt_lbnd, ucands_amp_lt_ubnd))]

            v_cands_amp_gt_lbnd = np.greater(np.array(v_cands['amp']), float(edge['vamp'])/bright_thresh_factor)
            v_cands_amp_lt_ubnd = np.less(np.array(v_cands['amp']), float(edge['vamp'])*bright_thresh_factor)
            v_cands = v_cands.loc[list(np.logical_and(v_cands_amp_gt_lbnd, v_cands_amp_lt_ubnd))]
            
            uv_inverted = False

        elif edge['u_row'] > edge['v_row']:
            # if the edgelist u dot is the v dot here
            u_cands_gt_lbnd = np.greater(np.array(u_cands['amp']), float(edge['vamp'])/bright_thresh_factor)
            u_cands_lt_ubnd = np.less(np.array(u_cands['amp']), float(edge['vamp'])*bright_thresh_factor)
            u_cands = u_cands.loc[list(np.logical_and(u_cands_gt_lbnd, u_cands_lt_ubnd))]

            v_cands_amp_gt_lbnd = np.greater(np.array(v_cands['amp']), float(edge['uamp'])/bright_thresh_factor)
            v_cands_amp_lt_ubnd = np.less(np.array(v_cands['amp']), float(edge['uamp'])*bright_thresh_factor)
            v_cands = v_cands.loc[list(np.logical_and(v_cands_amp_gt_lbnd, v_cands_amp_lt_ubnd))]
            
            uv_inverted = True


        u_cands.sort_values(by='row', inplace=True)
        v_cands.sort_values(by='row', inplace=True)

        minV = 0
        #print('n ucands:', len(u_cands))
        #print('n vcands:', len(v_cands))

        for i, udot in u_cands.iterrows():
            #print(minV, 'of', len(v_cands))
            for j in range(minV, len(v_cands)):
                if not uv_inverted:
                    rdist = v_cands.iloc[j]['row'] - udot['row']
                    cdist = v_cands.iloc[j]['col'] - udot['col']
                    zdist = v_cands.iloc[j]['z'] - udot['z']
                elif uv_inverted:
                    rdist = udot['row'] - v_cands.iloc[j]['row'] 
                    cdist = udot['col'] - v_cands.iloc[j]['col']
                    zdist = udot['z'] - v_cands.iloc[j]['z']
                    
                if rdist < edge['row_dist'] - maxXYSearchError:
                    minV += 1
                    #print('Increment minV to', minV)
                    
                elif rdist > edge['row_dist'] + maxXYSearchError:
                    #print('rdist to vcand out of bounds. Break')
                    break
                
                elif (rdist- edge['row_dist'])**2 + (cdist - edge['col_dist'])**2 < maxXYErrorsq and np.abs(zdist - edge['z_dist']) < maxZSearchError:
                    # found match! run graph traversal
                    match = [(tuple(edge['u_row':'u_z']), tuple(edge['v_row':'v_z']))]
                    self.travRefFindMatches(self.compPSFsTree, match, udot, v_cands.iloc[j], self.matchDict, self.edgeMatchDict,
                                            maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError, minMatches, stopRefProp)


    def travRefFindMatches(self, compPSFsTree, match, cPSF1, cPSF2, matchDict, edgeMatchDict, maxXYSearchError,
                           maxZSearchError, maxXYMatchError, maxZMatchError, minMatches, stopRefProp):
        '''
        Given two dots in the readout image matched to dots in the reference graph by an edge, traverse the reference graph to find
        the rest of the fiducial markers in the readout image.
        :param compPSFsTree: Search tree
        :param match: 1 element list containing edge tuple from reference graph (that contains tuples of coordinates of each node)
            example: [((1,1,2), (4,4,2))]
        :param cPSF1: Series of information describing the first dot in the readout image matching the reference edge
        :param cPSF2: Series of information describing the second dot in the readout image matching the reference edge
        :param matchDict: Dictionary of matches between reference and comparison graph. Keys are two element tuples
            where the first element is a tuple of the coordinates fo the reference dot, and the second element is a tuple
            of the coordinates of the comparison dot. The value is tne number of times the dots were matched by position
            vectors.
        :param edgeMatchDict:
        :param maxXYSearchError:
        :param maxZSearchError:
        :param maxXYMatchError:
        :param maxZMatchError: Maximum allowed z error between matched position vectors
        :param minMatches: Minimum number of edges in which a reference readout dot must be matched for the pair
             to be included in the offset calculation
        :param stopRefProp:   the proportion of nodes in the reference graph, that when matched with at least minMatches
             edges to dots in the readout image, to consider enough to stop searching for more.
        :return:
        '''
        print('Starting Ref Graph Traversal')

        # add 2 tuples where the first element is a dot in the reference image, and the second is a dot matched to it
        # in the comparison image
        dotQueue = []
        unSearchedEdgeQueue = []
        nAmbiguous = 0
        nMatchedAll = 0

        # add first matching PSF to match dictionary
        # the psf with smaller row value is matched to the first node in the edge
        if cPSF1['row'] < cPSF2['row']:
            match1Key = (match[0][0], (cPSF1['row'], cPSF1['col'], cPSF1['z']))
            match2Key = (match[0][1], (cPSF2['row'], cPSF2['col'], cPSF2['z']))
            #dotQueue.append(match1Key)
            #dotQueue.append(match2Key)
        else:
            match1Key = (match[0][0], (cPSF2['row'], cPSF2['col'], cPSF2['z']))
            match2Key = (match[0][1], (cPSF1['row'], cPSF1['col'], cPSF1['z']))
            #dotQueue.append(match1Key)
            #dotQueue.append(match2Key)
        self._placeMatchInDictQueue(match1Key, match2Key, matchDict, dotQueue, maxXYMatchError, maxZMatchError)
        self._placeMatchInDictQueue(match2Key, match1Key, matchDict, dotQueue, maxXYMatchError, maxZMatchError)

        # Now traverse the graph
        while len(dotQueue) > 0:
            #print(len(dotQueue))
            # first element of matched dot is reference coordinates, second element is comparison image coordinates
            matchedDot = dotQueue.pop(0)
            neighbor_list = list(self[matchedDot[0]].copy())

            nUnmatched = 0
            nMatched = 0
            # look for each of the dots known neighbors in the reference graph in the comparison image
            #for refNeighbor in self[matchedDot[0]]:
            while neighbor_list:
                refNeighbor = neighbor_list.pop(0)
                # if edge not already matched:
                if (matchedDot[0], refNeighbor) not in edgeMatchDict and (refNeighbor, matchedDot[0]) not in edgeMatchDict:
                    matched, ambiguous = self._findEdgeInReadOut(matchedDot, refNeighbor, compPSFsTree, matchDict,
                                        dotQueue, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError)
                    nAmbiguous += ambiguous
                    nMatched += matched
                    nUnmatched += (not matched and not ambiguous)
                if nUnmatched > 10 and nMatched == 0:
                    # give up, this isn't a real match
                    break
                if matchedDot in matchDict and matchDict[matchedDot] >= minMatches:
                    # found enough matches, add un searched edges to queue in case if not enough nodes are found in the traversal without them.
                    unSearchedEdgeQueue.append((matchedDot, neighbor_list))
                    break
            nMatchedAll += nMatched

        n_well_matched = self._find_n_well_matched(matchDict, minMatches)
        if n_well_matched < self.number_of_nodes() * stopRefProp:
            print('Looking through reserved edges.')
            # look through unSearched edges
            while unSearchedEdgeQueue and self._find_n_well_matched(matchDict, minMatches) < self.number_of_nodes() * stopRefProp:
                matchedDot, neighbor_list = unSearchedEdgeQueue.pop(0)
                for refNeighbor in neighbor_list:
                    if (matchedDot[0], refNeighbor) not in edgeMatchDict and (refNeighbor, matchedDot[0]) not in edgeMatchDict:
                        matched, ambiguous = self._findEdgeInReadOut(matchedDot, refNeighbor, compPSFsTree, matchDict, dotQueue, maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError)
                        nAmbiguous += ambiguous

        if nMatchedAll > 0:
            #self._placeMatchInDictQueue(match1Key, match2Key, matchDict, dotQueue, maxXYMatchError, maxZMatchError)
            #self._placeMatchInDictQueue(match2Key, match1Key, matchDict, dotQueue, maxXYMatchError, maxZMatchError)
            print('Number of matched Dots:', self._find_n_well_matched(matchDict, minMatches), 'of %d in reference.' % self.number_of_nodes(), 'In traversal,', nMatchedAll,
                  'edges unambiguously matched,', nAmbiguous, 'ambiguous edge matches.')
            if self._find_n_well_matched(matchDict, minMatches) > self.number_of_nodes():
                raise Exception('Found more more matches than reference dots.')
        else:
            self.matchDict.pop(match1Key)
            self.matchDict.pop(match2Key)
            print("No matches found.")

    def _findEdgeInReadOut(self, matchedDot, refNeighbor, compPSFsTree, matchDict, dotQueue, maxXYSearchError,
                           maxZSearchError, maxXYMatchError, maxZMatchError):
        '''
        Given a dot matched in the reference graph and readout image, searches for a known neighbor of the reference dot
        in the readout image.
        :param matchedDot:
        :param refNeighbor:
        :param compPSFsTree: latKDTree of dots in the readout image.
        :param matchDict: Dictionary of matches between reference and comparison graph. Keys are two element tuples
            where the first element is a tuple of the coordinates fo the reference dot, and the second element is a tuple
            of the coordinates of the comparison dot. The value is tne number of times the dots were matched by position
            vectors.
        :param dotQueue: list (used as queue) of matched dots for which to search for additional edges.
        :param maxXYMatchError: Maximum allowed lateral error between position vectors when matching
        :param maxZMatchError: Maximum allowed z error between position vectors when matching.
        :return:
        Boolean matched
        Boolean ambiguous - true if more than one match
        '''
        r_to_neighbor = refNeighbor[0] - matchedDot[0][0]
        c_to_neighbor = refNeighbor[1] - matchedDot[0][1]
        z_to_neighbor = refNeighbor[2] - matchedDot[0][2]
        matches, matchInds = compPSFsTree.query_cylinder_point(
            (matchedDot[1][0] + r_to_neighbor, matchedDot[1][1] + c_to_neighbor, matchedDot[1][2] + z_to_neighbor),
            maxXYSearchError, maxZSearchError)
        # check if only one match, ignore if otherwise (ambiguous)

        # print('nMatches:', len(matches))
        if len(matches) == 1:
            newMatchedPSF = matches  # compPSFs.iloc[matches]
            # print('matchedDot:', matchedDot, ', refNeighbor:', refNeighbor)
            # print('newMatchedPSF:', newMatchedPSF)
            # toDo: 
            newMatchedKey = (
                refNeighbor, (float(newMatchedPSF['row']), float(newMatchedPSF['col']), float(newMatchedPSF['z'])))
            self._placeMatchInDictQueue(newMatchedKey, matchedDot, matchDict, dotQueue, maxXYMatchError, maxZMatchError)
            matched = True
            ambiguous = False
        elif len(matches) == 0:
            matched = False
            ambiguous = False
        elif len(matches) > 1:
            ambiguous = True
            matched = False
        return matched, ambiguous

    def _find_n_well_matched(self, matchDict, minMatches):
        '''
        Returns the number of dots that have been matched by the minimum number of edges to allow inclusion
        in offset calculation.
        :param matchDict: Dictionary with matched dot pairs between reference and readout image as keys, and number of
        edges matching them as values.
        :param minMatches: The minimum allowed number of matches for a reference-readout dot matching pair to be included
        in the offset calculation
        :return:
        The number of reference-readout dot match pairs that have  passed the threshold to be included in the
        offset calculation.
        '''
        # Find number of well matched dots
        n_well_matched = len([match for match in matchDict if matchDict[match] > minMatches])
        return n_well_matched

    def _placeMatchInDictQueue(self, newMatch, oldMatch, matchDict, dotQueue, maxXYMatchError, maxZMatchError):
        '''
        Place matched psfs from the reference graph and comparison image in the match dictionary and Queue to look for neighbors.
        The match dictionary keeps track of how many each pair of dots in the reference and comparision image have been matched by edges.
        :param newMatch: two element tuple where the first element is a three element tuple of the coordinates (row, col, z) of the dot in the
            reference image, and the second element is a three element tuple of the dot in the comparison image' coodinates.
            Example: ((1,1,1), (2,2,1); Example 2: ((row1, col1, z1), (row2, col2, z2))
        :param matchDict: Dictionary using tuples of matched dots as described above as keys, that keeps track of the number of
            times such pairs are matched.
        :param dotQueue: Queue in which to place new matches to look for more matching neighbors.
        :return:
        '''
        if newMatch in matchDict:
            if matchDict[newMatch] == -1:
                if self.withinError(newMatch, oldMatch, maxXYMatchError, maxZMatchError):
                    #print('now within error')
                    matchDict[newMatch] = 1
                    if oldMatch in matchDict:
                        matchDict[oldMatch] += 1
                    else:
                        matchDict[oldMatch] = 1
                #else:
                    #print('still not within error')
            else:
                matchDict[newMatch] += 1
                if oldMatch in matchDict:
                    matchDict[oldMatch] += 1
                else:
                    matchDict[oldMatch] = 1
        else:
            # print("Adding to dotQueue:", match)
            if self.withinError(newMatch, oldMatch, maxXYMatchError, maxZMatchError):
                matchDict[newMatch] = 1
                # ToDo: figure out more rigorous solution to this issue that exception catching
                try:
                    matchDict[oldMatch] += 1
                except KeyError:
                    matchDict[oldMatch] = 1
            else:
                matchDict[newMatch] = -1
            dotQueue.append(newMatch)

    def withinError(self, newMatch, oldMatch, maxXYError, maxZError):
        #XYErrorSq = (newMatch[1][0] - oldMatch[1][0]) ** 2 + (newMatch[1][1] - oldMatch[1][1]) ** 2
        #ZErrorAbs = np.abs(newMatch[1][2] - oldMatch[1][2])
        #print('XYErrorSq:', XYErrorSq, 'ZErrorAbs:', ZErrorAbs)
        #withinError = XYErrorSq < maxXYError**2 and ZErrorAbs < maxZError
        roRowDiff = newMatch[1][0] - oldMatch[1][0]
        RefRowDiff = newMatch[0][0] - oldMatch[0][0]

        roColDiff = newMatch[1][1] - oldMatch[1][1]
        RefColDiff = newMatch[0][1] - oldMatch[0][1]

        roZDiff = newMatch[1][2] - oldMatch[1][2]
        RefZDiff = newMatch[0][2] - oldMatch[0][2]
        XYErrorSq = (roRowDiff - RefRowDiff)**2 + (roColDiff - RefColDiff)**2
        ZErrorAbs = np.abs(roZDiff - RefZDiff)
        withinError = XYErrorSq < maxXYError ** 2 and ZErrorAbs < maxZError
        return withinError

    def estOffset(self, matchDict, minMatches = 3):
        '''
        Given a dictonary of matching nodes between the reference and comparison graph, estimates the offset
        between the two images.
        :param matchDict: dictionary of matches between reference and comparison graph. Keys are two element tuples
            where the first element is a tuple of the coordinates fo the reference dot, and the second element is a tuple
            of the coordinates of the comparison dot. The value is tne number of times the dots were matched by position
            vectors.
        :param minMatches: the minimum number of times dots must have been matched to be included in offset calculation
        :return:
        r_mean - mean row offset between matches
        c_mean- mean column offset between matches
        r_mean_se - standard error of row offset between matches
        c_mean_se - standard error of column offset between matches
        r_weight_mean - weighted mean (weighted by number of matches) of row offsets
        c_weight_mean - weighted mean (weighted by number of matches) of column offsets
        '''
        #Gather matches of coordinates into arrays
        rdisp = []
        cdisp = []
        zdisp = []
        nmatches = []
        for match in matchDict:
            # only consider dots that were matched more than the minimum allowed number of times
            if matchDict[match] >= minMatches:
                rdisp.append(match[1][0] - match[0][0])
                cdisp.append(match[1][1] - match[0][1])
                zdisp.append(match[1][2] - match[0][2])
                nmatches.append(matchDict[match])

        dot_r_offsets = np.array(rdisp)
        dot_c_offsets = np.array(cdisp)
        dot_z_offsets = np.array(zdisp)
        n_dot_edges = np.array(nmatches)

        # remove any outliers in XY that my have been introduced by mismatches. Don't consider z since z measurements are poor.
        outlier = True
        while outlier:
            # compute summary stats
            r_mean = np.mean(dot_r_offsets)
            c_mean = np.mean(dot_c_offsets)
            rstdv = np.std(dot_r_offsets)
            cstdv = np.std(dot_c_offsets)

            # check if there are lateral outliers
            outliers = np.logical_or(abs(dot_r_offsets - r_mean) > 3 * rstdv, abs(dot_c_offsets - c_mean) > 3 * cstdv)
            outlier = np.any(outliers)

            # remove outliers
            dot_r_offsets = dot_r_offsets[np.logical_not(outliers)]
            dot_c_offsets = dot_c_offsets[np.logical_not(outliers)]
            n_dot_edges = n_dot_edges[np.logical_not(outliers)]

        z_mean = np.mean(dot_z_offsets)
        zstdv = np.std(dot_z_offsets)

        r_mean_se = rstdv / np.sqrt(len(dot_r_offsets))
        c_mean_se = cstdv / np.sqrt(len(dot_c_offsets))
        z_mean_se = zstdv / np.sqrt(len(dot_z_offsets))

        r_weight_mean = np.sum(dot_r_offsets * n_dot_edges) / np.sum(n_dot_edges)
        c_weight_mean = np.sum(dot_c_offsets * n_dot_edges) / np.sum(n_dot_edges)

        return r_mean, c_mean, z_mean, r_mean_se, c_mean_se, z_mean_se, r_weight_mean, c_weight_mean, len(dot_r_offsets)

    def getMatchDF(self, matches, offsets, minMatches, hyb_num):
        '''
        returns dataframe of fidicuial coordinates in reference image, matched coordinates in hyb images,
        and aligned coordinates off in hyb images. Columns indexed by hyb and then number in hyb
        :param matches: dictionary of matches
        :param offsets: the output of estOffset (see above)
        :return:
        Dataframe as described above
        '''
        r_offset, c_offset, z_offset, r_mean_se, c_mean_se, z_mean_se, r_weight_mean, c_weight_mean, n_dots = offsets
        ref_row = []
        ref_col = []
        ref_z = []
        comp_row = []
        comp_col = []
        comp_z = []
        for match in matches:
            if matches[match] > minMatches:
                ref_row.append(match[0][0])
                ref_col.append(match[0][1])
                ref_z.append(match[0][2])
                comp_row.append(match[1][0])
                comp_col.append(match[1][1])
                comp_z.append(match[1][2])

        #index_arrays = [list(range(len(ref_row))), [hyb_num]*len(ref_row)]
        #index_tuples = list(zip(*index_arrays))
        #index = pd.MultiIndex.from_tuples(index_tuples, names=['Hyb', 'Dot'])
        matchesDF = pd.DataFrame()#index=index)
        matchesDF['Hyb'] = [hyb_num]*len(ref_row)
        matchesDF['ref_row'] = ref_row
        matchesDF['ref_col'] = ref_col
        matchesDF['ref_z'] = ref_z
        matchesDF['comp_row'] = comp_row
        matchesDF['comp_col'] = comp_col
        matchesDF['comp_z'] = comp_z
        matchesDF['aligned_row'] = matchesDF['comp_row'] - r_offset
        matchesDF['aligned_col'] = matchesDF['comp_col'] - c_offset
        matchesDF['aligned_z'] = matchesDF['comp_z'] - z_offset

        return matchesDF

class latKDTree(KDTree):
    '''
    Allows look up of vectors by the lateral components using a KDTree, filters Z components later
    '''
    def __init__(self, data):
        '''
        Builds KDTree from lateral components of input vectors
        :param data: Dataframe containing information describing 3D position vectors. Must mave two consecutive columns: "row" and "col"
        must have additional z column.
        '''
        self.input_data = data
        super(latKDTree, self).__init__(self.input_data.loc[:, 'row':'col'])


    def query_cylinder_point(self, v, errorXYMax, errorZMax):
        '''
        Searches component vectors using first a KDTree for the lateral components, then checks whether match's z
        component is within error.
        :param v: Three dimensional vector to search for
        :param errorXYMax: maximum allowed lateral error
        :param errorZMax: maximum allowed z error
        :return:
        matches - subset dataframe including rows of position vectors matching v
        acceptInds - row indices of position vector dataframe  containing position vectors matching v
        '''
        v_xy = v[:2]
        v_z = v[2]
        latMatchesInds = self.query_ball_point(v_xy, errorXYMax)
        #latMatches = self.input_data.iloc[latMatchesInds, :]
        acceptInds = []
        #for i, match in latMatches.iterrows():
        for i in latMatchesInds:
            if np.abs(v_z - self.input_data.iloc[i]['z']) < errorZMax:
                acceptInds.append(i)
        matches = self.input_data.iloc[acceptInds, :]
        return matches, acceptInds



class spatMatDotHolder:
    '''
    maintains a sorted list of barcodes, sorted by their taxi cab distance from (0,0)
    '''

    def __init__(self):
        self.spatMatDots = []
        self.keys = []
        self.curr_hyb = 0

    def set_current_hyb(self, hyb):
        self.curr_hyb = hyb

    def insert_new_spatDot(self, spatDot):
        # insert spatial dot in order by taxiCabDist_ave
        pos = bisect(self.keys, spatDot.taxiCabDist_ave)
        self.keys.insert(pos, spatDot.taxiCabDist_ave)
        self.spatMatDots.insert(pos, spatDot)

    def matchPSF(self, psf):
        # do binary search for PSF with closest taxicab distance
        closest = bisect(self.keys, psf['row'])

        # check whether PSFs within sqrt(2) maxError match. If more than 1 matches, set aside/toss
        # first look for spatial dots with centers with longer taxicab distances than the psf's
        pos = closest
        found = False
        while pos < len(self.spatMatDots) and np.abs(
                self.spatMatDots[pos].taxiCabDist_ave - psf['row']) < spatMatDot.maxError * np.sqrt(2):
            if self.spatMatDots[pos].checkMatch(psf):
                self._addPSFToSpatDot(pos, psf)
                found = True
                break
            pos += 1

        # now check for barcodes with centers with shorter taxicab distances than the psf's
        if not found:
            pos = closest - 1
            while pos >= 0 and np.abs(
                    self.spatMatDots[pos].taxiCabDist_ave - psf['row']) < spatMatDot.maxError * np.sqrt(2):
                if self.spatMatDots[pos].checkMatch(psf):
                    self._addPSFToSpatDot(pos, psf)
                    found = True
                    break
                pos -= 1
        return found

    def _addPSFToSpatDot(self, pos, psf):
        # first add the psf to the barcode
        self.spatMatDots[pos].addHybPSF(self.curr_hyb, psf)

        # now check that the barcode is still sorted properly
        new_row_ave = self.spatMatDots[pos].row_ave

        # if new taxicab distance is now larger than what used to be the next largest
        if pos != len(self.keys) - 1 and self.keys[pos + 1] < new_row_ave:
            print('larger')
            # find how many positions this barcode should move up in the list
            n = 2
            while self.keys[pos + n] < new_row_ave and pos + n < len(self.spatMatDots):
                n += 1
            # remove and reinsert spatial dot and key into their lists (inserts to left of insertion position)
            barcode = self.spatMatDots.pop(pos)
            self.keys.pop(pos)
            self.keys.insert(pos + n - 1, new_row_ave)
            self.spatMatDots.insert(pos + n - 1, barcode)

        # if new taxicab distance is now smaller than what used to be the next smallest
        elif pos != 0 and self.keys[pos - 1] > new_row_ave:
            print('smaller')
            # find how many positions this spatial dot should move down in the list
            n = -1
            while pos + n and self.keys[pos + n] > new_row_ave:
                n -= 1
            # remove and reinsert barcode and key into their lists (inserts to left of insertion position)
            barcode = self.spatMatDots.pop(pos)
            self.keys.pop(pos)
            self.keys.insert(pos - n + 1, new_row_ave)
            self.spatMatDots.insert(pos - n + 1, barcode)


class spatMatDot:
    '''
    Stores a spatially matched Dot. What hybs each dot is from and information about each dot.
    '''
    maxError = 1  # pixels
    maxErrorSq = maxError ** 2
    nHybs = None  # set in constructor of psfDecoder

    def __init__(self, hyb, psf):
        self.hybs = [0] * spatMatDot.nHybs
        self.psfs = [None] * spatMatDot.nHybs
        self.row_ave = psf['row']
        self.col_ave = psf['row']
        self.hybs[hyb] = 1
        self.psfs[hyb] = psf
        self.nPSFs = 1

        # for sorting
        self.taxiCabDist_ave = self.row_ave + self.col_ave

    def checkMatch(self, psf):
        '''
        Checks if a psf center is within error of the average center of psfs in this spatially matched dot.
        :param psf: pandas.Series object containing psf parameters
        :return:
        '''
        return (psf['row'] - self.row_ave) ** 2 + (psf['col'] - self.col_ave) ** 2 < spatMatDot.maxErrorSq

    def addHybPSF(self, hyb, psf):
        '''
        Adds a psf in a given hyb to the spatially matched dot
        :param hyb: the hyb in which the spatially matched dot was found
        :param psf: pandas.Series object containing psf parameters
        :return:
        raises Exception if a psf has already been added from the same hyb
        '''
        if self.hybs[hyb]:
            print("Warning: A PSF has already been matched in hyb %d" % hyb)
            # raise Exception("A PSF has already been matched in hyb %d" % hyb)
        self.hybs[hyb] = 1
        self.psfs[hyb] = psf
        self.row_ave = (self.row_ave * self.nPSFs + psf['row']) / (self.nPSFs + 1)
        self.col_ave = (self.col_ave * self.nPSFs + psf['col']) / (self.nPSFs + 1)
        self.nPSFs += 1
        # self.taxiCabDist_ave = self.row_ave + self.col_ave
