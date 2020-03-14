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
from matplotlib import py
import pandas as pd
import seaborn as sns
import copy
from scipy.spatial import cKDTree as KDTree
from bisect import bisect_left

def index(a, x):
    'Locate the leftmost value exactly equal to x'
    i = bisect_left(a, x)
    #if i != len(a) and a[i] == x:
    return i
    #raise ValueError


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
        #self.hybs = []
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
            #self.hybs.append(mod)
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

        #os.chdir(self.directory)
        if filename[-3:] == 'pkl':
            with open(filename, 'rb') as f:
                self.PSFs = pickle.load(f)
        else:
            self.psfsDF = pd.read_csv(filename, sep=',', header=0, index_col=('Hyb', 'Dot'))

    '''
    def saveAlignment(self, filename):
        with open(filename, 'w') as f:
            json.dump(self.graphs_Sorted_Nodes,f)

    
    def loadOffsets(self, filename):
        #Loads saved alignment, possibly from Dapi Aligning matlab scripts.
        #:param filename:
        #string filename in which the offsets are stored.
        #:return:
        #saves pandas dataframe of offsets into self.offsets
        
        offsets = pd.read_csv(filename, header=1, index_col=0)
        self.offsets = pd.DataFrame()
        self.offsets['x'] = -(offsets['x'] - np.mean(offsets['x']))
        self.offsets['y'] = -(offsets['y'] - np.mean(offsets['y']))
    '''

    def find_matching_dots2(self, maxError = 0.5, nHybs=40):
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
                        self.spatMatchedDots.insert_new_spatDot(spatMatDot(hybNum,psf))



    def find_matching_dots(self, refHyb = 0, search_dist = 0.2):
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
            #self.sortedPSFs = self.alignedPSFs.sort_values(by='row')
            self.sortPSFCoord()
            self.psfSorted = True
        
        # Now match dots.
        
        # For each dot in the refHyb
        for i, psfi in enumerate(self.sortedPSFs[refHyb]):
        #for i, psfi in self.sortedPSFs.loc[0].iterrows():
            # Initialize sum of coordintates with which the average
            # of the coordinates will be calculated on each round
            #row_ave = psfi.row
            #col_ave = psfi.col
            
            imatches = [(i, psfi)]
            
            # for all of the other hybs:
            for j, hyb in enumerate(self.sortedPSFs):#[1:]):
            #for j, hyb in self.sortedPSFs.groupby(level=0):
                if j == refHyb:
                #if j == 0:
                    continue
                #find the highest and lowest index of dots within the search distnace of psfi
                lo = index(self.sortedPSFKeys[j], (psfi.row - search_dist, psfi.col))
                hi = index(self.sortedPSFKeys[j], (psfi.row + search_dist, psfi.col))
                #lo = index(self.sortedPSFKeys[j], psfi['row']-search_dist)
                #hi = index(self.sortedPSFKeys[j], psfi['row']+search_dist)
                #lo = index(np.array(hyb['row']), float(psfi['row']) - search_dist)
                #hi = index(np.array(hyb['row']), float(psfi['row']) + search_dist)
                
                #initialize list of dots within 1 pixel width of psfi
                candidates = []
                
                # for each dot within 1 row of psf 0
                for k, dot in enumerate(hyb[lo:hi]):
                #for ind, dot in hyb.iloc[lo:(hi+1)].iterrows():
                    #if the dot is within the search distance of psfi, add it to the candidates list
                    if (dot.row - psfi.row)**2 + (dot.col - psfi.col)**2 < search_dist:
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
                        else: ic += 1
                        
                    #if this shortens the list of candidates to 1, then add it, mark as matched
                    if len(candidates) == 1: 
                        if candidates[0][1].matched:
                            print("Warning: PSF at %d, %d already matched." % (candidates[0][1].row, candidates[0][1].col))
                        imatches.append(candidates[0])
                        candidates[0][1].matched = True
                        
                    # if their are still multiple candidates, add the list
                    elif len(candidates) > 1: imatches.append(candidates)
                    
                    # if all were removed, then just add originial candidates
                    else: 
                        imatches.append(removed)
                        print("Warning: multiple already matched PSFs matched to psf at %d, %d." % (psfi.row, psfi.col))

                # Otherwise, no matches :(
                else: imatches.append(None)
            #add matches to list of all matched dots
            self.matchedDots.append(imatches)

        # if there are multiple dots within 1 pixel of each other, try clustering by coordinates and brightness

        return
    '''
    plots matches, unmatched dots from hyb0, unmatched dots from hybn over the image of hybn
    '''
    def plotMatchedUnmatchedDots(self, directory=None, nhybs=40, bnds = 2048):
        if not directory:
            if not self.directory:
                raise ValueError('Please supply a directory of the images to plot.')
            else:
                directory = self.directory


        os.chdir(directory)
        # for each hyb > 0
        for hyb in range(nhybs):
            # plot matched dots
            os.chdir('HybCycle_%d' % hyb)
            im = skimage.io.imread('MMStack_Pos0.ome.tif')
            #im = im[0][:bnds, :bnds, 1]
            im = im[self.z][:bnds, :bnds, self.channel]
            os.chdir('..')
            plt.imshow(np.log2(im), cmap='Greys', label='Log2 Intensity')
            matchedx = []
            matchedy = []
            unmatchedx0 = []
            unmatchedy0 = []
            for matchedDots in self.matchedDots:
                if matchedDots[hyb] == None:
                    x = matchedDots[0][1].col
                    y = matchedDots[0][1].row

                    #shift from average coordinate to actual coordinate in this hyb
                    #x += self.mean_col_shifts[hyb]
                    #y += self.mean_row_shifts[hyb]
                    x += self.offsets['x'].iloc(0)[hyb]
                    y += self.offsets['y'].iloc(0)[hyb]

                    unmatchedx0.append(x)
                    unmatchedy0.append(y)
                elif type(matchedDots[hyb][1]) != tuple:
                    x = matchedDots[hyb][1].col
                    y = matchedDots[hyb][1].row

                    # shift from average coordinate to actual coordinate in this hyb
                    #x += self.mean_col_shifts[hyb]
                    #y += self.mean_row_shifts[hyb]
                    x += self.offsets['x'].iloc(0)[hyb]
                    y += self.offsets['y'].iloc(0)[hyb]
                    matchedx.append(x)
                    matchedy.append(y)

            # get dots matched by graph:
            graphAlignedx = []
            graphAlignedy = []
            if self.graphs_Sorted_Nodes:
                # for matchedDot in hybGraph:
                for matchedDot in self.graphs_Sorted_Nodes[hyb]:
                    #these do not need to be shifted, because they are being plotted over the image
                    # to which they were fitted
                    mdrow, mdcol = matchedDot[0]
                    #mdrow -= self.mean_row_shifts[hyb]
                    #mdcol -= self.mean_col_shifts[hyb]
                    graphAlignedx.append(mdcol)
                    graphAlignedy.append(mdrow)



            plt.plot(matchedx, matchedy, 'g.', label='Matched Dots')
            plt.plot(graphAlignedx, graphAlignedy, 'g*', label='Graph Aligned Dot')
            plt.plot(unmatchedx0, unmatchedy0, 'rx', label='Unmatched Hyb0 dots')

            unmatchedxhyb = []
            unmatchedyhyb = []
            #for dot in self.PSFs[hyb]:
            for dot in self.sortedPSFs[hyb]:
                if not dot.matched:
                    #x = dot.col + self.mean_col_shifts[hyb]
                    #y = dot.row + self.mean_row_shifts[hyb]
                    x = dot.col + self.offsets['x'].iloc(0)[hyb]
                    y = dot.row + self.offsets['y'].iloc(0)[hyb]

                    unmatchedxhyb.append(x)
                    unmatchedyhyb.append(y)

            plt.plot(unmatchedxhyb, unmatchedyhyb, 'm+', label='Unmatched hybN dots')
            plt.legend()
            cbar = plt.colorbar()
            cbar.set_label('Log2 Intensity')
            plt.show()

            # plot unmatched dots from hyb
            # plot unmatched dots from hyb0
        return

    def sortPSFCoord(self):
        '''
        sorts PSFs, makes shift corrections
        :return:
        self.sortedPSFs- list that contains list of PSFs in each hyb. The PSF lists for each hyb are sorted using
            CoordComp, that is by row than column coordinate. The coordinates are are aligned using the graph based
            alignnment.
        self.sortedPSFKeys - tuples of row, column for each PSF in sorted PSFs.
        '''


        # first, sort dots in each hyb by coordinate
        self.sortedPSFs = []
        self.sortedPSFKeys = []
        #self.graphAlignedDots = []
        for hybi, hybPSFs in enumerate(self.PSFs):
            # now, sort
            shybPSFs = sorted(hybPSFs, key = lambda psf: (psf.row, psf.col))

            # make shift correction
            #for i, hyb in enumerate(shybPSFs):
            for dot in shybPSFs:
                #dot.row -= self.mean_row_shifts[hybi]
                #dot.col -= self.mean_col_shifts[hybi]
                dot.row -= self.offsets['row_offset'].iloc(0)[hybi]
                dot.col -= self.offsets['col_offset'].iloc(0)[hybi]

            # build list of keys and initialize matched status of each psf
            shybPSFsKeys = [(psf.row, psf.col) for psf in shybPSFs]
            for psf in shybPSFs: psf.matched = False

            '''
            if self.graphs_Sorted_Nodes:
                # remove dots matched during graph alignment
                for matchedDot in self.graphs_Sorted_Nodes[hybi]:
                    #for matchedDot in hybGraph:
                    mdrow, mdcol = matchedDot[0]

                    # make shift correction from graph_Sorted_Nodes
                    #mdrow -= self.mean_row_shifts[hybi]
                    #mdcol -= self.mean_col_shifts[hybi]
                    mdrow -= self.offsets['row_offset'].iloc(0)[hybi]
                    mdcol -= self.offsets['col_offset'].iloc(0)[hybi]

                    ind = index(shybPSFsKeys, (mdrow, mdcol))
                    shybPSFsKeys.pop(ind)
                    shybPSFs.pop(ind)
            '''


            self.sortedPSFs.append(shybPSFs)
            self.sortedPSFKeys.append(shybPSFsKeys)


    def evaluateWellMatchedDots(self, directory=None, bnds=2048):
        '''
        Looks at subset of dots that were matched across all 40 hybs. Plots image of 0th hyb with fit center in 0th hyb
        and matched centers in all other hybs
        :return:
        '''

        dotGroups = []
        #graphAlignedDotGroups = [[] for i in range(len(self.sortedPSFs[0]))]
        '''
        hyb0GraphMatchedx = []
        hyb0GraphMatchedy = []
        graphMatchedx = []
        graphMatchedy = []
        if self.graphs_Sorted_Nodes:
            # get dots matched by graph:
            for hyb in range(len(self.sortedPSFs)):
                for dot, matchedDot in enumerate(self.graphs_Sorted_Nodes[hyb]):
                    # for matchedDot in hybGraph:
                    graphAlignedDotGroups[dot].append(matchedDot[0])
                    mdrow, mdcol = matchedDot[0]
    
                    #shift to average alignment
                    #mdrow -= self.mean_row_shifts[hyb] - self.mean_row_shifts[0]
                    #mdcol -= self.mean_col_shifts[hyb] - self.mean_col_shifts[0]
                    mdrow -= self.offsets['y'].iloc(0)[hyb] - self.offsets['y'].iloc(0)[0]
                    mdcol -= self.offsets['x'].iloc(0)[hyb] - self.offsets['x'].iloc(0)[0]
    
                    if hyb == 0:
                        hyb0GraphMatchedx.append(mdcol)
                        hyb0GraphMatchedy.append(mdrow)
                    else:
                        graphMatchedx.append(mdcol)
                        graphMatchedy.append(mdrow)
        '''


        # get matched by alignment in all hybs
        hyb0AlignedDotsx = []
        hyb0AlignedDotsy = []
        AlignMatchedDotsx = []
        AlignMatchedDotsy = []
        AlignMatchedDotsxGrouped = []
        AlignMatchedDotsyGrouped = []
        for dot in range(len(self.matchedDots)):
            AllMatched = True
            xdots = []
            ydots = []
            dots = []
            for hyb in range(len(self.matchedDots[0])):
                # check if for match failures
                if self.matchedDots[dot][hyb] == None or type(self.matchedDots[dot][hyb]) == list:
                    AllMatched = False
                    xdots = []
                    ydots = []
                    dots = []
                    x = None
                    y = None
                    break
                else:
                    #x = self.matchedDots[dot][hyb][1].col + self.mean_row_shifts[0]
                    #y = self.matchedDots[dot][hyb][1].row + self.mean_col_shifts[0]
                    #x = self.matchedDots[dot][hyb][1].col + self.offsets['y'].iloc(0)[0]
                    #y = self.matchedDots[dot][hyb][1].row + self.offsets['x'].iloc(0)[0]
                    if hyb == 0:
                        hyb0dotx = x
                        hyb0doty = y
                        dots.append(self.matchedDots[dot][hyb][1])
                    #else:
                    xdots.append(x)
                    ydots.append(y)
                    dots.append(self.matchedDots[dot][hyb][1])

            if AllMatched:
                hyb0AlignedDotsx.append(hyb0dotx)
                hyb0AlignedDotsy.append(hyb0doty)
                AlignMatchedDotsx += xdots
                AlignMatchedDotsy += ydots
                AlignMatchedDotsxGrouped.append(xdots)
                AlignMatchedDotsyGrouped.append(ydots)
                dotGroups.append(dots)



        
        '''
        hyb0_alignMatched = []
        for matches in AllMatchedDots:
            hyb0_alignMatched.append(matches.pop())
            
        hyb0_graphMatched = []
        for matches in AllMatchedDots:
            hyb0_alignMatched.append(matches.pop())
        '''
        #print('Number of graph matched:', len(hyb0GraphMatchedx))
        print("Number of matched:", len(hyb0AlignedDotsx)) #+ len(hyb0GraphMatchedx))

        #find statistics on each dot group
        xs = []
        ys = []
        amps = []
        xstdv = []
        ystdv = []
        xMean = []
        yMean = []
        ampstdv = []
        ampmean = []
        RMSESTDV = []
        RMSEMean = []
        minAmp = np.inf
        for dotGroup in dotGroups:
            x = []
            y = []
            amp = []
            RMSE = []
            for dot in dotGroup:
                x.append(dot.col)
                y.append(dot.row)
                amp.append(dot.amp)
                RMSE.append(dot.RMSE)
                if dot.amp < minAmp: minAmp = dot.amp
            ampmean.append(np.mean(amp))
            ampstdv.append(np.std(amp))
            xs.append(x)
            ys.append(y)
            amps.append(amp)
            xstdv.append(np.std(x))
            ystdv.append(np.std(y))
            xMean.append(np.mean(x))
            yMean.append(np.mean(y))
            RMSESTDV.append(np.std(RMSE))
            RMSEMean.append(np.mean(RMSE))

        ampmean_div_RMSE = np.array(ampmean)/np.array(RMSEMean)

        summaryStats = pd.DataFrame({
            'XSTDV': xstdv,
            'YSTDV': ystdv,
            'XMean': xMean,
            'YMean': yMean,
            'ampSTDV': ampstdv,
            'Mean_Brightness': ampmean,
            'RMSE_Mean': RMSEMean,
            'RMSE_STDV': RMSESTDV,
            'Mean_Brightness_div_Mean_RMSE': ampmean_div_RMSE
        })

        print(minAmp)
        sns.scatterplot(x='XSTDV', y='YSTDV', hue='Mean_Brightness', size='Mean_Brightness_div_Mean_RMSE',
                        sizes=(10,100), data=summaryStats)
        plt.xlabel('Standard Deviation of X Center Fits')
        plt.ylabel('Standard Deviation of Y Center Fits')
        plt.title('Summary Statistics of dots matched in all 40 Hybridizations')
        plt.show()

        amp_zscores = []
        for i, dotgroup in enumerate(dotGroups):
            for dot in dotgroup:
                amp_zscores.append((dot.amp-ampmean[i])/ampstdv[i])


        if not directory:
            if not self.directory:
                raise ValueError('Please supply a directory of the images to plot.')
            else:
                directory = self.directory

        os.chdir(directory)
        os.chdir('HybCycle_0')

        #im = skimage.io.imread('MMStack_Pos0.ome.tif')
        im = skimage.io.imread(self.imfname)
        im = im[0][:bnds, :bnds, 1]
        plt.imshow(np.log2(im), cmap='Greys')


        for i in range(len(AlignMatchedDotsyGrouped)):
            plt.plot(AlignMatchedDotsxGrouped[i], AlignMatchedDotsyGrouped[i], '-k')
        plt.plot(AlignMatchedDotsx, AlignMatchedDotsy, '.g', label='align Matched peaks')
        '''
        if self.graphs_Sorted_Nodes:
            plt.plot(graphMatchedx, graphMatchedy, 'g*', label='graph Matched peaks')
        plt.plot(hyb0AlignedDotsx, hyb0AlignedDotsy, '.b', label='Hyb0 align Matched peaks')
        if self.graphs_Sorted_Nodes:
            plt.plot(hyb0GraphMatchedx, hyb0GraphMatchedy, 'r*', label='Hyb0 graph Matched peaks')
        '''

        cbar = plt.colorbar()
        cbar.set_label('Log2 Intensity')
        plt.legend()
        plt.show()
        return amp_zscores, summaryStats, np.array(xs), np.array(ys), np.array(amps)

    def plotHybPSFs(self, hyb, log = True):
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
        ypsfs= []
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
        os.chdir(self.directory)
        if df:
            self.refPSFs = pd.read_csv(refPSFs)
        else:
            with open(refPSFs, 'rb') as f:
                refPSFs = pickle.load(f)
                self.refPSFs = self.convtPSFList2DF([refPSFs])


    def find_offsets_from_ref(self, maxError=0.1, refFirst = False, trav=False, brightThresh= None):
        '''
        Align images using bead
        :param maxError: maximum allowed difference between matched displacement vectors between beads in reference
        image and in hyb image
        :return:
        data frame of offsets and offset standard errors
        '''
        #self.convtPSFList2DF()
        print('Initiaizing Bead Graph')
        if refFirst:
            self.refGraph = refDotGraph(self.psfsDF.loc[0], brightThresh=brightThresh, maxError=maxError)
            self.matchedRefDots = self.psfsDF.loc[0].copy
        else:
            self.refGraph = refDotGraph(self.refPSFs, brightThresh=brightThresh, maxError=maxError)
            self.matchedRefDots = self.refPSFs.copy()
        #self.refGraph = refDotGraph(self.refPSFs)
        if not trav:
            self.refGraph.sortEdgesRows()


        row_offset = []
        col_offset = []
        row_offset_SE = []
        col_offset_SE = []
        #for i, hybPSFs in enumerate(self.PSFs):
        for hyb, hybPSFs in self.psfsDF.groupby(level=0):
            #print('Initializing Hyb %d Graph' % i)
            #hybGraph = dotGraph(self.PSFs[i])
            print('Aligning hyb %d' % hyb)

            if trav:
                r_mean, c_mean, r_mean_se, c_mean_se, r_weight_mean, c_weight_mean = self.refGraph.alignWith2(hybPSFs,
                                                                                                         maxError)
            else:
                hybGraph = dotGraph(hybPSFs, brightThresh) #self.PSFs[i])
                hybGraph.sortEdgesRows()
                r_mean, c_mean, r_mean_se, c_mean_se, r_weight_mean, c_weight_mean = self.refGraph.alignWith(hybGraph, maxError)
            print('Offsets and SEs:', r_mean, c_mean, r_mean_se, c_mean_se)
            row_offset.append(r_mean)
            col_offset.append(c_mean)
            row_offset_SE.append(r_mean_se)
            col_offset_SE.append(c_mean_se)
        self.offsets = pd.DataFrame()
        self.offsets['row_offset'] = row_offset
        self.offsets['col_offset'] = col_offset
        self.offsets['row_offset_SE'] = row_offset_SE
        self.offsets['col_offset_SE'] = col_offset_SE
        return self.offsets

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
        self.offsets = pd.read_csv(filename, header = 0)

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
                #row.append(psf.row + self.offsets['row_offset'][hyb])
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
            pd.concat((self.alignedPSFs, hybPSFs), axis=0)


    def convtPSFList2DF(self, PSFs = None):
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
        self.alignedPSFs = pd.read_csv(filename, sep=',', header=0, index_col=('Hyb', 'Dot'))

    def compareHybAlignmentToRef(self, hyb, ref0=False):
        '''
        Plots the bead image with the bead peaks with aligned peaks from the given hyb number
        :param hyb: int number of the hyb to compare alignment to beads
        :return:
        '''
        os.chdir('beads')
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
        #ref_rows = 2048 - np.array(ref_rows)
        plt.plot(ref_cols, ref_rows, 'xb', label='Reference Image Fit Centers')


        #self.alignedPSFs.loc[hyb].plot(x='row', y='col', )
        plt.legend()
        cbar = plt.colorbar()
        cbar.set_label('Log2 Intensity')
        #plt.plot()
        plt.show()








'''
MatchedDots holds stores the dots from each hybs groups according to the ones that correspond
         each other

class MatchedDots:
    
    def __init__(self):
        self.dots = []
        self.kdots = {}
        
    def match(self, psf, hyb):
        
    def addKnown(self, key, psf):
        if key in self.kdots:
            self.kdots[key].append(psf)
        else:
            self.kdots[key] = [psf]
'''        
                       

def edgesEquivalent(e1, e2, maxError):
    '''
    Tests whether the edge displacement vectors are equal to an error of maxError pixels.
    :param e1: dictionary of length attributes of first edge {l: a, rd: b, cd, c}
    :param e2: dictionary of length attributes of second edge {l: a, rd: b, cd, c}
    :param maxError: maximum allowed difference between edge displacement vectors
    :return:
    boolean; true if equal to within error, false otherwise
    '''
    #if xdist, ydist, and absolute distances are all within 1 pixel, return True
        
    #if np.abs(e1['l'] - e2['l']) < 1.4 and np.abs(e1['rd'] - e2['rd']) < 1 and np.abs(e1['cd'] - e2['cd']) < 1:
    if  np.abs(e1['rd'] - e2['rd'])**2 + (e1['cd'] - e2['cd'])**2 < maxError:
        return True
    # otherwise, return false
    else:
        return False

# assigns a sortable value to an edge
'''
Inputs:
    dists - tuple (row_dist, col_dist)
'''
def edgeValue(dists):
    # row distance is always positive, so make value negative if column distance is negative
    if dists[1] >= 0:
        return dists[1] + dists[0]*2048
    else:
        return dists[1] - dists[0]*2048
    

'''
#def index(a, x, lo=(0, 0), hi=None):
def index(a, x, lo=0, hi=None):
    '''
#find PSF in sorted list based on example code from python docs:
#https://docs.python.org/2/library/bisect.html

'''
    #if coordComp(lo, (0,0)):#lo[0] < 0 and ; lo[1] < 0:
    if lo < 0:
        raise ValueError('components of lo must be non-negative')
    if hi is None:
        hi = len(a)#(2048, 2048)
    while lo < hi:
        mid = (lo+hi)//2
        #if a[mid] < x: lo = mid+1
        if coordComp(a[mid], x): lo = mid+1
        else: hi = mid
    #a.insert(lo, x)
    return lo
'''
    
'''
If c1 is less (row than column) than c2, returns true, otherwise, returns false.
'''
def coordComp(c1, c2):
    if c1[0] < c2[0]:
        return True
    elif c1[0] == c2[0] and c1[1] < c2[1]:
        return True
    else:
        return False


class dotGraph(nx.Graph):
    '''
    Inherits from nx.Graph. ALso contains data and methods for sorting and matching.
    '''
    def __init__(self, PSFs=None, brightThresh=None, nBrightest=None):
        super(dotGraph,self).__init__()
        self.edgesXSort = []
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
        sPSFs = PSFs.sort_values(by='amp', ascending=False)
        if brightThresh:
            sPSFs = sPSFs[sPSFs['amp'] > brightThresh]

        self.graph = nx.Graph()

        # print('Sorted PSFs, going on to iterative adding...')

        # Add add the dots to graph as a node
        # including its brightness and coordinates
        # all nodes are connected by edges holding the displacement vector between them.
        self.disp_vecs = []
        self.dispVecEdges = []
        # for psf in sPSFs:

        # for i, psf in self.sPSFs.iterrows():
        for i, psf in sPSFs.iterrows():
            # self.add_node((psf.row, psf.col), amp=psf.amp)
            self.add_node((psf['row'], psf['col']), amp=psf['amp'])
            # self.graph.add_node((psf['row'], psf['col']), amp=psf['amp'])
            # add edges between the new node and the already added nodes
            # for each node in the graph
            for node in self.__iter__():
                # if it is not the node that was just added
                # if node != (psf.row, psf.col):
                if node != (psf['row'], psf['col']):
                    # calculate the length of the node. First the euclidean length
                    # length = np.sqrt((node[0] - psf.row) ** 2 + (node[1] - psf.col) ** 2)
                    # then the row and column distances between the node in a smaller row and the node in a larger row
                    # to ensure the row component of the displacement vector is always positive
                    # row_dist, col_dist = self.findDispVec(node[0], node[1], psf.row, psf.col)
                    row_dist, col_dist = self.findDispVec(node[0], node[1], psf['row'], psf['col'])
                    '''
                    if node[0] < psf.row:
                        row_dist = psf.row - node[0]
                        col_dist = psf.col - node[1]
                    else:
                        row_dist = node[0] - psf.row
                        col_dist = node[1] - psf.col
                    '''
                    # add an edge between the PSF added and the the node that was already in the graph
                    # self.add_edge(node, (psf.row, psf.col), l=length, rd=row_dist, cd=col_dist)
                    self.add_edge(node, (psf['row'], psf['col']), rd=row_dist, cd=col_dist)
                    # self.graph.add_edge(node, (psf['row'], psf['col']), rd=row_dist, cd=col_dist)
                    self.disp_vecs.append((row_dist, col_dist))
                    # self.dispVecEdges.append([(node[0], node[1]), (psf.row, psf.col)])
                    # add
                    if node[0] < psf['row']:
                        self.dispVecEdges.append(((node[0], node[1]), (psf['row'], psf['col'])))
                    else:
                        self.dispVecEdges.append(((psf['row'], psf['col']), (node[0], node[1])))

        self.disp_vecs = np.array(self.disp_vecs)

    def buildKDTree(self):
        print('Building KDTree')
        self.edgeKDTree = KDTree(self.disp_vecs)


    def sortEdgesRows(self):
        # make list of edges sorted by row distance
        edgeList = []
        for (u, v, eattrs) in self.edges.data():
            edgeList.append((u, v, eattrs))

        self.rsortedEdges = sorted(edgeList, key=lambda edge: edge[2]['rd'])
        
    def findDispVec(self, r1, c1, r2, c2):
        '''
        Find dispacement vector such that the row component is always positive.
        :param r1: row coordinate of first psf
        :param c1: column coordinate of first psf
        :param r2: row coordinate of second psf
        :param c2: column coordinate of second psf
        :return:
        row - the row component of the displacement vector between two psfs (always positive)
        col - the column component of the displacement vector between two psfs
        '''
        #if  < psf.row:
        if r1 < r2:
            row_dist = r2 - r1 #psf.row - node[0]
            col_dist = c2 -c1 #psf.col - node[1]
        else:
            row_dist = r1 - r2 #node[0] - psf.row
            col_dist = c1 - c2 #node[1] - psf.col
        return row_dist, col_dist

    def sortDotsBrightness(self, PSFs, brightThresh = None):
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
                edgeList = tuple(sorted(edgeList))#, key=edgeValue))
                self.sortedNodes.append((node, edgeList))
        self.sortedNodes = sorted(self.sortedNodes, key=itemgetter(1))
        #self.graphs_Sorted_Nodes.append(nodeList)
        for node in toRemove:
            self.remove_node(node)
        return self.sortedNodes

class refDotGraph(dotGraph):
    '''
    Inherits from dotGraph. Holds the dot graph for the reference image for alignment.
    Contains methods for aligning with other images.
    '''
    def __init__(self, PSFs, brightThresh=3000, nBrightest=None, maxError = 0.1):
        print('Building Reference Graph')
        super(refDotGraph, self).__init__(PSFs, brightThresh, nBrightest)
        self.maxError = maxError



    def lookUpDispVec(self, rdist, cdist):
        matches = self.edgeKDTree.query_ball_point((rdist, cdist), self.maxError)
        disp_vec_matches = self.disp_vecs[matches]
        edges = []
        for match in matches:
            edges.append(self.dispVecEdges[match])
        return edges, disp_vec_matches

    

    def alignWith(self, dGraph, maxError = 0.1):
        '''
        Finds alignment with another dot Graph
        :param dGraph:
        :return:
        '''
        print('Finding intersection')
        matchDict = self.findEdgeIntersection(dGraph, maxError=maxError)
        offset = self.estOffset(matchDict)
        #refMatchGraph, compMatchGraph = self.buildMatchedGraphs(matches)
        #offset = self.estimateOffset(refMatchGraph, compMatchGraph, maxError)
        #return offset, matches, refMatchGraph, compMatchGraph, matchDict
        return offset

    def alignWith2(self, compPSFs, brightThresh=None):
        '''
        Reference dots with dots from a comparison image
        :param compPSFs: DataFrame of dots for a hybridization formatatted like output of convtPSFList2DF
        :param brightThresh: Ignore dots dimmer than this
        :return:
        Offset of hyb from reference
        '''
        self.buildKDTree()
        print('Starting Alignment')
        thrCompPSFs = compPSFs[compPSFs['amp'] > brightThresh]
        matchDict = self.findMatches(thrCompPSFs, brightThresh)
        offset = self.estOffset(matchDict)
        return offset

    def findEdgeIntersection(self, dGraph, maxError = 0.2, n=None):
        '''
        Looks for matching edges in reference and comparison graphs. Keeps track of nodes that are matched in matching
        edges in the matchDict dictionary. See return for description of matchDict.
        :param dGraph: another dot graph with which to find the intersection
        :param n: the number of the dot graph comparing to
        :return:
        dictionary of matched nodes. Keys are tuples containing first a tuple of the coordinates of the dot in the
        reference graph, then a tuple of the coordinents in the comparison dot. The value is the number of edges
        each pair was matched in.
        '''

        ref_rsortedEdges = copy.copy(self.rsortedEdges)
        comp_rsortedEdges = copy.copy(dGraph.rsortedEdges)

        maxErrorSq = maxError**2

        matches = []
        multimatches = []
        matchDict = {}
        print('n_refEdges:', len(ref_rsortedEdges), 'n_compEdges:', len(comp_rsortedEdges))

        # for each edge in dGraph
        for edge in comp_rsortedEdges:
            # look for matches within error
            redge = 0
            edge_ref_matches = []
            #print('Matching Edge:', edge)
            while redge < len(ref_rsortedEdges):
                # if the reference edge is mre than the maxError rows lower than the edge under consideration
                if ref_rsortedEdges[redge][2]['rd'] < edge[2]['rd'] - maxError:
                    # remove the reference edge from consideration
                    ref_rsortedEdges.pop(0)

                # if the reference edge is more than the maxError rows higher than the edge under consideration
                elif ref_rsortedEdges[redge][2]['rd'] > edge[2]['rd'] + maxError:
                    # no more matches will be found, move onto the next edge
                    break

                # if the reference edge is within error of the edge under consideration
                elif (ref_rsortedEdges[redge][2]['rd'] - edge[2]['rd'])**2 + (ref_rsortedEdges[redge][2]['cd'] - edge[2]['cd'])**2 < maxErrorSq:
                    # add to matches
                    edge_ref_matches.append(ref_rsortedEdges[redge])
                redge += 1
                
            if len(edge_ref_matches) > 1:
                # There is more than one match! consider separately
                multimatches.append([edge_ref_matches, edge])

            elif len(edge_ref_matches) == 1:
                matches.append([edge_ref_matches[0], edge])
                if (edge_ref_matches[0][0], edge[0]) in matchDict:
                    matchDict[(edge_ref_matches[0][0], edge[0])] += 1
                else:
                    matchDict[(edge_ref_matches[0][0], edge[0])] = 1

            # otherwise, no matches

        return matchDict

    def findMatches(self, compPSFs, brightThresh):
        '''
        Find dots from an image that match those in the reference image. To do this, first find a position
        vector between dots in the image that matches a position vector in the reference image to match the dots
        on either end of the the position vector. Then, look for matches to the known neighbors of the dots in the reference
        image by following the displacement vectors stored in its edges.
        :param PSFs:
        :param brightThresh:
        :return:
        '''
        matchDict = {}
        edgeMatchDict = {}
        #sPSFs = self.sortDotsBrightness(PSFs, brightThresh)
        sCompPSFs = compPSFs.sort_values(by='amp', ascending=False)
        compPSFsTree = KDTree(sCompPSFs.loc[:,'row':'col'])
        nPSFs = len(sCompPSFs)
        found_enough = False
        # First, look for a position vector that matches a position vector in the reference image
        for i in range(nPSFs):
            for j in range(i):
                iPSF = sCompPSFs.iloc[i]
                jPSF = sCompPSFs.iloc[j]
                rdist, cdist = self.findDispVec(iPSF['row'], iPSF['col'], jPSF['row'], jPSF['col'])
                #check whether there is a matching position in the reference image
                #match = self.edgeKDTree.query_ball_point((rdist, cdist), self.maxError)
                edgeMatch, dispVecMatches = self.lookUpDispVec(rdist, cdist)
                if len(edgeMatch) == 1 and edgeMatch[0] not in edgeMatchDict:
                    # Start graph traversal which adds matches to matchDict
                    edgeMatchDict[edgeMatch[0]] = 1
                    self.travRefFindMatches(sCompPSFs, compPSFsTree, edgeMatch, dispVecMatches, iPSF, jPSF, matchDict, edgeMatchDict)
                    if len(matchDict) > self.number_of_nodes()/2:
                        found_enough = True
                        break
            if found_enough:
                break
        return matchDict


    def travRefFindMatches(self, compPSFs, compPSFsTree, match, dispVecMatches, cPSF1, cPSF2, matchDict, edgeMatchDict):
        '''

        :param compPSFs:
        :param compPSFsTree:
        :param match:
        :return:
        '''
        print('Starting Ref Graph Traversal')

        #add 2 tuples where the first element is a dot in the reference image, and the second is a dot matched to it
        # in the comparison image
        dotQueue = []

        #the psf with smaller row value is matched to the first node in the edge
        if cPSF1['row'] < cPSF2['row']:
            #dotQueue.append((match[0][0], cPSF1))
            #dotQueue.append((match[0][1], cPSF2))
            match1Key = (match[0][0], (cPSF1['row'], cPSF1['col']))
            match2Key = (match[0][1], (cPSF2['row'], cPSF2['col']))
            self._placeMatchInDict(match1Key, matchDict, dotQueue)
            self._placeMatchInDict(match2Key, matchDict, dotQueue)
        else:
            #dotQueue.append((match[0][0], cPSF2))
            #dotQueue.append((match[0][1], cPSF1))
            match1Key = (match[0][0], (cPSF2['row'], cPSF2['col']))
            match2Key = (match[0][1], (cPSF1['row'], cPSF1['col']))
            self._placeMatchInDict(match1Key, matchDict, dotQueue)
            self._placeMatchInDict(match2Key, matchDict, dotQueue)

        # Now traverse the graph
        while len(dotQueue) > 0:
            matchedDot = dotQueue.pop(0)

            nUnmatched = 0
            nMatched = 0
            # look for each of the dot's known neighbors in the reference graph in the comparison image
            for refNeighbor in self[matchedDot[0]]:
                # if edge not already matched:
                if (matchedDot[0], refNeighbor) not in edgeMatchDict and (refNeighbor, matchedDot[0]) not in edgeMatchDict:
                    r_to_neighbor = refNeighbor[0] - matchedDot[0][0]
                    c_to_neighbor = refNeighbor[1] - matchedDot[0][1]
                    matches = compPSFsTree.query_ball_point((matchedDot[1][0] + r_to_neighbor, matchedDot[1][1] + c_to_neighbor), self.maxError)
                    # check if only one match, ignore if otherwise (ambiguous)
                    
                    #print('nMatches:', len(matches))
                    if len(matches) == 1:
                        newMatchedPSF = compPSFs.iloc[matches]
                        #print('matchedDot:', matchedDot, ', refNeighbor:', refNeighbor)
                        #print('newMatchedPSF:', newMatchedPSF)
                        newMatchedKey = (refNeighbor, (float(newMatchedPSF['row']), float(newMatchedPSF['col'])))
                        self._placeMatchInDict(newMatchedKey, matchDict, dotQueue)
                        nMatched += 1
                    elif len(matches) == 0:
                        nUnmatched += 1
                if nUnmatched > 10 and nMatched == 0:
                    #give up, this isn't a real match
                    break
                if matchDict[matchedDot] > minMa
                    # found enough matches
                    break

                        

        print('Number of matched Dots:', len(matchDict), 'of %d in reference.' % self.number_of_nodes())


    def _placeMatchInDict(self, match, matchDict, dotQueue):
        if match in matchDict:
            matchDict[match] += 1
        else:
            #print("Adding to dotQueue:", match)
            matchDict[match] = 1
            dotQueue.append(match)


    def estOffset(self, matchDict):
        '''
        Given a dictonary of matching nodes between the reference and comparison graph, estimates the offset
        between the two images.
        :param matchDict:
        :return:
        r_mean - mean row offset between matches
        c_mean- mean column offset between matches
        r_mean_se - standard error of row offset between matches
        c_mean_se - standard error of column offset between matches
        r_weight_mean - weighted mean (weighted by number of matches) of row offsets
        c_weight_mean - weighted mean (weighted by number of matches) of column offsets
        '''
        rdisp = []
        cdisp = []
        nmatches = []
        for match in matchDict:
            if matchDict[match] > 3:
                rdisp.append(match[1][0] - match[0][0])
                cdisp.append(match[1][1] - match[0][1])
                nmatches.append(matchDict[match])

        dot_r_offsets = np.array(rdisp)
        dot_c_offsets = np.array(cdisp)
        n_dot_edges = np.array(nmatches)

        # remove any outliers that my have been introduced by mismatches
        outlier = True
        while outlier:
            # compute summary stats
            r_mean = np.mean(dot_r_offsets)
            c_mean = np.mean(dot_c_offsets)
            rstdv = np.std(dot_r_offsets)
            cstdv = np.std(dot_c_offsets)

            # check if there are outliers
            outliers = np.logical_or(abs(dot_r_offsets - r_mean) > 3 * rstdv, abs(dot_c_offsets - c_mean) > 3 * cstdv)
            outlier = np.any(outliers)

            # remove outliers
            dot_r_offsets = dot_r_offsets[np.logical_not(outliers)]
            dot_c_offsets = dot_c_offsets[np.logical_not(outliers)]
            n_dot_edges = n_dot_edges[np.logical_not(outliers)]

        r_mean_se = rstdv / np.sqrt(len(dot_r_offsets))
        c_mean_se = cstdv / np.sqrt(len(dot_c_offsets))
        r_weight_mean = np.sum(dot_r_offsets * n_dot_edges) / np.sum(n_dot_edges)
        c_weight_mean = np.sum(dot_c_offsets * n_dot_edges)/ np.sum(n_dot_edges)

        return r_mean, c_mean, r_mean_se, c_mean_se, r_weight_mean, c_weight_mean


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
        while pos < len(self.spatMatDots) and np.abs(self.spatMatDots[pos].taxiCabDist_ave - psf['row']) < spatMatDot.maxError*np.sqrt(2):
            if self.spatMatDots[pos].checkMatch(psf):
                self._addPSFToSpatDot(pos, psf)
                found = True
                break
            pos += 1

        #now check for barcodes with centers with shorter taxicab distances than the psf's
        if not found:
            pos = closest - 1
            while pos >= 0 and np.abs(self.spatMatDots[pos].taxiCabDist_ave - psf['row']) < spatMatDot.maxError*np.sqrt(2):
                if self.spatMatDots[pos].checkMatch(psf):
                    self._addPSFToSpatDot(pos, psf)
                    found = True
                    break
                pos -= 1
        return found

    def _addPSFToSpatDot(self, pos, psf):
        # first add the psf to the barcode
        self.spatMatDots[pos].addHybPSF(self.curr_hyb, psf)

        #now check that the barcode is still sorted properly
        new_row_ave = self.spatMatDots[pos].row_ave

        # if new taxicab distance is now larger than what used to be the next largest
        if pos != len(self.keys)-1 and self.keys[pos+1] < new_row_ave:
            print('larger')
            # find how many positions this barcode should move up in the list
            n = 2
            while self.keys[pos+n] < new_row_ave and pos+n < len(self.spatMatDots):
                n += 1
            # remove and reinsert spatial dot and key into their lists (inserts to left of insertion position)
            barcode = self.spatMatDots.pop(pos)
            self.keys.pop(pos)
            self.keys.insert(pos + n - 1, new_row_ave)
            self.spatMatDots.insert(pos + n - 1, barcode)

        # if new taxicab distance is now smaller than what used to be the next smallest
        elif pos != 0 and self.keys[pos-1] > new_row_ave:
            print('smaller')
            # find how many positions this spatial dot should move down in the list
            n = -1
            while pos + n and self.keys[pos+n] > new_row_ave:
                n -= 1
            # remove and reinsert barcode and key into their lists (inserts to left of insertion position)
            barcode = self.spatMatDots.pop(pos)
            self.keys.pop(pos)
            self.keys.insert(pos-n+1, new_row_ave)
            self.spatMatDots.insert(pos - n + 1, barcode)

class spatMatDot:
    '''
    Stores a spatially matched Dot. What hybs each dot is from and information about each dot.
    '''
    maxError = 1 #pixels
    maxErrorSq = maxError**2
    nHybs = None # set in constructor of psfDecoder
    def __init__(self, hyb, psf):
        self.hybs = [0] * spatMatDot.nHybs
        self.psfs = [None] * spatMatDot.nHybs
        self.row_ave = psf['row']
        self.col_ave = psf['row']
        self.hybs[hyb] = 1
        self.psfs[hyb] = psf
        self.nPSFs = 1

        #for sorting
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
            #raise Exception("A PSF has already been matched in hyb %d" % hyb)
        self.hybs[hyb] = 1
        self.psfs[hyb] = psf
        self.row_ave = (self.row_ave*self.nPSFs + psf['row'])/(self.nPSFs+1)
        self.col_ave = (self.col_ave*self.nPSFs + psf['col'])/(self.nPSFs+1)
        self.nPSFs += 1
        #self.taxiCabDist_ave = self.row_ave + self.col_ave
