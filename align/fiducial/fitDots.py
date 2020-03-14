# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 12:10:16 2019

@author: jonat
"""

import numpy as np
from scipy.optimize import least_squares
from fit_gaussian import gaussian2D
from skimage.feature import peak_local_max
from matplotlib import pyplot as plt
import skimage
import copy
import pickle
#import cv2 as cv


'''
The MixMod class fits a guassian mix model to an image and acts as in interface to the model
results and its parameters
'''
class dotFinder:
    
    '''
    The MixMod constructor takes and save the image to fit, the width parameter
    of the imaging system's Point Spread function, and the type of PSF used.
    Inputs:
        im - the image to fit
        PSFwidth - the width of the PSF of the system
        PSF_type - the type of PSF used in the model. Gaussian by default.
    '''
    def __init__(self, im, bg, PSFwidth, PSF_type = 'Gaussian', meanBgPropThresh = 3/4):
        self.im = im
        self.bg = bg
        self.PSFwidth_guess = PSFwidth
        self.meanBgPropThresh = meanBgPropThresh
        
        fitPSF.width_guess = PSFwidth
        
        if PSF_type == 'Gaussian':
            self.PSF = gaussian2D
        else:
            raise ValueError(PSF_type + ' is not a valid PSF.')
        
        # initialize list of single PSFs that fit well on their own without 
        # inclusion in a mixture model
        self.allAcceptedPSFs = []
        self.singlePSFs = []
        self.singlePSFsMod = np.zeros(np.shape(im))
        self.singlePSFAmps = []
        
        # initialize list of peaks that do not fit well on their own to be
        # to be included in mixture model
        self.candidate_overlapping_peaks = []
        
        # initialize list of PSFs added to fill hole in model where nearby
        # peak marked as overlapping does not fit well in initial mixture model
        self.nonPeakedPSFs = []
        
        
        self.wfholes = []
        self.olholes = []
        self.mixModsWHoles = []
    '''
    fitSinglePSFs fits a single PSF around each maximum found in the image poorly
        fitting PSFs are rejected. PSFs that fit well are kept, and those that 
        would wander more than 1 pixel from their maximum a flagged for the mix model.
        
    Inputs:
        im - if None, performs algorithm on self.im. Otherwise, on the passed image
        threshold - rejects local maxima less than this
    
    Returns:
        model_im - a 2D numpy array of the same dimensions as the orignial image
                with the fitted Gaussians plotted on it.
        singlePSFs - a list of well fitting fitPSF objects containing information
                    about gaussians accepted to the model
        overlapping_peaks - a list of PSF objects that are bright and wandered
                        into their bounds during fitting. Will be given a longer
                        leash in the Mixture model.
        
    '''
    def fitSinglePSFs(self, im=None, threshold=None):
        if type(im) != np.ndarray:# and im == None:
            im = self.im
            model_im = self.singlePSFsMod
            singlePSFs = self.singlePSFs
            overlapping_peaks = self.candidate_overlapping_peaks
            significance_factor = 0.25 # 0.25 seemed decent
            candidate_amps = self.singlePSFAmps
            first = True
            
        else:
            model_im = np.zeros(np.shape(im))
            singlePSFs = []
            overlapping_peaks = []
            significance_factor = 1 # higher standard for accepting holes
            candidate_amps = []
            first = False
            print('not first')

        self.subtractBackground()

        im_filt = self.im_filt
        
        
        # find peaks
        loc_maxes = peak_local_max(im_filt, indices = True)
        if threshold:
            print('filtering')
            loc_maxes = [peak for peak in loc_maxes if im[peak[0], peak[1]] > threshold]
        
        if first:
            self.loc_maxes = loc_maxes
        
        candidate_peaks = []
        
        
        
        #for each peak
        for loc_max in loc_maxes:
            #print(loc_max)
            

            # fit psf around peak. Only pixels within 1.3 PSF width (rounded up) 
            # of center should be considered when evaluating the fit
            peak_fit = fitPSF(self.imbgsub, *loc_max)
            #peak_fit = fitPSF(im_filt, *loc_max)
            
            
            # if amplitude insignificant, do not include in model
            #if not peak_fit.amplitude_significant:
            
            # be pretty lenient with the fit since we do not expect
            if peak_fit.RMSE < 0.3*peak_fit.amp:
                # then peak is a candidate for including in the model
                #self.thrown_out.append(peak_fit)
                candidate_peaks.append(peak_fit)
                candidate_amps.append(peak_fit.amp)
        
        print('Number of Candidate Peaks:', len(candidate_peaks))
        candidate_amps = np.array(candidate_amps)
        mean_amp = np.mean(candidate_amps)
        for candidate in candidate_peaks:
            # if any bounds are active, flag for mix model
            if np.any(candidate.fit.active_mask[2:4] != 0) and candidate.RMSE < candidate.amp:
                #print('active bound')
                if not first:
                    print('hole active bound')
                    # don't accept dim holes
                    if (candidate.amp > np.percentile(self.singlePSFAmps, 25) and candidate.amp > 500):
                        # don't accept holes that are within a PSF width of a previously fitted peak
                        cleared = True
                        for psf in self.allAcceptedPSFs:
                            dist = np.sqrt((candidate.row - psf.row)**2 + (candidate.col - psf.col)**2)
                            print(dist)
                            if dist < 1.5:
                                print('too close')
                                cleared = False
                                break
                        if cleared:
                            singlePSFs.append(candidate)
                            candidate.paint(model_im)
                            self.allAcceptedPSFs.append(candidate)
                # if the boundary is active, only keep if fairly bright.
                elif mean_amp/4 < candidate.amp:
                    overlapping_peaks.append(candidate)
                    self.allAcceptedPSFs.append(candidate)
                #else:
                #    continue
            # otherwise no bounds are active and the fit is at least decent: 
            # accept to single PSFs list
            elif candidate.RMSE < significance_factor*candidate.amp:
                if first:
                    singlePSFs.append(candidate)
                    candidate.paint(model_im)
                    self.allAcceptedPSFs.append(candidate)
                else:
                    print('hole')
                    # don't accept dim holes
                    if (candidate.amp > np.percentile(self.singlePSFAmps, 70)):# and candidate.amp > 800):
                        # don't accept holes that are within a PSF width of a previously fitted peak
                        cleared = True
                        for psf in self.allAcceptedPSFs:
                            dist = np.sqrt((candidate.row - psf.row)**2 + (candidate.col - psf.col)**2)
                            #print(dist)
                            if dist < 3:
                                print('too close')
                                cleared = False
                                break
                        if cleared:
                            singlePSFs.append(candidate)
                            candidate.paint(model_im)
                            self.allAcceptedPSFs.append(candidate)
        return (model_im, singlePSFs, overlapping_peaks)
                

    def subtractBackground(self):
        bg = np.copy(self.bg)
        im = np.copy(self.im)

        mean_bg = np.mean(bg)
        #stdv_bg = np.std(bg)

        # subtract background
        bgsub = mean_bg * self.meanBgPropThresh
        im[im < bgsub] = 0
        im[im != 0] -= np.floor(bgsub).astype(np.int32)

        bg[bg < bgsub] = 0
        bg -= np.floor(bgsub).astype(np.int32)
        im[im < bg] = 0
        bg[im == 0] = 0
        im -= bg
        self.imbgsub = im
        im_filt = np.copy(im)


        # create mask of opening to remove speckles
        mask = np.zeros(np.shape(im_filt))
        mask[im_filt > 0] = 1
        selem = skimage.morphology.disk
        mask = skimage.morphology.opening(mask, selem(2))
        im_filt[mask == 0] = 0
        self.im_filt = im_filt

    def savePSFs(self, filename):
        '''
        Saves PSfs to pickle file
        :param filename: name of file to save to
        :return:
        '''
        with open(filename, 'wb') as f:
            pickle.dump(self.PSFs, f)


    '''
    plots Single PSFs, if peaks_0 set True, will mark peaks successfully fit
        on image by setting them to zero
    '''
    def plotSinglePSFs(self, logScale = False):
        fig, ax = plt.subplots(1,5)
        #fig, ax = plt.subplots(1,4)
        
        ax[0].set_title('Original Image')
        if logScale:
            ax[0].imshow(np.log2(self.im))#, cmap = 'Greys')
            ax[1].imshow(np.log2(self.imbgsub+1))
            ax[2].imshow(np.log2(self.imbgsub+1))
            ax[3].imshow(np.log2(self.im_filt))
            ax[4].imshow(np.log2(self.singlePSFsMod+1))
        else:
            ax[0].imshow(self.im)
            ax[1].imshow(self.imbgsub)
            ax[2].imshow(self.imbgsub)
            ax[3].imshow(self.im_filt)
            ax[4].imshow(self.singlePSFsMod)

        ax[1].set_title('Background Subtracted Image With Local Maxima and Accepted Peaks')
            #pkim = np.copy(self.im)
        frows = []
        fcols  = []
        for psf in self.singlePSFs:
            frows.append(psf.row)
            fcols.append(psf.col)
            #pkim[int(np.round(psf.row)), int(np.round(psf.col))] = 0
        prows = []
        pcols = []
        for peak in self.loc_maxes:
            prows.append(peak[0])
            pcols.append(peak[1])
        orows = []
        ocols = []
        for psf in self.candidate_overlapping_peaks:
            orows.append(psf.row)
            ocols.append(psf.col)
        ax[1].plot(pcols, prows, 'rx', label='Rejected Local Maxima')
        ax[1].plot(fcols, frows, 'w.', label='Accepted PSFs')
        ax[1].plot(ocols, orows, 'k+', label='Potentialy Overlapping PSFs')


        ax[2].set_title('Background Subtracted Image')
        ax[3].set_title('Background Subtraction, Morphological Opening')


        ax[4].set_title('Accepted PSF Plot')

        
        '''
        canvas2 = np.zeros(np.shape(self.im))
        for psf in self.overlapping_peaks: psf.paint(canvas2)
        ax[2].set_title('Significant PSFs')
        ax[2].imshow(canvas2)
        
        ax[3].set_title('Difference')
        ax[3].imshow(self.singlePSFsMod + canvas2 - self.im)
        '''
        plt.show()
        
    '''
    fitAmpMix models the whole image including all the PSFs added to the singlePSFs list
        to get amplitudes for all of them with the same background.
        
    Inputs:
        regionBounds - tuple region bounds of form: (row_min, row_max, col_min, col_max)
                used directly for indexing, so the indice of the max bounds
                will not be included.
        additionalPSFs - a list of more PSFs to add to the model. Can be used to add holes back in.
        save - Boolean, whether  Saves the fit object to self.ampMixMod and a
            plot of it to self.ampMixModIm
    Output: 
        ampMixModFit - the mixmodel fit object
        ampMixModIm - image of the mixmodel predictions
    
    '''
    def fitAmpMixMod(self, regionBounds = None, additionalPSFs = [], save = False):
        #if additionalPSFs:
        PSFs = self.singlePSFs + self.candidate_overlapping_peaks + additionalPSFs
        #else:
            #PSFs = self.singlePSFs + self.overlapping_peaks
        
        if regionBounds:
            [row_min, row_max, col_min, col_max] = regionBounds
            region = self.im[row_min:row_max, col_min:col_max]
            
            PSFs = [psf for psf in PSFs if psf.row >= row_min and psf.row <= row_max and psf.col >= col_min and psf.col <= col_max]
        else:
            row_min = 0
            col_min = 0
            region = self.im
                
        ampMixModIm = np.zeros(np.shape(region))
        
        #estimate starting values for the background and amplitude of each PSF
        nPSFs = len(PSFs)
        amps = np.zeros(nPSFs)
        backgrounds = np.zeros(nPSFs)
        #row_centers = np.zeros(nPSFs)
        #col_centers = np.zeros(nPSFs)
        for i, psf in enumerate(PSFs):
            backgrounds[i] = psf.b
            amps[i] = psf.amp
            #row_centers[i] = psf.row
            #col_centers[i] = psf.col
        b_ave = np.mean(backgrounds)
        
        # ready model starting parameters and bounds
        amps_start = amps + backgrounds - b_ave
        amps_start[amps_start < 0] = 0
        params_start = tuple(np.insert(amps_start, 0, b_ave))
        lbounds = tuple([np.min(backgrounds)] + [0]*nPSFs)
        ubounds = amps+backgrounds
        ubounds = tuple(np.insert(ubounds, 0, np.max(region)))
        bounds = (lbounds, ubounds)
        
        def mixAmpResiduals(params):
            canvas = np.zeros(np.shape(region))
            canvas += params[0]
            for i, psf in enumerate(PSFs):
                psf.paint(canvas, params[i+1], row_min, col_min)
                
            residuals = canvas - region
            
            return residuals.flatten()
        
        ampMixMod = least_squares(mixAmpResiduals, params_start, bounds = bounds)
        
        ampMixModIm += ampMixMod.x[0]
        for i, psf in enumerate(PSFs):
            psf.paint(ampMixModIm, ampMixMod.x[i+1], row_min, col_min)
            
        if save:
            self.ampMixMod = ampMixMod
            self.ampMixModIm = ampMixModIm
            
        return ampMixMod, ampMixModIm
            
    '''
    plot ampMixMod
    '''
    def plotAmpMixMod(self, peaks_0 = False):
        fig, ax = plt.subplots(1,4)
        
        ax[0].set_title('Image')
        if peaks_0:
            pkim = np.copy(self.im)
            allholes = self.wfholes + self.olholes
            allholes = [hole for iterholes in allholes for hole in iterholes]
            psfs = self.singlePSFs + self.candidate_overlapping_peaks #+ allholes
            for psf in psfs:
                pkim[int(np.round(psf.row)), int(np.round(psf.col))] = 0
            ax[0].imshow(pkim)
        else:
            ax[0].imshow(self.im)
        
        ax[1].set_title('Single PSFs Model')
        ax[1].imshow(self.singlePSFsMod)
        
        ax[2].set_title('Amplitude Mixture Model')
        ax[2].imshow(self.ampMixModIm)
        
        ax[3].set_title('Difference')
        ax[3].imshow(self.ampMixModIm - self.im)
        
        plt.show()

    '''
    selectOverlapping - runs mix models to allow peaks marked as overlapping to move farther and find their fit.
    '''
    def selectOverlapping(self):
        # first, find regions for each mix model by dilating around candidate 'overlapping' peaks, and accepted peaks
        regions = np.zeros(np.shape(self.im))

        for psf in self.singlePSFs:
            psf.paint(regions)

        for peak in self.candidate_overlapping_peaks:
            peak.paint(regions)

        regions[regions > 0] = 1


        # find each connected region that contains a candiate overlapping peak
        labeled_regions, n_regions = skimage.measure.label(regions)

        # for each connected region that contains a candidate overlapping peak, set up a mix model
        op_candidates = copy.deepcopy(self.candidate_overlapping_peaks)

        # for each region
        for i in range(n_regions+1):
            # look for overlapping peak candidates in region
            cand = 0
            region_cands = []
            while cand < len(op_candidates):
                row = op_candidates[cand].row
                col = op_candidates[cand].col
                if labeled_regions[row,col] == i:
                    region_cands.append(op_candidates.pop(cand))
                else:
                    cand +=1

            # if at least 1 overlapping peak candidate was found, set up a mix model
            if len(region_cands) > 0:
                accepted_psfs = []
                # find accepted peaks in region
                for psf in self.singlePSFs:
                    if labeled_regions[psf.row, psf,col] == i:
                        accepted_psfs.append(psf)

                # set up mix model bounds
                for psf in accepted_psfs:
                    width = psf.width
                    row = psf.row
                    col = psf.col
                    amp = psf.amp




        # allow accepted peaks to vary amplitude and width, but not position. Allow canidate overlapping
        # peaks to vary amplitude, width, and position by a few pixels.

        return
    '''
    findHoles looks for places in the image where the individually fit PSFs in an amplitude
        mixture model underpredicts intensity in the image. To do this, 
    '''
    def findHoles(self):
        # fit a tiling of small amplitude mixture models so that holes in the
        # central region of each fit will be very well defined and that each individual
        # model will have a managable runtime.
        # Images are 2048x2048 pixels, so let's use tiles of 64x64, interleafed 
        # so a new tile starts every 32 pixels
        
        # initialize dictionary for storing information on holes keys are the bounds
        # of the tile (row_min, col_min, col_max)
        tiles = {}
        
        tilewidth= 32
        
        for pAcross in range(64):
            for pDown in range(64):
                continue
        
        return

    '''
    findHolesInTile looks for holes that an amplitude mixture model of fit
    single PSFs underpredicts in a subregion 'tile' of the whole image. This
    the number number of PSFs and therefore parameters that must be fit in the 
    mixture model, and allows us to search the entire image more quickly than
    fitting it all at once, and better accounts for local variations in 
    background.
    
    Inputs:
        regionBounds - tuple region bounds of form: (row_min, row_max, col_min, col_max)
                used directly for indexing, so the indice of the max bounds
                will not be included.
    Returns:
        wellFittingHoles - list of lists containing well fitting holes found in each iteration
        overlappingHoles - list of lists containing holes that hit fitting bounds in each iteration
    '''
    def findHolesInTile(self, regionBounds):
        
        wellFittingHoles = []
        overlappingHoles = []
        #iteration = 1
        
        
        #while(True):
        for iteration in range(1):
            print('iteration:', iteration+1)
            
            allholes = wellFittingHoles + overlappingHoles
            allholes = [hole for iterholes in allholes for hole in iterholes]
            # fit mixMod to tile
            mixmodfit, mixmodim = self.fitAmpMixMod(regionBounds = regionBounds, additionalPSFs = allholes)
            self.mixModsWHoles.append((mixmodfit, mixmodim))
            [rmin, rmax, cmin, cmax] = regionBounds
            
            #find difference between the model and the original image
            difference = self.im[rmin:rmax, cmin:cmax] - self.ampMixModIm
            
            # anything that is less than 0 is well accounted for by the model, so 
            # let's set it to 0
            difference[difference < 0] = 0
            
            
            model_im, singlePSFs, overlapping_peaks = self.fitSinglePSFs(difference)
            
            if len(singlePSFs) == 0 and len(overlapping_peaks) == 0:
                break
            
            wellFittingHoles.append(singlePSFs)
            overlappingHoles.append(overlapping_peaks)
            print('Total Holes Found:', len(singlePSFs) + len(overlapping_peaks))
            iteration += 1
            # todo: adjust overlapping peaks coordinates with mix model before proceeding
        
        self.wfholes += wellFittingHoles
        self.olholes += overlappingHoles
        return wellFittingHoles, overlappingHoles
        
    '''
    '''
    def plotHoles(self, holefit = 0):
        fig, ax = plt.subplots(1,4)
        
        ax[0].set_title('Image with accepted peaks')
        pkim = np.copy(self.im)
        psfs = self.singlePSFs + self.candidate_overlapping_peaks
        for psf in psfs:
            pkim[int(np.round(psf.row)), int(np.round(psf.col))] = 0
        ax[0].imshow(pkim)
        
        ax[1].set_title('Image with accepted peaks and holes')
        pkim = np.copy(self.im)
            #allholes = self.wfholes + self.olholes
            #allholes = [hole for iterholes in allholes for hole in iterholes]
            #psfs = self.singlePSFs + self.overlapping_peaks + allholes
        for psf in self.allAcceptedPSFs:
            pkim[int(np.round(psf.row)), int(np.round(psf.col))] = 0
        ax[1].imshow(pkim)
        
        ax[2].set_title('Mixmod Difference before holes')
        ax[2].imshow(self.ampMixModIm - self.im)
        
        ax[3].set_title('Mixmod Difference with holes')
        ax[3].imshow(self.mixModsWHoles[holefit][1] - self.im)
        
        plt.show()
        
    '''
    '''
    def plotStages(self):
        
        fig, ax = plt.subplots(1,4)
        
        ax[0].set_title('Image')
        ax[0].imshow(self.im)
        
        imlocmax = np.copy(self.im)
        for lmax in self.loc_maxes: imlocmax[lmax[0], lmax[1]] = 0
        ax[1].set_title('Local Maxima')
        ax[1].imshow(imlocmax)
        
        pkim = np.copy(self.im)
        psfs = self.singlePSFs + self.candidate_overlapping_peaks
        for psf in psfs:
            pkim[int(np.round(psf.row)), int(np.round(psf.col))] = 0
        ax[2].set_title('Accepted Peaks')
        ax[2].imshow(pkim)
        
        imallpsf = np.copy(self.im)
        for psf in self.allAcceptedPSFs:
            imallpsf[int(np.round(psf.row)), int(np.round(psf.col))] = 0
        ax[3].set_title('Accepted Peaks and holes')
        ax[3].imshow(imallpsf)
        
        plt.show()
    
    '''
    fitMixModels identifies PSFs flagged for the mix model by fitSinglePSFs
    '''
    '''
    def fitMixModels(self):
        # Plot each peak marked for the mixture model and dilate by 3 PSF widths
        clusters = np.zeros(np.shape(self.im))
        for psf in self.overlapping_peaks: clusters[psf.peak_row, psf.peak_col] = 255
        dilation_radius = int(self.PSFwidth*3)
        elem = cv.getStructuringElement(cv.MORPH_ELLIPSE,(dilation_radius, dilation_radius))
        clusters = cv.dilate(clusters, elem)
        cv.imshow('clusters', clusters)
        cv.waitKey(0)
                   
        return
    '''
'''
The fitPSF class fits a pointspread function in an image
'''        
class fitPSF:
    # class member function PSF and class data member PSF_width are set
    # by the constructor of the MixMod class. Therefore, this class can only
    
    PSF = gaussian2D
    width = None
    multiplier = 2
    
    '''
    fitPSF fits the PSF in its constructor.
    Inputs:
        peak_row - the row of the peak to fit
        peak_col - the column of the beak to fit
        im - the image to fit to
    '''
    def __init__(self, im, peak_row, peak_col):
        
        self.peak_row = peak_row
        self.peak_col = peak_col
        self.nrows, self.ncols = np.shape(im)
        
        # we will the bounds to the immediate region of the image that we
        # can draw the PSF later
        #row_min = int(np.floor(peak_row - self.width*4))
        #col_min = int(np.floor(peak_col - self.width*4))
        #row_max = int(np.ceil(peak_row + self.width*4))
        #col_max = int(np.ceil(peak_col + self.width*4))
        multiplier = fitPSF.multiplier
        [row_min, row_max, col_min, col_max, ntrimmed1] = self.getRegionBounds(peak_row, peak_col, multiplier)
        #print('Base Region:')
        #print([row_min, row_max, col_min, col_max, ntrimmed1])
        #print('------------------')
        
        region = im[row_min:(row_max + 1), col_min:(col_max+1)]
        base_size = int((2*np.ceil(fitPSF.width_guess*multiplier)+2)**2) #np.size(region) + ntrimmed1

        b_start = np.min(region)
        b_lower = 0
        b_upper = np.max(region)

        width_start = self.width_guess
        width_lower = width_start / 2
        width_upper = width_start * 2
        
        r_lower = peak_row-1
        r_upper = peak_row+1
        
        c_lower = peak_col-1
        c_upper = peak_col+1
        
        amp_start = b_upper - b_start
        amp_lower = 0
        amp_upper = b_upper
        
        start = (b_start,width_start,  peak_row, peak_col, amp_start)
        lbounds = (b_lower,width_lower,  r_lower, c_lower, amp_lower)
        ubounds = (b_upper, width_upper, r_upper, c_upper, amp_upper)
        bounds = (lbounds, ubounds)
        
        # define function to calculate residuals to pass to fitting function
        def Gauss2DResiduals(params):
            #find subregion to fit to
            #row_fit_min = int(np.floor(params[1] - self.width*1.3))
            #col_fit_min = int(np.floor(params[2] - self.width*1.3))
            #row_fit_max = int(np.ceil(params[1] + self.width*1.3))
            #col_fit_max = int(np.ceil(params[2] + self.width*1.3))
            
            # get baseline  region bounds to ensure that residuals function always
            # even when the region moves around the image and gets cut off the edge.
            regionBounds = self.getRegionBounds(params[2], params[3], multiplier)
            [row_fit_min, row_fit_max, col_fit_min, col_fit_max, ntrimmed] = regionBounds
            
            grid = np.mgrid[row_fit_min:(row_fit_max + 1), col_fit_min:(col_fit_max+1)]
            region = im[row_fit_min:(row_fit_max + 1), col_fit_min:(col_fit_max+1)]
            
            #calculate residuals
            residuals = region - fitPSF.PSF(grid, *params)

            residuals = residuals.flatten()
            #print([row_fit_min, row_fit_max, col_fit_min, col_fit_max, ntrimmed])
            #print('center:', params[1], params[2])
            #print('base_size:', base_size, 'residual size:', np.size(residuals))
            
            if base_size > np.size(residuals):#ntrimmed > 0:
                residuals = np.concatenate((residuals, [np.mean(residuals)]*(base_size-np.size(residuals))))
            
            return residuals
        
        self.fit = least_squares(Gauss2DResiduals, start, bounds = bounds)
        self.b, self.width, self.row, self.col, self.amp = self.fit.x
        self.RMSE = np.sqrt(np.mean(self.fit.fun**2))
        
        # checks whether the the fit amplitude of the PSF is larger than the 
        # RMSE, or if fit_amp is 0
        #if self.RMSE >= 3*self.amp or self.amp == 0:
        #    self.amplitude_significant = False
        #else:
        #    self.amplitude_significant =  True
        
    '''
    getRegionBounds returns the bound of a region of an image in which to 
        consider a PSF. It find the bounds by multiplying the system width by
        the given multiplier parameter, adding/subtracting this width to the
        center coordinates, rounding such that they are farther away, and finally
        checking to see if these leave the bounds of the image. The number of
        trimmed elements is returned.
    Inputs:
        center_row - the center row coordinate of the region to consider
        center_col - the center column coordinate of the region to consider
        multiplier - defines the size of the region by multiplying the width parameter
                    of the PSF to give how far the region should stretch from the
                    center coordinates.
    Outputs:
        row_min - the minimum row of the region
        row_max - the maximum row of the region
        col_min - the minimum row of the region
        col_max - the maximum row of the region
        uncut_size - the size of the region before trimming according to the 
                    bounds of the actual image.
        nrows - if trying to find a region of a tile, the number of rows in the tile
        ncols - if trying to find a region of a tile, the number of columns in the tile
    '''
    def getRegionBounds(self, center_row, center_col, multiplier, nrows = None, ncols = None):
        if self.width:
            width = self.width
        else:
            width = self.width_guess
        row_min = int(np.floor(center_row - width*multiplier))
        col_min = int(np.floor(center_col - width*multiplier))
        row_max = int(np.ceil(center_row + width*multiplier))
        col_max = int(np.ceil(center_col + width*multiplier))
        uncut_size = (row_max - row_min)*(col_max - col_min)
        
        if not nrows and not ncols:
            nrows = self.nrows
            ncols = self.ncols
        
        if row_min < 0:
            row_min = 0
        if col_min < 0:
            col_min = 0
        if row_max >= nrows: #self.nrows:
            row_max = nrows - 1 #self.nrows - 1
        if col_max >= ncols: #self.ncols:
            col_max = ncols - 1 #self.ncols - 1
            
            
        ntrimmed = uncut_size - (row_max - row_min)*(col_max - col_min)
            
        return [row_min, row_max, col_min, col_max, ntrimmed]
        
    
    '''
    paint paints the PSF onto the past canvas, a 2D numpy array. the canvas
        should be the same size as the original image passed to the PSF object.
    
    Inputs:
        alt_amp - alternative amplitude with which to paint PSF other than its fitted
                PSF. Useful for fitting follow up models.
        row_shift - how many rows (positive) to shift the center by (for painting tiles)
        col_shift - how many columns (positive) to shift the center by (for painting tiles)
    '''
    def paint(self, canvas, alt_amp = None, row_shift = 0, col_shift = 0):
        #grid = np.mgrid[self.row_min:(self.row_max+1), self.col_min:(self.col_max+1)]
        nrows, ncols = np.shape(canvas)
        [row_min, row_max, col_min, col_max, ntrimmed] = self.getRegionBounds(self.row- row_shift, self.col - col_shift, fitPSF.multiplier, nrows, ncols)
        try:
            grid = np.mgrid[row_min:(row_max + 1), col_min:(col_max+1)]
        except:
            raise(Exception('Uh oh!'))
        b, width, row_center, col_center, amp = self.fit.x
        if alt_amp:
            amp = alt_amp
        PSF_grid = gaussian2D(grid, 0, width, row_center, col_center, amp)
        #row_min -= row_shift
        #col_min -= col_shift
        row_max_ind = row_max + 1 #- row_shift
        col_max_ind = col_max + 1 #- col_shift
        #print(row_min, row_max, col_min, col_max)
        try:
            canvas[row_min:row_max_ind, col_min:col_max_ind] += PSF_grid
        except:
            raise(Exception('Uh oh!'))
            
        return canvas
            
    
        