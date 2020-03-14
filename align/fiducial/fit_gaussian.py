# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 20:28:18 2019

@author: jonat
"""
import numpy as np
from scipy.optimize import least_squares

# give bounds of maximum and minimum allowable sigma of the PSF as a global variable:
sigma_min = 1 # 1 pixel
sigma_max = 2 # 2 pixels

'''
gaussian2D - Draws a 2D Gaussian curve on a grid

Inputs:
    grid - numpy array giving x and y values a grid given by numpy.mgrid
    center_col - location along the grid to center the curve on the axis of columns
    center_row - location along the grid to center the curve on the axis of rows
    amp - amplitude of the curve
    sigma - the width of the curve
    b - the base level of the curve
Returns:
    2D grid of the same dimensions as the passed grid, with the 2D Gaussian drawn on it.    
'''
def gaussian2D(grid, b, sigma,  center_row, center_col, amp):
    row_grid = grid[0]-center_row
    col_grid = grid[1]-center_col
    rsqgrid = col_grid**2 + row_grid**2
    return amp*np.exp(-rsqgrid/sigma**2)+ b

'''
fitGaussian2D fits a single gaussian to the given image using
 the Levenberg-Marquart algorithm as implemented in scipy.optimize.least_squares

Inputs:
    im - the image containing pixels of measured fluorescence intensity values
    params_start - tuple containing the starting values for the fit parameters
    (bguess, center_col_guess, center_row_guess, amp_guess):
        b_guess - the starting background
        sigma - the sigma parameter of the gaussian approximation to the system PSF
        center_row_guess - the starting value for the gaussian's center along the axis of rows
        center_col_guess - the starting value for the gaussian's center along the axis of columns
        amp_guess - the starting value for gaussian amplitude
    params_bounds -  optional: tuple of two tuples containing the bounds for the search 
        space for each fitting parameter. The first subtuple contains lower bound,
        the second contains upper bounds. If not given, default bounds are set.
        See scipy.optimize.leastsquares
 
Returns: list containing the fit parameters and other relevant information 
            returned by scipy.optimize.least_squares       
 '''   
def fitGaussian2D(im, params_start, params_bounds = None):

    # if bounds for fitting parameters not given,
    if not params_bounds:
        # then set bounds
        params_bounds = Gaussian2DBounds(params_start)
        
    # get grid for calculating gaussian curve
    nrows, ncols = np.shape(im)
    grid = np.mgrid[0:nrows, 0:ncols]
    
    # define function to calculate residuals to pass to fitting function
    def Gauss2DResiduals(params):
        residuals = im - gaussian2D(grid, *params)
        return residuals.flatten()
    
    return least_squares(Gauss2DResiduals, params_start, params_bounds)
        
'''
GaussianMixModel2D draws a Gaussian mixture model on a 2D grid
    Inputs:
        grid - grid giving row and column coordinates as returned by numpy.mgrid
        sigma - the width of the Gaussian approximation the the system PSF
        params - tuple of size at least 2 containing the parameters for the 
        mixture model. 0th element contains the base background level of the 
        curve, b. The 1st through nth elements are tuples containing the
        parameters for 2D gaussian curves: (center_row, center_col, amp).
    Returns:
        a 2D numpy array of the dimensions specified by the passed grid with the
        mixture model drawn on it.
'''
def GaussianMixModel(grid, params):
    throwaway, nrows, ncols = np.shape(grid)
    out = np.zeros([nrows, ncols])
    
    for i, gparams in enumerate(params):
        # draw each gaussian on the output grid
        if i > 0:
            out += gaussian2D(grid, 0, *gparams)
            
        # if at the 0th parameter (the background), add it then break
        elif i == 0:
            out += gparams
        
    return out
            
    
'''
fitGaussianMixMod fits a gaussian mixture model onto the given image based on
the given starting parameters
    im - the image containing pixels of measured fluorescence intensity values
    sigma - the known sigma parameter of the gaussian approximation to the system PSF
    params_start - n tuple containing starting values for the mix model 
        parameters. The 0th position contains the starting value for the 
        baseline/background value. The 1st through nth positions contain tuples
        containing the the fit parameters of each Gaussian to be fit.
        For example: params_start = (b, g1, g2, ..., gn)
        where gi = (center_row_guess, center_col_guess, amp_guess)
        and
            center_row_guess - the starting value for the gaussian's center along the axis of rows
            center_col_guess - the starting value for the gaussian's center along the axis of columns
            amp_guess - the starting value for gaussian amplitude
    params_bounds -  n tuple containing the bounds for the search 
        space for each fitting parameter of each gaussian to fit. The 0th
        position contains a tuple containing the lower and upper bounds for the
        base/background parameter, b, search space. The 1st through nth positions
        contain 2 tuples denoting the bounds for the search space
        of each Gaussian in the mixture mode's parameters. The first tuple in the 
        2 tuple has the lower bounds of the gaussian fit parameters, and the 
        second tuple of the 2-tuple contains the upper bounds of the gaussian's fit parameters.
        See scipy.optimize.leastsquares
'''
def fitGaussianMixMod(im, sigma, params_start, bounds = None):
    print('fitting Mixture model of %d Gaussians' % (len(params_start)-1))
    
    if not bounds:
        # make default bounds
        bounds = gMixModDefBounds(params_start)
    
    # flatten params so that it is ready to feed into scipy.optimize.least_squares
    params_start_flat = []
    for i, fparams in enumerate(params_start):
        #if at the 1st through nth position, add the gaussian parameters in a tuple
        if i > 0:
            for param in fparams:
                params_start_flat.append(param)
        # if at the 0th position, add the background as float
        elif i == 0:
            params_start_flat.append(fparams)
            
    
    # flatten bounds so it is ready to feed into scipy.optimize.least_squares
    lbound = [bounds[0][0]]
    ubound = [bounds[0][1]]
    #for each function to be bounded
    for fbounds in bounds[1:]:
        #get the bounds of that function
        for i, fbound in enumerate(fbounds):
            # for each parameter in that function
            for param in fbound:
                # add to flat lower bounds list if it is a lower bound
                if i == 0:
                    lbound.append(param)
                #add to flat upper bounds list if it is an upper bound
                elif i == 1:
                    ubound.append(param)
                    
    #get grid for calculating mix model         
    nrows, ncols = np.shape(im)
    grid = np.mgrid[0:nrows, 0:ncols]
    
    def GaussianMixModelResiduals(params):
        # will need to regroup parameters into nested list before passing into
        # GaussianMixMod
        nested_params = nestMixModParams(params)
        residuals = im - GaussianMixModel(grid, sigma, nested_params)
        return residuals.flatten()
    
    return least_squares(GaussianMixModelResiduals, params_start_flat, bounds = (lbound, ubound))
    
    
def findGaussianMixMod():
    return
    # Should fit at least as many points as there are maxes. If there
    # are more than 16 positive points in the mask for every peak, the model
    # likely needs to include more PSFs than peaks

'''
Gaussian2DBounds sets default bounds for the search space for the fit parameters
of a 2D Gaussian.
    Inputs:
        params_guess - the starting values for the fitting parameters in a tuple:
            params_guess = (b_guess, center_row_guess, center_col_guess, amp_guess)
                b_guess - the guess for the background/baseline of the Gaussian
                center_row_guess - guess for the x coordinate center of the Gaussian
                center_col_guess - guess for tye y coordinate center of the Gaussian
                amp_guess - guess for the amplitude of the Gaussian.
    Returns:
        Two tuple containing two tuples containing the lower bounds of each 
        fit parameter in the first tuple and the upper bounds in the second tuple:
            ((b_lower, xlower, ylower, amp_lower),
                        (b_upper, xupper, yupper, amp_upper))
'''
def Gaussian2DBounds(params_guess):
    # unpack parameters
        b_guess, center_row_guess, center_col_guess, amp_guess = params_guess
        
        # set bounds
        b_lower = 0
        b_upper = amp_guess+b_guess
        lgbounds, ugbounds = _g2DBound(center_row_guess, center_col_guess, amp_guess)
        
        params_bounds = (tuple([b_lower] + lgbounds),
                        tuple([b_upper] + ugbounds))
        return params_bounds

'''
_g2DBounnds is a helper function that returns default upper and lower bounds
    for fitting parameters of a gaussian (except the background/baseline) based
    on their starting guess.
'''
def _g2DBounds(center_row_guess, center_col_guess, amp_guess):
    rlower = center_row_guess-1
    rupper = center_row_guess+1
    clower = center_col_guess-1
    cupper = center_col_guess+1
    amp_lower = 0
    amp_upper = amp_guess*1.3
    return [[rlower, clower, amp_lower], [rupper, cupper, amp_upper]]

'''
gMixModDefBounds gives default bounds on the search space for Gaussian mix
    model parameters based on their starting guesses.
Input:
    params_guess - starting guesses for the mix model's parameters given in the
        same for as should be given to the 'params' in GaussianMixModel2D or 
        'params_start' in fitGaussianMixMod.
Returns: bounds on the search space for gaussian mix model parameters in the 
    form accepted by fitGaussianMixMod
'''
def gMixModDefBounds(params_guess):
    bounds = [(0, 1000)]
    for i, fparams in enumerate(params_guess):
        if i != 0:
            flbounds, fubounds = _g2DBounds(*fparams)
            bounds.append((tuple(flbounds), tuple(fubounds)))
    return tuple(bounds)

'''
Converts flat array of mix model parameters to nested array form.
'''
def nestMixModParams(flat_params):
    if len(flat_params)%3 != 1:
        print('flat params:', flat_params)
    nested_params = []
    #add the base line
    #nested_params.append(flat_params.pop(0))
    nested_params.append(flat_params[0])
    #now add the gaussians
    i = 1
    while len(flat_params) > i:
        g = []
        for j in range(3):
            #print(i, len(flat_params))
            g.append(flat_params[i])
            i += 1
        nested_params.append(tuple(g))
    return tuple(nested_params)