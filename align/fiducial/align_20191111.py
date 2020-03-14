# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:29:26 2019

@author: jonat
"""

import sys
sys.path.append('C:\\Users\\jonat\\PycharmProjects\\seqfishimageprocessing')
sys.path.append('/Users/jonathanwhite/PycharmProjects/seqfishimageprocessing')

import os
from dotAligner3D import dotAligner


os.chdir('..')

def align(maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError, 
          minMatch, stopRefProp, nLongestCheck, ref_fname, hyb_fname, position):
    
    saveOffsetsName = 'offsets_ch%d_mXYsEr%s_mZsEr%d_mXYmEr%s_mZmEr%d_minMatch%d_minMatch_nLongCheck%d_' + hyb_fname
    smaxXYSearchError = str(maxXYSearchError).replace('.', 'p')
    smaxXYMatchError = str(maxXYMatchError).replace('.', 'p')
    #sstopRefProp = str(stopRefProp).replace('.', 'p')
    ref_fname = 'ch_%d_' + ref_fname
    hyb_fname = 'ch_%d_' + hyb_fname
    
    for ch in (1,2,3):
        os.chdir('ch%d' % ch)
        os.chdir('ch%d_points' % ch)
        
        dal = dotAligner()
        
        #load points
        dal.loadRefPSFs(ref_fname % (ch, position))
        dal.loadPSFs(hyb_fname % (ch, position))
        
        #return to channel directory
        os.chdir('..')
        
        #find offsets
        #dal.find_offsets_from_ref(maxXYError=maxXYError, maxZError=maxZError, minMatches=minMatch, stopRefProp=stopRefProp, nBrightestCheck=nLongestCheck)
        dal.find_offsets_from_ref(maxXYSearchError=maxXYSearchError, maxZSearchError=maxZSearchError, maxXYMatchError=maxXYMatchError, maxZMatchError=maxZMatchError, minMatches=minMatch, stopRefProp=stopRefProp, nLongestCheck=nLongestCheck)
        #dal.find_offsets_from_ref_par(maxXYSearchError=maxXYSearchError, maxZSearchError=maxZSearchError, 
            #maxXYMatchError=maxXYMatchError, maxZMatchError=maxZMatchError, minMatches=minMatch,
            #stopRefProp=stopRefProp, nLongestCheck=nLongestCheck)
        
        # set up directory to save points
        offsets_dir_name = 'ch%d_offsets' % ch
        if offsets_dir_name not in os.listdir():
            os.mkdir(offsets_dir_name)
        os.chdir(offsets_dir_name)
        
        dal.saveOffsets(saveOffsetsName % (ch, smaxXYSearchError, maxZSearchError,
            smaxXYMatchError, maxZMatchError, minMatch, nLongestCheck, position))
        #dal.saveMatchesDF(saveMatchesDFName % (ch, maxXYError, maxZError, minMatch, nLongestCheck))
        
        #dal.saveAlignedPSFs(saveAlignedName % (ch, maxXYError, maxZError, minMatch))
        
        #return to positions directory
        os.chdir('..')
        os.chdir('..')
        
        