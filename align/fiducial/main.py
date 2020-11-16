# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:25:33 2019

@author: jonat


Aligns all hybs in all positions
"""

import os
import sys
from reformat_reference import reformat_reference
from reformat_read_out import reformat_read_out
from align import align


preformatted_dir = 'pre_formated_aligned2hyb0'
ref_fname = 'ref-points-ForBeadAlignment-pos%d-radial3d-finalBeads-raw-intensity.csv'
ro_fname = 'hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d.csv'

maxXYError = 1.3
maxZError = 3
minMatch = 10
stopRefProp = 0.1
nLongestCheck = 100


home_dir = '/Users/jonathanwhite/Box/Jonathan/2019-09-09-rep2-2-DNAFISH/phase2_many_hybs'
#os.chdir('..'
#home_dir = os.getcwd()
os.chdir(home_dir)

maxPosition = 4


# for each  position 
for pos in range(maxPosition+1):
    # format reference points
    os.chdir(home_dir)
    reformat_reference(ref_fname, pos, preformatted_dir)

    # format read out points
    os.chdir(home_dir)
    reformat_read_out(ro_fname, pos, preformatted_dir)

# do alignments
os.chdir(home_dir)
os.chdir('positions')
ref_name = ref_fname
hyb_name = ro_fname
for pos in range(maxPosition + 1):
    os.chdir('pos%d' % pos)
    
    align(maxXYError, maxZError, minMatch, stopRefProp, nLongestCheck, ref_name, hyb_name, pos)
        

# assemble queue of positions to align

# create processes to align queue




