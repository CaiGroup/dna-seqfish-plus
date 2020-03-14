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
from multiprocessing import Pool


preformatted_dir = 'pre_formated'
ref_fname = 'ref-points-ForBeadAlignment-pos%d-radial3d-finalBeads-raw-intensity-2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH.csv'
ro_fname = 'hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d-2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH.csv'

maxXYMatchError = 0.5
maxXYSearchError = 3
maxZMatchError = 1
maxZSearchError = 5
minMatch = 5
stopRefProp = 0.2
nLongestCheck = 3000

def alignPos(pos):
    os.chdir('pos%d' % pos)
    
    align(maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError, 
          minMatch, stopRefProp, nLongestCheck, ref_name, hyb_name, pos, 'final')

#home_dir = '/Users/jonathanwhite/Box/Jonathan/2019-09-09-brain-rep2-2-DNAFISH/phase2_many_hybs'
home_dir = 'C:\\Users\\jonat\Box\\Jonathan\\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\\20191114'
#os.chdir('..'
#home_dir = os.getcwd()
os.chdir(home_dir)

maxPosition = 0


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

#if __name__ == '__main__':
    #pool = Pool(processes = 1)
    #alignment = pool.map(alignPos, range(1))
alignPos(0)
#for pos in range(maxPosition + 1):
    #os.chdir('pos%d' % pos)
    
    #align(maxXYError, maxZError, minMatch, stopRefProp, nLongestCheck, ref_name, hyb_name, pos)
        

# assemble queue of positions to align

# create processes to align queue




