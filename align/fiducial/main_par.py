#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 15:25:33 2019

@author: jonat


Aligns all hybs in all positions
"""

import os 
import sys #, getopt
from reformat_reference import reformat_reference
from reformat_read_out import reformat_read_out
from align import align
from multiprocessing import Pool

#def main(argv)
preformatted_dir = 'pre_formated'
#home_dir = sys.arg[1]
#pos = sys.arg[2]
#ref_fname = sys.arg[3]
#ro_fname = sys.arg[4]
home_dir = 'I:\\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\\points'
ref_fname = 'ref-points-ForBeadAlignment-pos%d-radial3d-final_fiducial_markers-raw-intensity.csv'
ro_fname = 'hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped.csv'
pos = 1

#home_dir = 'C:\\Users\\jonat\Box\\Jonathan\\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\\20191114'
#ref_fname = 'ref-points-ForBeadAlignment-pos%d-radial3d-initialBeads-raw-intensity-2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH.csv'
#ro_fname = 'hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d-2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH.csv'

# The Search Error Dictates 
maxXYSearchError = 1.5
maxZSearchError = 4

maxZMatchError = 1
maxXYMatchError = 0.3

minMatch = 10
stopRefProp = 0.2
nLongestCheck = 5000

def alignPos(pos):
    os.chdir('pos%d' % pos)
    
    align(maxXYSearchError, maxZSearchError, maxXYMatchError, maxZMatchError, 
          minMatch, stopRefProp, nLongestCheck, ref_name, hyb_name, pos, 'init')


#os.chdir('..'
#home_dir = os.getcwd()
os.chdir(home_dir)

#maxPosition = 0


# for each  position 
#for pos in range(maxPosition+1):
    # format reference points
os.chdir(home_dir)
reformat_reference(ref_fname, pos, preformatted_dir)

# format read out points
os.chdir(home_dir)
reformat_read_out(ro_fname, pos, preformatted_dir)

# do alignments
os.chdir(home_dir)
if 'positions' not in os.listdir():
    os.mkdir('positions')
os.chdir('positions')

ref_name = ref_fname
hyb_name = ro_fname

#if __name__ == '__main__':
    #pool = Pool(processes = 1)
    #alignment = pool.map(alignPos, range(1))
alignPos(pos)
#for pos in range(maxPosition + 1):
    #os.chdir('pos%d' % pos)
    
    #align(maxXYError, maxZError, minMatch, stopRefProp, nLongestCheck, ref_name, hyb_name, pos)
        

# assemble queue of positions to align

# create processes to align queue

#if __name__ == "__main__":
#   main(sys.argv[1:])



