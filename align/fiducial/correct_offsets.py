# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 14:31:14 2019

@author: jonat
"""

import os
import pandas as pd
import numpy as np


initial_aligned_positions = 'C:\\Users\\jonat\\Box\\Jonathan\\2019-09-09-brain-rep2-2-DNAFISH\\phase2_many_hybs\\alignInitialFiducials\\positions'
final_aligned_positions = 'C:\\Users\\jonat\\Box\\Jonathan\\2019-09-09-brain-rep2-2-DNAFISH\\phase2_many_hybs\\alignFinalFiducials\\positions'
#final_aligned_positions = 'C:\\Users\\jonat\\Box\\Jonathan\\2019-09-09-brain-rep2-2-DNAFISH\\phase2_many_hybs\\positions'

initDir = 'C:\\Users\\jonat\\Box\\Jonathan\\2019-09-09-brain-rep2-2-DNAFISH\\phase2_many_hybs\\alignInitialFiducials'
finalDir = 'C:\\Users\\jonat\\Box\\Jonathan\\2019-09-09-brain-rep2-2-DNAFISH\\phase2_many_hybs\\alignFinalFiducials'

initial_aligned_offsets_0p3 = pd.DataFrame()
initial_aligned_offsets_0p5 = pd.DataFrame()
final_aligned_offsets0p5 = pd.DataFrame()
final_aligned_offsets0p5 = pd.DataFrame()

maxXYError = '0p3'
maxZError = 3
minMatch = 10
stopRefProp = 0.1
nLongestCheck = 300

hyb_fname_final = 'hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d.csv'
#hyb_fname_initial = 
saveOffsetsName = 'offsets_ch%d_mXYEr%s_maxZEr%d_minMatch%d_minMatch_nLongCheck%d_'
fname = 'offsets_%sch%d_mXYEr%s_maxZEr%d_minMatch10_minMatch_nLongCheck1000_hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d.csv'

#tfname1 ='offsets_%sch%d_mXYEr%s_maxZEr%d_minMatch10_minMatch_nLongCheck10000_hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d.csv'
#tfname = 'offsets_%sch%d_mXYEr%s_maxZEr%d_minMatch10_minMatch_nLongCheck1000_hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d.csv'


def findOffsets(directory, baseFname, ref, XYerror, Zerror):
    os.chdir(directory)
    aligned_offsets = pd.DataFrame()
    
    for position in range(5):
        os.chdir('pos%d' % position)
        for channel in (1,2,3):#(1,2,3):
            os.chdir('ch%d' % channel)
            os.chdir('ch%d_offsets' % channel)
            fname = baseFname%(ref, channel, XYerror, Zerror, position)
            poschoffsets = pd.read_csv(fname)
            poschoffsets['ch'] = [channel]*len(poschoffsets)
            poschoffsets['pos'] = [position]*len(poschoffsets)
            aligned_offsets = aligned_offsets.append(poschoffsets)
            os.chdir('..')
            os.chdir('..')
        os.chdir('..')
    aligned_offsets.set_index(['pos','ch'], inplace=True)
    return aligned_offsets

def refineOffsetsAligned(prefOffsets, looser_offsets):
    #prefnp = np.array(prefOffsets)
    #loosnp = np.array(looser_offsets)
    #prefnan = np.logical_or(np.isnan(prefnp), np.equal(prefnp, 0))
    #prefnp[prefnan] = loosnp[prefnan]
    
    lowMatches = prefOffsets.n_offsets_calculated < 10
    prefOffsets.loc[lowMatches] = looser_offsets.loc[lowMatches]
    return prefOffsets 
    

final_0p3_offsets = findOffsets(final_aligned_positions, fname, '', '0p3', 3)
final_0p5_offsets = findOffsets(final_aligned_positions, fname, '', '0p5', 3)

final_0p3_offsets.drop('Unnamed: 0', axis=1, inplace=True)
final_0p5_offsets.drop('Unnamed: 0', axis= 1, inplace=True)

finalfillOffs = refineOffsetsAligned(final_0p3_offsets, final_0p5_offsets)
os.chdir(finalDir)
finalfillOffs.to_csv('final_offsets_0p3_fill_0p5.csv')


init_0p3_offsets = findOffsets(initial_aligned_positions, fname, 'initial_', '0p3', 2)
init_0p5_offsets = findOffsets(initial_aligned_positions, fname, 'initial_', '0p5', 3)

init_0p3_offsets.drop('Unnamed: 0', axis=1, inplace=True)
init_0p5_offsets.drop('Unnamed: 0', axis= 1, inplace=True)

initfillOffs = refineOffsetsAligned(init_0p3_offsets, init_0p5_offsets)
os.chdir(initDir)
initfillOffs.to_csv('init_offsets_0p3_fill_0p5.csv')

final_fill_init = refineOffsetsAligned(finalfillOffs, initfillOffs)
os.chdir('..')
final_fill_init.to_csv('final_offsets_init_fill.csv')

#final_0p5_offsets = findOffsets(final_aligned_positions, tfname, '', '0p5', 3)

#initial_0p3_offsets = findOffsets(initial_aligned_positions, 'initial_', '0p3', 2)

#initial_0p5_offsets = findOffsets(initial_aligned_positions, 'initial_', '0p5', 3)

