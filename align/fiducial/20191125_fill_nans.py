# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 09:38:47 2019

@author: jonat
"""

import pandas as pd
import os
from align_spec_hyb import align

directory = 'C:\\Users\\jonat\\Box\\Jonathan\\2019-09-09-brain-rep2-2-DNAFISH\\phase2_many_hybs\\positions'
file_to_fill = '20191124_test_par.csv'
filled_file = '20191124_test_par_filled.csv'

os.chdir(directory)
to_fill = pd.read_csv(file_to_fill)

to_fill.set_index(['pos','ch','hyb'], inplace=True)

not_aligned = to_fill.isnull().any(axis=1)

to_align = to_fill.loc[not_aligned]

n_longest_edges = 1000
xy_search_error = 3
z_search_error = 4
xy_traversal_error = 2
z_traversal_error = 3
xy_match_error = 0.7
z_match_error = 1
min_edge_match = 10
min_dot_matches = 50
min_bright_prop = 0.5
max_bright_prop = 1.5
n_unmatch_give_up = 50
ref_fname = 'ref-points-ForBeadAlignment-pos%d-radial3d-finalBeads-raw-intensity.csv'
ro_fname = 'hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d.csv'
offsets = []

for i in to_align.index:
    pos, ch, hyb = i
    os.chdir('pos%d' % pos)
    offset = align(
          n_longest_edges,
          xy_search_error, 
          z_search_error,
          xy_traversal_error,
          z_traversal_error,
          xy_match_error, 
          z_match_error, 
          min_edge_match, 
          min_dot_matches, 
          min_bright_prop,
          max_bright_prop,
          n_unmatch_give_up,
          ref_fname, 
          ro_fname, 
          pos,
          ch,
          hyb
            )
    to_fill.loc[i, 'row':] = offset
    offsets.append(offset)
    os.chdir('..')
    
to_fill.to_csv(filled_file)
    

