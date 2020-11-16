# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:36:54 2019

@author: jonat
"""

import os
import sys
from reformat_reference import reformat_reference
from reformat_read_out import reformat_read_out
from align_20191121 import align
from multiprocessing import Pool
#sys.path.append('C:\\Users\\jonat\\PycharmProjects\\seqfishimageprocessing')
sys.path.append('C:\\gitlab\\seqfishimageprocessing')
#sys.path.append('C:\\Users\\jonat\\PycharmProjects\\seqfishimageprocessing')
from refAligner import RefAligner
import pandas as pd

#preformatted_dir = 'pre_formated_aligned2hyb0'
#ref_fname = 'ref-points-ForBeadAlignment-pos%d-radial3d-finalBeads-raw-intensity.csv'
#ro_fname = 'hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d.csv'
#savefname = '20191124_test_par.csv'
preformatted_dir = 'pre_formated'
home_dir = 'I:\\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\\points'
ref_fname = 'ref-points-ForBeadAlignment-pos%d-radial3d-initial_fiducial_markers-raw-intensity.csv'
ro_fname = 'hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped-UpdatedThresholdPointRatio-logdots-Pos0vs3.csv'
savefname = '20191125_pos%d_offsets_initial_roithresholdfiltered_UpdatedThresholdPointRatio-logdots-Pos0vs3.csv'


# search error should be longer than traversal then match error
n_longest_edges = 500
xy_search_error = 3 #strict: 2
z_search_error = 4 #strict: 3
xy_traversal_error = 3 
z_traversal_error = 5
xy_match_error = 0.7 #strict: 0.5
z_match_error = 1
min_edge_match = 5
min_dot_matches = 20
min_bright_prop = 0.5
max_bright_prop = 1.5
n_unmatch_give_up = 30

reformat = False

dals = []
alignments = []

position = int(sys.argv[1])
maxPosition = 2
channels = (1,2) #(1,2,3)
hybs = range(1, 81)
processes = 12
savefname = savefname % position


os.chdir(home_dir)
reformat_reference(ref_fname, position, preformatted_dir)

# format read out points
os.chdir(home_dir)
reformat_read_out(ro_fname, position, preformatted_dir)

os.chdir('positions')

ref_fname = 'ch_%d_' + ref_fname
ro_fname = 'ch_%d_' + ro_fname

for pos in (position,): #range(maxPosition+1):
    os.chdir(home_dir)
    os.chdir('positions')
    os.chdir('pos%d' % pos)
    chdals = []
    for ch in channels:
        os.chdir('ch%d' % ch)
        os.chdir('ch%d_points' % ch)
        print(os.listdir())
        ref_dots = ref_fname % (ch, pos)
        ro_dots = ro_fname % (ch, pos)
        dal = RefAligner(ro_dots, ref_dots)
        dal.set_n_longest_edges(n_longest_edges)
        dal.set_xy_search_error(xy_search_error)
        dal.set_z_search_error(z_search_error)
        dal.set_xy_traversal_error(xy_traversal_error)
        dal.set_z_traversal_error(z_traversal_error)
        dal.set_xy_match_error(xy_match_error)
        dal.set_z_match_error(z_match_error)
        dal.set_min_edge_match(min_edge_match)
        dal.set_min_dot_matches(min_dot_matches)
        dal.set_min_bright_prop(min_bright_prop)
        dal.set_max_bright_prop(max_bright_prop)
        dal.set_n_unmatch_give_up(n_unmatch_give_up)
        chdals.append(dal)
        os.chdir('..')
        os.chdir('..')
        for hyb in hybs:
            alignments.append((pos, ch, hyb))
dals.append(chdals)
os.chdir('..')
#os.chdir('..')
        
    



def alignPosChHyb(alignment):
    pos, ch, hyb = alignment
    
    print('pos', pos, ', ch', ch, ', hyb', hyb)
    offsets = dals[0][ch-1].align_hyb(hyb)

        
    to_return = [pos, ch, hyb] + offsets
    return to_return

    alignPosChHyb(alignment)
#alignPosChHyb(alignments[0])
#for alignment in alignments:
    '''
offsets = [alignPosChHyb(alignment) for alignment in alignments]

save_df = pd.DataFrame(data=offsets, columns=('pos','ch', 'hyb', 'row', 'column', 'z', 'row_se', 'column_se', 'z_se', 'n_matches'))
save_df.set_index(['pos','ch', 'hyb'])
save_df.to_csv(savefname)

'''    
if __name__ == '__main__':
    pool = Pool(processes = processes)
    print('Starting Alignments')
    #alignment = pool.map(alignPos, range(5))
    alignment = pool.map(alignPosChHyb, alignments)
    pool.close()
    pool.join()
    save_df = pd.DataFrame(data=alignment, columns=('pos','ch', 'hyb', 'row', 'column', 'z', 'row_se', 'column_se', 'z_se', 'n_matches'))
    save_df.set_index(['pos','ch', 'hyb'])
    save_df.to_csv(savefname)
    

    
