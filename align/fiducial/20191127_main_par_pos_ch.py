# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 14:36:54 2019

@author: jonat
"""

import os
import sys
from reformat_reference import reformat_reference
from reformat_read_out import reformat_read_out
#from align_20191121 import align
from multiprocessing import Pool
#sys.path.append('C:\\Users\\jonat\\PycharmProjects\\seqfishimageprocessing')
sys.path.append('C:\\gitlab\\seqfishimageprocessing')
from refAligner import RefAligner
import pandas as pd
from datetime import datetime as dt


preformatted_dir = 'pre_formated'
home_dir = 'I:\\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\\points'
ref_fname = 'ref-points-ForBeadAlignment-pos%d-radial3d-initial_fiducial_markers-raw-intensity.csv'
ro_fname = 'hyb-points-0_80-ForBeadAlignment-pos%d-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH-UpdatedThreshold-12-10-2019-MedianIntRatio-logdots-Pos3-UseThresholdPos0.csv'
savefname = '20191125_pos%d_offsets_initial_roithresholdfiltered_UpdatedThresholdPointRatio-logdots-Pos0vs3_UpdatedThresholdPos0-12-11-2019'
#savefname = '20191202_'
datestr = dt.strftime(dt.now(), '%Y%m%d')
#paramstr = 'xyse%s_zse%s_xyte%s_zte%s_xyme%s_zme%s.csv'

n_longest_edges = 200
xy_search_error = 3
z_search_error = 3
xy_traversal_error = 2
z_traversal_error = 2
xy_match_error = 1
z_match_error = 1
min_edge_match = 5
min_dot_matches = 10
min_bright_prop = 0.2
max_bright_prop = 2
n_unmatch_give_up = 20

xyse_str = str(xy_search_error).replace('.', 'p')
zse_str = str(z_search_error).replace('.', 'p')
xyte_str = str(xy_traversal_error).replace('.', 'p')
zte_str = str(z_traversal_error).replace('.', 'p')
xyme_str = str(xy_match_error).replace('.', 'p')
zme_str = str(z_match_error).replace('.', 'p')
params = (xyse_str, zse_str, xyte_str, zte_str, xyme_str, zme_str)



reformat = False

dals = []
alignments = []

position = int(sys.argv[1])
savefname = savefname % position
maxPosition = 1
channels = (1,2)  # usually (1,2,3) channels
hybs = range(1, 81)
n_processes = 3

param_name = '_xyse%s_zse%s_xyte%s_zte%s_xyme%s_zme%s.csv' % params
savefname = savefname + datestr + param_name


os.chdir(home_dir)
reformat_reference(ref_fname, position, preformatted_dir)

# format read out points
os.chdir(home_dir)
reformat_read_out(ro_fname, position, preformatted_dir)

os.chdir('positions')

ref_fname = 'ch_%d_' + ref_fname
ro_fname = 'ch_%d_' + ro_fname

alignments = []
for pos in (position,):
    os.chdir('pos%d' % pos)
    chdals = []
    for ch in channels:
        os.chdir('ch%d' % ch)
        os.chdir('ch%d_points' % ch)
        print(os.listdir())
        ref_dots = ref_fname % (ch, pos)
        ro_dots = ro_fname % (ch, pos)
        dal = RefAligner(ro_dots, ref_dots)
        dal.set_n_longest_edges(n_longest_edges, all_pairs=False)
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
        alignments.append((pos, ch))
    dals.append(chdals)
    os.chdir('..')
#os.chdir('..')
        
    



def alignPosCh(alignment):
    pos, ch= alignment
    
    print('pos', pos, ', ch', ch)
    #offsets = dals[pos][ch-1].align(saveMatches=True)
    #ToDo: Change back to above when debugged
    offsets = dals[0][ch-1].align(saveMatches=True)
    offsets.insert(0, 'ch', ch)
    offsets.insert(0, 'pos', pos)
    dals[0][ch-1].save_matches(datestr + '_matches_pos%d_ch%d' % (pos, ch) + param_name)
     
    return offsets

#offsets = [alignPosCh(alignment) for alignment in alignments]
#offsets = []
#for alignment in alignments:

    #offsets.append(alignPosChHyb(alignment))


if __name__ == '__main__':
    pool = Pool(processes = n_processes)
    print('Starting Alignments')
    alignment = pool.map(alignPosCh, alignments)
    pool.close()
    pool.join()
    #save_df = pd.DataFrame(data=alignment, columns=('pos','ch', 'hyb', 'row', 'column', 'z', 'row_se', 'column_se', 'z_se', 'n_matches'))
    save_df = pd.concat(alignment)
    #save_df.set_index(['pos','ch', 'hyb'])
    save_df.to_csv(savefname)
    '''
    for ipos, pos in enumerate(dals):
        for ich, ch in enumerate(pos):
            ch.save_matches('matches_pos%d_ch%d.csv' % (ipos, ich))

    '''
    

    
