# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 16:51:50 2019
@author: jonat
"""

import os
import sys
from reformat_reference import reformat_reference
from reformat_read_out import reformat_read_out
#from align_20191121 import align
from multiprocessing import Pool
#sys.path.append('C:\\Users\\jonat\\PycharmProjects\\seqfishimageprocessing')
#sys.path.append('C:\\gitlab\\seqfishimageprocessing')
from refAligner import RefAligner
import pandas as pd
from datetime import datetime as dt

preformatted_dir = 'pre_formated'
home_dir = sys.argv[1] #'I:\\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\\points'
print('home dir: ' + home_dir)
ref_fname = sys.argv[2] #'ref-points-ForBeadAlignment-pos%d-radial3d-initial_fiducial_markers-raw-intensity.csv'
print('refname: ' + ref_fname)
ro_fname = sys.argv[3] #'hyb-points-0_80-ForBeadAlignment-pos%d-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH-UpdatedThreshold-12-10-2019-MedianIntRatio-logdots-Pos3-UseThresholdPos0.csv'
print('readout: ' + ro_fname)
savefname = sys.argv[4] #'20191125_pos%d_offsets_initial_roithresholdfiltered_UpdatedThresholdPointRatio-logdots-Pos0vs3_UpdatedThresholdPos0-12-11-2019'
print('savename: ' + savefname)
position = int(sys.argv[5])
print('position: ' + sys.argv[5])
datestr = dt.strftime(dt.now(), '%Y%m%d')
#paramstr = 'xyse%s_zse%s_xyte%s_zte%s_xyme%s_zme%s.csv'

n_longest_edges = int(sys.argv[6]) # E14 20; brain = 200


xy_search_error = 5
z_search_error = 2
xy_traversal_error = 3
z_traversal_error = 2
xy_match_error = 1
z_match_error = 1
min_edge_match = int(sys.argv[7])  # no more than 10%-20% NOT < 3; brain = 10; E14 = 3
min_dot_matches = int(sys.argv[8]) # E14 4; brain = 20 or 50 (half amount of fiducial markers)
min_bright_prop = 0.2
max_bright_prop = 2
n_unmatch_give_up = 30


xyse_str = str(xy_search_error).replace('.', 'p')
zse_str = str(z_search_error).replace('.', 'p')
xyte_str = str(xy_traversal_error).replace('.', 'p')
zte_str = str(z_traversal_error).replace('.', 'p')
xyme_str = str(xy_match_error).replace('.', 'p')
zme_str = str(z_match_error).replace('.', 'p')
params = (xyse_str, zse_str, xyte_str, zte_str, xyme_str, zme_str)

param_name = '_xyse%s_zse%s_xyte%s_zte%s_xyme%s_zme%s.csv' % params
savefname = savefname + datestr + param_name
save_match_name = datestr + 'matches_pos%d_ch%d' + param_name

reformat = True

dals = []
alignments = []

maxPosition = 0
numch = int(sys.argv[9])
channels = range(1,numch+1)
hybs = range(1, int(sys.argv[10]))
n_processes = 15

os.chdir(home_dir)

if reformat:
    print('Reformatting')
    # for each  position 
    for pos in (0,):
        # format reference points
        os.chdir(home_dir)
        reformat_reference(ref_fname % pos, pos, preformatted_dir, numch)
    
        # format read out points
        os.chdir(home_dir)
        reformat_read_out(ro_fname % pos, pos, preformatted_dir, numch)

os.chdir('positions')

ref_fname = 'ch_%d_' + ref_fname
ro_fname = 'ch_%d_' + ro_fname

alignments = []
for pos in (0,):#range(maxPosition + 1):
    os.chdir('pos%d' % pos)
    chdals = []
    for ch in channels:
        os.chdir('ch%d' % ch)
        os.chdir('ch%d_points' % ch)
        print(os.listdir())
        ref_dots = ref_fname % (ch, pos)
        ro_dots = ro_fname % (ch, pos)
        #ref_dots = ref_fname
        #ro_dots = ro_fname
        try:
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
        except FileNotFoundError as e:
            print(e)
            raise(e)
        os.chdir('..')
        os.chdir('..')
        alignments.append((pos, ch))
    dals.append(chdals)
    os.chdir('..')
#os.chdir('..')
        
    



def alignPosCh(alignment):
    pos, ch= alignment
    print('aligning')
    
    print('pos', pos, ', ch', ch)
    try:
        offsets = dals[pos][ch-1].align(saveMatches=True)
        offsets.insert(0,'ch', ch)
        offsets.insert(0,'pos', pos)
        fname = save_match_name % (position, ch)
        print('saving', fname)
        dals[pos][ch-1].save_matches(fname)
        return offsets
    except IndexError as e:
        print(e)
        
    

if __name__ == '__main__':
    pool = Pool(processes = n_processes)
    print('Starting Alignments')
    alignment = pool.map(alignPosCh, alignments)
    #alignment = map(alignPosCh, alignments)
    #alignment = [alignPosCh(a) for a in alignments[0:3]]
    pool.close()
    pool.join()
    #save_df = pd.DataFrame(data=alignment, columns=('pos','ch', 'hyb', 'row', 'column', 'z', 'row_se', 'column_se', 'z_se', 'n_matches'))
    save_df = pd.concat(alignment)
    #save_df.set_index(['pos','ch', 'hyb'])
    #save_df.to_csv(savefname, index=False)
    
    # check for nans
        
    save_df.set_index(['pos','ch','hyb'], inplace=True)
    
    not_aligned = save_df.isnull().any(axis=1)
    
    to_align = save_df.loc[not_aligned]
    
    n_longest_edges = n_longest_edges * 5
    xy_search_error = 6
    z_search_error = 3
    xy_traversal_error = 5
    z_traversal_error = 3
    xy_match_error = 1.5
    z_match_error = 2
    min_edge_match = min_edge_match
    min_dot_matches = min_dot_matches
    min_bright_prop = 0.5
    max_bright_prop = 1.5
    n_unmatch_give_up = 50

    offsets = []
    
    
    for i in to_align.index:
        pos, ch, hyb = i
        os.chdir('pos%d' % pos)
        dal = dals[pos][ch-1]
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
        offset = dal.align_hyb(hyb)
        save_df.loc[i, 'row':] = offset
        offsets.append(offset)
        os.chdir('..')
        
    
        
    not_aligned = save_df.isnull().any(axis=1)
    to_get_channel_ave = save_df.loc[not_aligned]
    print('filling in with other channels...')
    for i in to_get_channel_ave.index:
        print(i)
        pos, ch, hyb = i
        other_channel = []
        for j in (1, 2):
            if j != ch:
                other_channels = (pos, j, hyb)
        print(save_df)
        save_df.loc[i,'row':'z'] = list(save_df.loc[other_channels,'row':'z'])
        
        
    print("Saving", savefname)
    save_df.to_csv(savefname)