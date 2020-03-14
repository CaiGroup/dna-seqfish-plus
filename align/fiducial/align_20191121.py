# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:29:26 2019

@author: jonat
"""

import sys
sys.path.append('C:\\gitlab\\seqfishimageprocessing')
#sys.path.append('/Users/jonathanwhite/PycharmProjects/seqfishimageprocessing')

import os
from refAligner import RefAligner


os.chdir('..')

def align(n_longest_edges,
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
          hyb_fname, 
          position):
    
    saveOffsetsName = 'offsets_ch%d_xyse%s_zse%s_xyte%s_zte%s_xyme%s_zme%s_mem%d_mdm%d_nle%d_' + hyb_fname
    s_xy_search_error = str(xy_search_error).replace('.', 'p')
    s_xy_traversal_error = str(xy_traversal_error).replace('.', 'p')
    s_xy_match_error = str(xy_match_error).replace('.', 'p')
    s_z_search_error = str(z_search_error).replace('.', 'p')
    s_z_traversal_error = str(z_traversal_error).replace('.', 'p')
    s_z_match_error = str(z_match_error).replace('.', 'p')
    #sstopRefProp = str(stopRefProp).replace('.', 'p')
    ref_fname = 'ch_%d_' + ref_fname
    hyb_fname = 'ch_%d_' + hyb_fname
    
    for ch in (1,2,3):
        os.chdir('ch%d' % ch)
        os.chdir('ch%d_points' % ch)
        ref_dots = ref_fname % (ch, position)
        ro_dots = hyb_fname % (ch, position)
        
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
        
        
        #return to channel directory
        os.chdir('..')
        
        #find offsets
        dal.align()

        # set up directory to save points
        offsets_dir_name = 'ch%d_offsets' % ch
        if offsets_dir_name not in os.listdir():
            os.mkdir(offsets_dir_name)
        os.chdir(offsets_dir_name)
        
        save_name = saveOffsetsName % (ch, s_xy_search_error, s_z_search_error, s_xy_traversal_error, s_z_traversal_error,
            s_xy_match_error, s_z_match_error, min_edge_match, min_dot_matches, n_longest_edges, position)
        dal.save_offsets(save_name)
        #dal.saveMatchesDF(saveMatchesDFName % (ch, maxXYError, maxZError, minMatch, nLongestCheck))
        
        #dal.saveAlignedPSFs(saveAlignedName % (ch, maxXYError, maxZError, minMatch))
        
        #return to positions directory
        os.chdir('..')
        os.chdir('..')
        
        