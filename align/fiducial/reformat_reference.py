#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 15:21:23 2019

@author: jonathanwhite
"""

import pandas as pd
import os


def reformat_reference(fname, position, preform_dir, num_ch):
    '''
    Reformats points from matlab code to that recognized by python code. Should be run after reformat_reference
    which initializes the directories.
    
    Inputs:
        fname - filename of points with string formating symbol, "%d", in place of position number
        position - the integer postion of the file of points to reformat
        preform_dir - the directory of points to reformat
    
    Outputs:
        Creates directory for 
    '''
    
    #fname = fname % position
    
    os.chdir(preform_dir)
    pnts = pd.read_csv(fname, index_col='ch')
    os.chdir('..')
    os.chdir('positions')
    
    # set up position directory
    pos_dir_name = 'pos%d' % position
    if pos_dir_name not in os.listdir():
        os.mkdir(pos_dir_name)
    os.chdir(pos_dir_name)


    for ch in range(1,num_ch+1):
        # set up directories
        chdir_name = 'ch%d' %ch
        if chdir_name not in os.listdir():
            os.mkdir(chdir_name)
        os.chdir(chdir_name)
        
        ch_points_dir_name =  'ch%d_points' % ch
        if ch_points_dir_name not in os.listdir():
            os.mkdir(ch_points_dir_name)
        os.chdir(ch_points_dir_name)
        
        # format points
        chPnts = pnts.loc[ch]
    
        df = pd.DataFrame()
        
        df['row'] = 2048 - chPnts['y']
        df['col'] = chPnts['x']
        df['z'] = chPnts['z']
        if 'int' in chPnts:
            df['amp'] = chPnts['int']
    
        df.to_csv('ch_%d_'%ch +fname, index = False)
        
        # return to position directory
        os.chdir('..')
        os.chdir('..')
        
    # return to home directory
    os.chdir('..')
    os.chdir('..')