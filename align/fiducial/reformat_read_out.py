#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 15:13:18 2019

@author: jonathanwhite
"""

import pandas as pd
import os

def reformat_read_out(fname, position, preform_dir, num_ch):
    '''
    reformats points from matlab code to that recognized by python code. Should be run after reformat_reference
    which initializes the directories.
    
    Inputs:
        fname - filename of points with string formating symbol, "%d", in place of position number
        position - the integer postion of the file of points to reformat
        preform_dir - the directory of points to reformat
    
    Outputs:
        Creates directory for 
    '''
    os.chdir(preform_dir)
    
    # insert position number into file name
    #fname = fname % position

    pnts = pd.read_csv(fname, index_col = 'ch')
    
    #print(pnts.iloc[1:5])
    
    # go back to home directory
    os.chdir('..')
    
    os.chdir('positions')
    os.chdir('pos%d' % position)
    
    for ch in range(1,num_ch+1): #(1,2,3):
        os.chdir('ch%d' % ch)
        os.chdir('ch%d_points' % ch)
        df = pd.DataFrame()
        chPnts = pnts.loc[ch]
    
        #make dot index
    
        #df = pd.DataFrame()
    
        df['Hyb'] = chPnts['hyb']
        #df['ch'] = pnts
        df['row'] = 2048 - chPnts['y']
        df['col'] = chPnts['x']
        df['z'] = chPnts['z']
        df['amp'] = chPnts['int']
        
        ch_fname = 'ch_%d_'%ch + fname
    
        df.to_csv(ch_fname, index=False)
        
        # return to position directory
        os.chdir('..')
        os.chdir('..')
        
    # return to home directory
    os.chdir('..')
    os.chdir('..')
