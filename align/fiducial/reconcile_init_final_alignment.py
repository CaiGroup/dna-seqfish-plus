# -*- coding: utf-8 -*-
"""
Created on Thu Nov 14 18:19:41 2019

@author: jonat
"""

import os
import pandas as pd

def aggregateOffsets(directory, baseFname, ref, XYerror, Zerror):
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