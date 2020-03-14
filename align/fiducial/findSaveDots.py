# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 10:07:02 2019

@author: Long Cai - 2
"""

from dotAligner import dotAligner

z = 1
channel = 1

fname ='NIH3T3_thr900_Pos1_z%d_ch%d_Dots_2048x2048.pkl' % (z, channel)
print("finding:", fname)

#dal = dotAligner('Data', 500)
#datadir = 'E:\\Jonathan\\40rounds_test\\Data'
#datadir = 'E:\\Jonathan\\Linus_10k_cleared_080918_NIH3T3'
datadir = '/Users/jonathanwhite/Documents/NIH3T3'

dal9 = dotAligner(directory=datadir,  z=z, channel=channel)
dal9.findPSFs(bnd=2048, threshold=900)

dal9.savePSFs(fname)