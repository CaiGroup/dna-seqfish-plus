#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 15:13:18 2019

@author: jonathanwhite
"""

import pandas as pd

pnts = pd.read_csv('rawPointsForBeadAlignment-ch1.csv')

print(pnts.iloc[1:5])

dots = []

hybLast = -1
dot= 0

for i, row in pnts.iterrows():
    if row['hyb'] != hybLast:
        dot = 0
    else:
        dot += 1
    dots.append(dot)

df = pd.DataFrame()

df['Hyb'] = pnts['hyb']
df['Dot'] = dots
df['row'] = 2048 - pnts['y']
df['col'] = pnts['x']
df['z'] = pnts['z']
df['amp'] = pnts['int']

df.to_csv('points4Alignment_ch1.csv', index=False)

        