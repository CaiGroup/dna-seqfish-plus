# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 08:33:29 2019

@author: jonat
"""

from dotAligner import dotAligner
from matplotlib import pyplot as plt
import numpy as np
from fitDots import fitPSF
import os

#datadir = 'E:\\Jonathan\\40rounds_test\\Data'
datadir = 'C:\\Users\\jonat\\Documents\\40rounds_test'#\\Data'

dal = dotAligner()
#dal = dotAligner()

#os.chdir('..')
os.chdir('C:\\Users\\jonat\\Documents\\40rounds_test')
dal.loadAlignedPSFs('20191021_alignedPSFs_v3_thr300.csv')#'Dots_400x400_bnds.pkl')


#dist = 2
os.chdir('Data\\HybCycle_0')
dal.compareHybAlignmentToRef(1, True)


'''
#dal.find_matching_dots(dist=dist)

#print(dal.matchedDots[0])

#find number of matches for each dot in the first hyb
lens = [len(list(filter(lambda x: x != None, md))) for md in dal.matchedDots]

# find the mean and sd of the amplitudes of matched PSFs
means = []
sds = []
amps = []
for md in dal.matchedDots:
    damps = []
    for dot in md:
        if type(dot) == list:
            print(dot)
        elif not dot == None:
            damps.append(dot[1].amp)
    sds.append(np.std(damps))
    means.append(np.mean(damps))
    amps.append(damps)

print('mean matches:', np.mean(lens), 'stdv:', np.std(lens))

plt.semilogy(np.array(lens), means, '.')
plt.xlabel('Number of Hybs (out of 40) in which a match to the dot was found')
plt.ylabel('Mean Amplitude of dot matches across hybs')
plt.title('Number of across dots in all hybs matched within 1 pixel to each dot in the first hyb at threshold 400')
plt.show()

dal.plotMatchedUnmatchedDots(directory='Data', nhybs=33, bnds=400)
'''

'''
fig, ax = plt.subplots(1,3)


ax[0].plot(lens)
ax[0].set_title(dist)

fitPSF.width = 10

canvas1 = np.zeros((2048, 2048))
print('len PSFs[0]:', len(dal.PSFs[0]))
for psf in dal.PSFs[0]: 
    canvas1  = psf.paint(canvas1, alt_amp = np.log(psf.amp))
    pass
ax[1].imshow(canvas1)
ax[1].set_title("Hyb 0 PSFs")



canvas2 = np.zeros((2048, 2048))
#for i, psf in enumerate(dal.sortedPSFs[0]):
for i, matched in enumerate(dal.matchedDots):
    psf = matched[0][1]
    # if (psf.row, psf.col) != dal.matchedDots[i][0][1]:
    #    print((psf.row, psf.col), 'unmatched')
    psf.paint(canvas2, alt_amp = lens[i])
ax[2].imshow(canvas2)
ax[2].set_title('Hyb 0 PSFs: Brightness proportional to nMatches')

plt.show()
'''