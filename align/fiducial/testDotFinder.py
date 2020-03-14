# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 11:59:57 2019

@author: jonat
"""

from fitDots import dotFinder
import skimage
import os
from matplotlib import pyplot as plt
import numpy as np

#os.chdir('C:\\Users\\jonat\\Dropbox\\_Jonathan\\cailab\\Dot_alignment')
#im = skimage.io.imread('test_im.tif')
os.chdir('C:\\Users\\jonat\\Documents\\40rounds_test\\Data\\HybCycle_0')
im = skimage.io.imread('MMStack_Pos0.ome.tif')
#im = skimage.io.imread('test_im.tif')
bnd = 200
im = im[0][:bnd,:bnd, 1]
imorig = np.copy(im)

os.chdir('..')

os.chdir('final_background')
bg = skimage.io.imread('MMStack_Pos0.ome.tif')
bg = bg[0][:bnd,:bnd, 1]
bgorig = np.copy(bg)


dfr = dotFinder(im, bg, 2, meanBgPropThresh = 3/4)
dfr.fitSinglePSFs()

fig, ax = plt.subplots(1,4)
ax[0].imshow(imorig)
ax[0].set_title('image')
ax[1].imshow(dfr.imbgsub)
ax[1].set_title('background subtracted image')
ax[2].imshow(dfr.im_filt)
ax[2].set_title('Background subtraction and Morphological opening')
ax[3].imshow(bgorig)
ax[3].set_title('background')
plt.show()


# background comparison plot
fig, ax = plt.subplots(1,2)
ax[0].imshow(im)
ax[1].imshow(bg)
plt.show()

#im = im[2,:60,:60]
#im = im[2,60:120,60:120]
#im = im[2,80:130,65:115]
#im = im[2,95:120,90:115]
#im = im[2,:200,:200]
#im = im[2,0:100,0:100]
#im = im[2,0:300,0:300]
#threshold = 00


dfr.plotSinglePSFs(True, logScale = False)
