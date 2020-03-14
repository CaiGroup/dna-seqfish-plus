from dotAligner import dotAligner
from matplotlib import pyplot as plt
import numpy as np
from fitDots import fitPSF
import os
import seaborn as sns
import pandas as pd

datadir = 'E:\\Jonathan\\40rounds_test\\Data'

'''
#dal = dotAligner(directory = datadir, nhybs = 3)
dal2 = dotAligner()

#os.chdir('..')

#dal.loadPSFs('Dots_400x400_bnds.pkl')
dal2.loadPSFs('Dots_2048x2048.pkl')

dal2.matchDots()

dal2.findOffsets()
'''
os.chdir('C:\\Users\\jonat\\Documents\\40rounds_test')

dal = dotAligner()
dal.loadPSFs('Dots_2048x2048.pkl')
dal.loadOffsets('dapi-tform-Linus40Test-2019-10-07.csv')
#diffxOffsets = np.array(dal.offsets['x']) - np.array(dal2.offsets['x'])
#diffyOffsets = np.array(dal.offsets['y']) - np.array(dal2.offsets['y'])
#print('Difference in X offsets', diffxOffsets)
#print('Difference in Y offsets', diffyOffsets)


dist = 2
dal.find_matching_dots(dist=dist)



ampzs, summaryStats, xs, ys, amps = dal.evaluateWellMatchedDots(directory='Data')

plt.hist(ampzs, bins = 80)

ss = summaryStats

logmb = np.log2(ss['Mean_Brightness'])
errorUpper = np.log2(ss['Mean_Brightness'] + ss['ampSTDV'])
errorLower = np.log2(ss['Mean_Brightness'] - ss['ampSTDV'])
logMBdRMSE = np.log2(ss['Mean_Brightness_div_Mean_RMSE'])
#plt.errorbar(ss['Mean_Brightness_div_Mean_RMSE'], logmb, yerr=errorBar, fmt='.')
for dot in range(len(logmb)):
    #x = [ss['Mean_Brightness_div_Mean_RMSE'][dot]]*2
    x = [logMBdRMSE[dot]] * 2
    y = [errorLower[dot], errorUpper[dot]]
    plt.plot(x,y,'-k')
plt.plot(logMBdRMSE, logmb,'.b')
plt.xlabel('Log2 Mean Brightness divided by mean RMSE')
plt.ylabel('Log2 Mean Brightness')
plt.title('Brightness and RMSE in dots matched across 40 hybs. Error bars 1 SD.')
plt.show()


#fig, ax = plt.subplots(4)
#sns.violinplot(x='Matched PSFs in each Hyb', y='Brightness', data=pd.DataFrame(amps))
#plt.show()


hybMeanBright = []
hybSTDVBright = []

for i in range(40):
    hybMeanBright.append(np.mean(amps[:,i]))
    hybSTDVBright.append(np.std(amps[:,i]))

plt.errorbar(np.arange(1, 41), hybMeanBright, yerr=hybSTDVBright)
plt.xlabel('Hyb Number')
plt.ylabel('Mean proghtess of dots')
plt.show()


ss2 = summaryStats[summaryStats['Mean_Brightness']>1000]

plt.hist(ss2['XSTDV']*110)
plt.xlabel('Standard Deviation of X localization (nm)')
plt.ylabel('Counts')
plt.title('Variation in X localizations for dots found in 40 hybs')
plt.show()

plt.hist(ss2['YSTDV']*110)
plt.xlabel('Standard Deviation of Y localization (nm)')
plt.ylabel('Counts')
plt.title('Variation in Y localizations for dots found in 40 hybs')
plt.show()

#xdisps = xs[:, 1:] - xs[:, 0]
xdisps = np.copy(xs)
xdisps = xdisps[summaryStats['Mean_Brightness']>1000,:]
x0 = np.copy(xdisps[:,0])
for i in range(np.shape(xdisps)[1]):
    xdisps[:,i] -= x0 #xs[:,0]

#ydisps = ys[:, 1:] - ys[:, 0]
ydisps = np.copy(ys)
ydisps = ydisps[summaryStats['Mean_Brightness']>1000,:]
y0 = np.copy(ydisps[:,0])
for i in range(np.shape(ydisps)[1]):
    ydisps[:,i] -= y0 #ys[:,0]
    
rdisps = np.sqrt(xdisps[:,1:]**2 + ydisps[:,1:]**2)

plt.hist(rdisps.flatten()*110, bins = 80)
plt.xlabel('Distance between initial dot localization and subsequent localizations (nm)')
plt.ylabel('Counts')
plt.title('Localization Accuracy of dots brighter than 1000 counts found in all 40 hybs')
plt.show()