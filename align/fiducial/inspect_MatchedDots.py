from dotAligner import dotAligner
from matplotlib import pyplot as plt
import numpy as np
from fitDots import fitPSF
import os
import seaborn as sns

#datadir = 'E:\\Jonathan\\40rounds_test\\Data'

#dal = dotAligner(directory = datadir, nhybs = 3)
dal = dotAligner()

#os.chdir('..')
os.chdir('C:\\Users\\jonat\\Documents\\40rounds_test')
#dal.loadPSFs('Dots_400x400_bnds.pkl')
dal.loadPSFs('Dots_2048x2048.pkl')

dal.graphMatchRecurringDots()

dal.findOffsets()

dist = 2
dal.find_matching_dots(dist=dist)

ampzs, summaryStats, xs, ys, amps = dal.evaluateWellMatchedDots(directory='Data')

plt.hist(ampzs, bins=80)
plt.title('Brightness distribution of dots between hybs (Pooled)')
plt.xlabel('Dot Brightness Z score')
plt.ylabel('Counts')
plt.show()

ss2 = summaryStats[summaryStats['Mean_Brightness']>3000]
ss = summaryStats


#plt.hist(ss['XSTDV']*110)
plt.hist(ss2['XSTDV']*110)
plt.xlabel('Standard Deviation of X localization (nm)')
plt.ylabel('Counts')
plt.title('Variation in X localizations for dots brighter than 3000 found in 40 hybs')
#plt.title('Variation in X localizations for dots found in 40 hybs')
plt.show()

#plt.hist(ss['YSTDV']*110)
plt.hist(ss2['YSTDV']*110)
plt.xlabel('Standard Deviation of Y localization (nm)')
plt.ylabel('Counts')
plt.title('Variation in Y localizations for dots brighter than 3000 found in 40 hybs')
#plt.title('Variation in Y localizations for dots found in 40 hybs')
plt.show()

'''
fig, ax = plt.subplots(4)
for i in range(4):
    start = i*144
    end = start+144
    ampsportion =amps[:, start:end]
    sns.violinplot(x='Matched PSFs in each Hyb', y='Brightness', data=ampsportion.transpose(), ax=ax[i])
'''
sns.violinplot(data=amps)
plt.xlabel("Hyb Number")
plt.ylabel('Dot Brightness Distribution')
plt.title('Distibution of dot brightness in each hyb')
#plt.yscale('log')
plt.show()

'''
#plt.errorbar(ss['RMSE_Mean'], ss['Mean_Brightness'], xerr=ss['RMSE_STDV'], yerr=ss['ampSTDV'], fmt='.')
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
plt.show()'''

#log 10 errorbar
mb = ss['Mean_Brightness']
errorUpper = ss['Mean_Brightness'] + ss['ampSTDV']
errorLower = ss['Mean_Brightness'] - ss['ampSTDV']
MBdRMSE = ss['Mean_Brightness_div_Mean_RMSE']
ss['localization_radius_sd'] = np.sqrt(ss['XSTDV']**2 + ss['YSTDV']**2)
#plt.errorbar(ss['Mean_Brightness_div_Mean_RMSE'], logmb, yerr=errorBar, fmt='.')
for dot in range(len(mb)):
    #x = [ss['Mean_Brightness_div_Mean_RMSE'][dot]]*2
    x = [MBdRMSE[dot]] * 2
    y = [errorLower[dot], errorUpper[dot]]
    plt.plot(x,y,'-k')
#plt.plot(MBdRMSE, mb,'.b')
sns.scatterplot(x='Mean_Brightness_div_Mean_RMSE', y='Mean_Brightness', hue='localization_radius_sd', data=ss, sizes=(100,100))
plt.xlabel('Mean Brightness divided by mean RMSE (A.U.)')
plt.ylabel('Mean Brightness (A.U.)')
plt.xscale("Log")
plt.yscale('Log')
plt.legend()
plt.title('Brightness and RMSE in dots matched across 40 hybs. Error bars 1 SD.')
plt.show()

#xdisps = xs[:, 1:] - xs[:, 0]
xdisps = np.copy(xs)
xdisps = xdisps[summaryStats['Mean_Brightness']>3000,:]
x0 = np.copy(xdisps[:,0])
for i in range(np.shape(xdisps)[1]):
    xdisps[:,i] -= x0 #xs[:,0]

#ydisps = ys[:, 1:] - ys[:, 0]
ydisps = np.copy(ys)
ydisps = ydisps[summaryStats['Mean_Brightness']>3000,:]
y0 = np.copy(ydisps[:,0])
for i in range(np.shape(ydisps)[1]):
    ydisps[:,i] -= y0 #ys[:,0]
    
rdisps = rdisps = np.sqrt(xdisps[:,1:]**2 + ydisps[:,1:]**2)

plt.hist(rdisps.flatten()*110, bins=80)
plt.xlabel('Distance between initial dot localization and subsequent localizations (nm)')
plt.ylabel('Counts')
plt.title('Localization Accuracy of dots brighter than 3000 counts found in all 40 hybs')
#plt.title('Localization Accuracy of dots found in all 40 hybs')
plt.show()