addpath('C:\github\streamline-seqFISH\src\preprocessing\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');

experimentDir = 'E:\2019-09-14-brain-rep3-2-DNA-FISH';
experimentName = '2019-09-14-brain-rep3-2-DNA-FISH';

folderArray = 0:89;
useBackground = true;
subtractBackground = false;
backgroundFolder = 'initial_background'; % default will be used
dapiRefPath = ''; % use HybCycle_0 in experiment directory
imageJBackSubtract = true;
saveProcessIms = false;
divideIms = false;
dim = '3d';

for position = 0
    preprocessimages(experimentName, experimentDir, position, folderArray, ...
        useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
        subtractBackground, saveProcessIms,divideIms, dim);
end