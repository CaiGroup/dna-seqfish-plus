experimentDir = 'I:\2019-09-09-brain-rep2-2-DNAFISH';
experimentName = '2019-09-09-brain-rep2-2-DNAFISH';

folderArray = 0:92;
useBackground = true;
subtractBackground = false;
backgroundFolder = 'initial_background'; % default will be used
dapiRefPath = ''; % use HybCycle_0 in experiment directory
imageJBackSubtract = true;
saveProcessIms = false;
divideIms = false;
dim = '2d';

for position = 3
    [I] = preprocessimages(experimentName, experimentDir, position, folderArray, ...
        useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
        subtractBackground, saveProcessIms,divideIms, dim);
end