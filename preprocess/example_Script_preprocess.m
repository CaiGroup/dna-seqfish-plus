% preprocess images using DAPI channel

experimentDir = 'I:\2019-09-04-rep2-2-RNAFISH';
experimentName = '2019-09-04-rep2-2-RNAFISH';
folderArray = 0:36; % 38 and 39 have 49 z-slices
useBackground = true;
subtractBackground = false;
imageJBackSubtract = true;
saveProcessIms = false;
divideIms = false;
dim = '3d'; % dimension to align images
finalFixedPath = 'I:\2019-09-09-rep2-2-DNAFISH\HybCycle_0';
finalMovingPath = 'I:\2019-09-04-rep2-2-RNAFISH\final_alignment';
alignCh = 3;

backgroundFolder = 'final_background';
refPath = 'I:\2019-09-04-rep2-2-RNAFISH\final_alignment';
dapiRefPath = 'I:\2019-09-04-rep2-2-RNAFISH\final_alignment';

for position = 0:4
    I = preprocessimmunoimages(experimentName, experimentDir, position, folderArray, ...
    useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
    subtractBackground, saveProcessIms, divideIms, dim, finalFixedPath, finalMovingPath, alignCh);
end