
addpath('C:\github\streamline-seqFISH\src\preprocessing\', '-end');
addpath('C:\github\streamline-seqFISH\src\preprocessing\bfmatlab\', '-end');

experimentDir = 'I:\2019-09-04-brain-rep2-2-RNAFISH';
experimentName = '2019-09-04-brain-rep2-2-RNAFISH';
folderArray = 0:36; % 38 and 39 have 49 z-slices
useBackground = true;
subtractBackground = false;
imageJBackSubtract = true;
saveProcessIms = false;
divideIms = false;
dim = '3d'; % use 2d for rounding z transformation if lessthan 0.5; use 3d if not...
finalFixedPath = 'I:\2019-09-09-brain-rep2-2-DNAFISH\HybCycle_0';
finalMovingPath = 'I:\2019-09-04-brain-rep2-2-RNAFISH\final_alignment';
%backgroundFolder = 'initial_background';
alignCh = 3;

backgroundFolder = 'final_background';
refPath = 'I:\2019-09-04-brain-rep2-2-RNAFISH\final_alignment';
dapiRefPath = 'I:\2019-09-04-brain-rep2-2-RNAFISH\final_alignment';

for position = 0%1:4
    I = preprocessimmunoimages(experimentName, experimentDir, position, folderArray, ...
    useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
    subtractBackground, saveProcessIms, divideIms, dim, finalFixedPath, finalMovingPath, alignCh);
end