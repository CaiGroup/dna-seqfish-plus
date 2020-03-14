PathName = 'E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\rehyb_visualization_20200210';
experimentName = 'E14-rep2-rehyb_visualization_20200210';
folderArray = 0:3;
useBackground = false;
subtractBackground = false;
backgroundFolder = []; % default will be used
dapiRefPath = [PathName '\HybCycle_0']; % use HybCycle_0 in experiment directory
imageJBackSubtract = false;
saveProcessIms = true;
divideIms = false;
dim = '3d';

experimentDir = PathName; 
listing2 = dir([experimentDir '\HybCycle_0\*Pos*']);

for position = 0:length(listing2)-1

    % This should correct the rotation.
    I = preprocessimages(experimentName, experimentDir, position, folderArray, ...
        useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
        subtractBackground, saveProcessIms,divideIms, dim);
end

