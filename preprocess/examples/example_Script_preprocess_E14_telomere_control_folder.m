%%%%
% This code didn't work well, somehow.
%%%%
PathName = 'G:\Yodai\DNA+\2019-11-14-Telomere-control\images\for supp figure\2-2 ligationfixation-700ms';
experimentName = '2019-11-14-Telomere-control';
folderArray = 0:1;
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
    I = preprocessimagesRot(experimentName, experimentDir, position, folderArray, ...
        useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
        subtractBackground, saveProcessIms,divideIms, dim);
end
    
