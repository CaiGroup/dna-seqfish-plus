PathName = 'G:\Yodai\DNA+\2019-11-14-Telomere-control\images\for supp figure';
experimentName = '2019-11-14-Telomere-control';
folderArray = 0:1;
useBackground = false;
subtractBackground = false;
backgroundFolder = []; % default will be used
imageJBackSubtract = false;
saveProcessIms = true;
divideIms = false;
dim = '3d';

listing = dir([PathName '\*-*']);



for condition = 1:length(listing)
    experimentDir = [PathName '\' listing(condition).name];
    listing2 = dir([experimentDir '\HybCycle_0\*Pos*']);
    dapiRefPath = [experimentDir '\HybCycle_0'];
    for position = 0:length(listing2)-1
        
        % This should correct the rotation.
        I = preprocessimagesRot(experimentName, experimentDir, position, folderArray, ...
            useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
            subtractBackground, saveProcessIms,divideIms, dim);
    end
    
end