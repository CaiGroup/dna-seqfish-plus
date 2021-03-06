% repeat of the alignment to do intron analysis with 3d localization
% information.
% chromatic aberration is not corrected with the preprocessed images. Need
% to correct later with points.
% 
% changed alignment strategy as follows.
% 1) align all RNA hybs dapi to RNA hyb1 dapi
% 2) compute tform between RNA hyb1 dapi and DNA hyb1 dapi, and apply to
% all hybs
% 3) compute rotation tform between aligned RNA hyb1 dapi and DNA hyb1 dapi
% and apply to all hybs.
% This will minimize alignment errors across RNA hyb by hyb. need to check
% hyb1 RNA and DNA are aligned well.

experimentDir = 'F:\Yodai\DNA+\2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH-IF - Swapped';
experimentName = '2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH-IF';
folderArray = 0:20; % up to the last hyb that contains Neat1.
useBackground = false;
subtractBackground = false;
backgroundFolder = []; % default will be used
dapiRefPath = 'F:\Yodai\DNA+\2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH-IF - Swapped\HybCycle_0';
dapiRefPathDNA = 'F:\Yodai\DNA+\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\HybCycle_0';
imageJBackSubtract = true;
saveProcessIms = false;
divideIms = false;
dim = '3d';
listing = dir([dapiRefPath '\*MMStack_Pos*']);

for position = 0:length(listing)-1
    I = preprocessimagesRot_v2(experimentName, experimentDir, position, folderArray, ...
        useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
        subtractBackground, saveProcessIms,divideIms, dim, dapiRefPathDNA);
end
    
