% repeat of the alignment to do intron analysis with 3d localization
% information.
% chromatic aberration is not corrected with the preprocessed images. Need
% to correct later with points.

experimentDir = 'F:\Yodai\DNA+\2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH-IF - Swapped';
experimentName = '2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH-IF';
folderArray = 0:17; %0:17 up to the last hyb that contains intron.
useBackground = false;
subtractBackground = false;
backgroundFolder = []; % default will be used
dapiRefPath = 'F:\Yodai\DNA+\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\HybCycle_0';
imageJBackSubtract = true;
saveProcessIms = false;
divideIms = false;
dim = '3d';
listing = dir([dapiRefPath '\*MMStack_Pos*']);

for position = 0:length(listing)-1
    I = preprocessimagesRot(experimentName, experimentDir, position, folderArray, ...
        useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
        subtractBackground, saveProcessIms,divideIms, dim);
end
    
