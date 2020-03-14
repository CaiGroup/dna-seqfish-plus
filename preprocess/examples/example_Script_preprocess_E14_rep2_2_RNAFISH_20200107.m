% repeat of the alignment to do intron analysis with 3d localization
% information.
% chromatic aberration is not corrected with the preprocessed images. Need
% to correct later with points.

experimentDir = 'E:\Yodai\DNAFISH+\2019-07-06-E14-DNA-seqFISH+rep2-2-RNAFISH-IF - Swapped';
experimentName = '2019-07-06-E14-DNA-seqFISH+rep2-2-RNAFISH-IF';
folderArray = 0:17; %0:17 up to the last hyb that contains intron.
useBackground = false;
subtractBackground = false;
backgroundFolder = []; % default will be used
dapiRefPath = 'E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\HybCycle_0';
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
    
