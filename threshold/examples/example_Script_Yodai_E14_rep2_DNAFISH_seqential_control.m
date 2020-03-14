% add paths

addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
experimentPath = 'E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH';
experimentLabel = '2020-01-16-sequential-control-hyb81-84-test';
experimentDir = fullfile(experimentPath,experimentLabel);
position = 0;
typedots = 'log';
numCh = 3; % number of channels to process and threshold
processedimages = I; % unaligned preprocessed images.
folderArray = 0:size(processedimages,1)-1;
saveImages = false; % option to save images by channel
backgroundFolderName = [];


threshold = thresholdbych(experimentDir, experimentName, position, ...
folderArray, typedots, saveImages, backgroundFolderName, numCh, processedimages);
