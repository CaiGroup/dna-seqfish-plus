% example script to get the threshold for J:\2019-10-21-E14-L22Rik-24hr-DNAFISHIF

% add paths
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
addpath('C:\github\streamline-seqFISH\src\FindThreshold\bfmatlab\', '-end');
addpath('C:\github\streamline-seqFISH\src\process_with_beads', '-end');

folderArray = 0:17;
experimentDir = 'E:\Yodai\DNAFISH+\2019-07-06-E14-DNA-seqFISH+rep2-2-RNAFISH-IF - Swapped';
experimentName = '2019-07-06-E14-DNA-seqFISH+rep2-2-RNAFISH';
position = 0;
typedots = 'log';
numCh = 2; % number of channels to process and threshold
processedimages = I; % If images are in the workspace
saveImages = false; % option to save images by channel
backgroundFolderName = [];


threshold = thresholdbych(experimentDir, experimentName, position, ...
folderArray, typedots, saveImages, backgroundFolderName, numCh, processedimages);
