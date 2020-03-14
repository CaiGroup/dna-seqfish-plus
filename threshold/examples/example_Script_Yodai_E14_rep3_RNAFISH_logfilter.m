% example script to get the threshold for J:\2019-10-21-E14-L22Rik-24hr-DNAFISHIF

% add paths
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
addpath('C:\github\streamline-seqFISH\src\FindThreshold\bfmatlab\', '-end');
addpath('C:\github\streamline-seqFISH\src\process_with_beads', '-end');

folderArray = 0:20;
experimentDir = 'F:\Yodai\DNA+\2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH-IF - Swapped';
experimentName = '2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH';
position = 0;
typedots = 'log';
numCh = 2; % number of channels to process and threshold
processedimages = I; % If images are in the workspace
saveImages = false; % option to save images by channel
backgroundFolderName = [];


threshold = thresholdbych(experimentDir, experimentName, position, ...
folderArray, typedots, saveImages, backgroundFolderName, numCh, processedimages);
