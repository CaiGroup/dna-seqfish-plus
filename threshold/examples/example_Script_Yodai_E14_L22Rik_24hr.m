% example script to get the threshold for J:\2019-10-21-E14-L22Rik-24hr-DNAFISHIF

% add paths
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
addpath('C:\github\streamline-seqFISH\src\FindThreshold\bfmatlab\', '-end');
addpath('C:\github\streamline-seqFISH\src\process_with_beads', '-end');

folderArray = 0:80;
experimentDir = 'E:\2019-10-21-E14-L22Rik-24hr-DNAFISHIF';
experimentName = '2019-10-21-E14-L22Rik-24hr-DNAFISHIF';
position = 0;
typedots = 'log';
numCh = 3; % number of channels to process and threshold
I = []; % If images are in the workspace
saveImages = false; % option to save images by channel
backgroundFolderName = 'initial_background';


threshold = thresholdbych(experimentDir, experimentName, position, ...
folderArray, typedots, saveImages, backgroundFolderName, numCh, I);
