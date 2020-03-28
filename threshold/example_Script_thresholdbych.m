% example to use thresholdbych.m

% load image
load('G:\image.mat');

folderArray = 0:29;
experimentDir = 'G:\Tim\102219\Experiment';
experimentName = 'Tim-102219';
position = 0;
typedots = 'log';
threshold = thresholdbych(experimentDir, experimentName, position, ...
folderArray, typedots, I);
