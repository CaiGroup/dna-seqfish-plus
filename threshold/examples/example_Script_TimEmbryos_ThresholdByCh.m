load('G:\Tim\102219\Experiment\processedimages\pos0\preProcessedData-pos0-Tim-Embryos-102219-2019-10-30.mat');

folderArray = 0:29;
experimentDir = 'G:\Tim\102219\Experiment';
experimentName = 'Tim-102219';
position = 0;
typedots = 'log';
threshold = thresholdbych(experimentDir, experimentName, position, ...
folderArray, typedots, I);
