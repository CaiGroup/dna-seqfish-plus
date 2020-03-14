% script to run the immunofluorescence for the dapi
% add paths examples
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\github\streamline-seqFISH\src\FindThreshold\', '-end');

experimentDir = 'K:\2019-11-01-E14-L22Rik-48hr-DNAFISHIF-plate2';

experimentName = '2019-11-01-E14-L22Rik-48hr-DNAFISHIF-plate2';
experimentLabel = '02-06-2019-bin1x1x1';
beadsFolderName = '2019-11-11-250nm-TetraSpeckBeads';
posArray = 0:7; % position array for images
chaPosArray = 0:3; % positions for the beads
z = 29;
numChannels = 4;
superres = 'radial3d';
binSize = 1; % will do both 1 and 2 bins currently/ change name of folder
searchradius = 3;
dnaFishDir = 'E:\2019-10-29-E14-L22Rik-48hr-DNAFISHIF';

[chaTform, initialThresholdTest] = immunodapi(experimentDir, experimentName, experimentLabel, ...
beadsFolderName, posArray, z, numChannels, superres, binSize, searchradius, dnaFishDir, chaPosArray);