addpath('C:\github\streamline-seqFISH\src\preprocessing\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');

experimentDir = 'M:\Yodai\left_confocal\20180207_Yodai_7_short_E14_serum_10kintron_rep4';
experimentName = '20180207_Yodai_7_short_E14_serum_10kintron_rep4';

folderArray = 0:29;
useBackground = true;
subtractBackground = false;
backgroundFolder = 'Blanks'; % default will be used
dapiRefPath = ''; % use HybCycle_0 in experiment directory
imageJBackSubtract = true;
saveProcessIms = false;
divideIms = false;
dim = '2d';

for position = 0:1
    [I] = preprocessimagesNewFolder(experimentName, experimentDir, position, folderArray, ...
        useBackground, backgroundFolder, dapiRefPath, imageJBackSubtract, ...
        subtractBackground, saveProcessIms,divideIms, dim);
end