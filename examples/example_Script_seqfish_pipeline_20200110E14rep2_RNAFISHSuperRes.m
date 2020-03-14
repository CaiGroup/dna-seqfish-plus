%% Variables
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');
experimentDir = 'E:\Yodai\DNAFISH+\2019-07-06-E14-DNA-seqFISH+rep2-2-RNAFISH-IF - Swapped';
experimentName = '2019-07-06-E14-DNA-seqFISH+rep2-2-RNAFISH-SuperRes';
experimentLabel = '2020-01-12-seqential-analysis';
superres = 'radial3d'; % 'radial3d', 'gaussian'
filtersigma = true; % to remove hot pixels.
posArray = 0:4;
folderArray = 0:17;
numCh = 2;
chArray = 1:numCh;
usechabboffsets = true;
usephysicaloffsets = false;
load('E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\analysis\2020-01-10-pipeline-test\points\variables\refpoints-chromaticaberrations-initial-final-2020-01-10.mat','chaTformFinal');
chaTform = chaTformFinal;
segPath = 'E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\segmentation';


%% save global variables
saveGlobalDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'variables');
if exist(saveGlobalDir, 'dir') ~= 7
    mkdir(saveGlobalDir);
end
dateStart = datetime;
formatDate = 'yyyy-mm-dd';
dateString = datestr(dateStart, formatDate);
saveGlobalName = ['global-variables-seqfish-pipeline_' experimentLabel '_' dateString];
saveGlobalPath = fullfile(saveGlobalDir, saveGlobalName);
%save(saveGlobalPath, 'experimentDir', 'experimentName', 'experimentLabel', 'superres', ...
%    'thresholdadjust', 'posArray', 'folderArray', 'numCh', 'sizeZ', 'numRounds', ...
%    'numChannels', 'sqrtradius', 'alloweddiff', 'usechabboffsets', 'usephysicaloffsets', ...
%    'filtersigma', 'numRefPoints', 'medRefIntensity', 'longestEdges', 'minEdgeMatch', ...
%    'minDotMatch');
       
for position = posArray
    listing = dir(['E:\Yodai\DNAFISH+\2019-07-06-E14-DNA-seqFISH+rep2-2-RNAFISH-IF - Swapped\processedimages\pos' num2str(position) '\preProcessedData*.mat']);
    load([listing(1).folder '\' listing(1).name]);
    for ch = chArray
        [points, pointsch, pointschpercell, intensity] ...
            = seqfishprocess_preprocessedimages(experimentDir, experimentName, experimentLabel, ...
            position, folderArray, superres, chaTform{ch} , usechabboffsets, [], usephysicaloffsets, filtersigma, I, segPath);
    end

end