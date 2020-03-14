%% Variables
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');
experimentDir = 'F:\Yodai\DNA+\2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH-IF - Swapped';
experimentName = '2019-07-16-E14-DNA-seqFISH+rep3-1-RNAFISH-SuperRes';
experimentLabel = '2020-01-28-seqential-analysis-logradial3d';
superres = 'radial3d'; % 'radial3d', 'gaussian'
filtersigma = true; % to remove hot pixels.
posArray = 0:5;
folderArray = 0:20;
numCh = 2;
chArray = 1:numCh;
usechabboffsets = true;
usephysicaloffsets = false;
load('F:\Yodai\DNA+\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\analysis\2020-01-28-full-decoding\points\variables\refpoints-chromaticaberrations-initial-final-2020-01-28.mat','chaTformFinal');
chaTform = chaTformFinal;
segPath = 'F:\Yodai\DNA+\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\segmentation';


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
    % update the following line per experiment
    listing = dir([experimentDir '\processedimages\pos' num2str(position) '\preProcessedData*.mat']);
    load([listing(1).folder '\' listing(1).name]);
    for ch = chArray
        [points, pointsch, pointschpercell, intensity] ...
            = seqfishprocess_preprocessedimages(experimentDir, experimentName, experimentLabel, ...
            position, folderArray, superres, chaTform{ch} , usechabboffsets, [], usephysicaloffsets, filtersigma, I, segPath);
    end

end