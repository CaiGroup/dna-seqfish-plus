% Example Script 
% running on experiment E14L22Rik24hr_RNAFISH with 2 channels
%
% Requirement: need to run the dnafish alignment firstc
%
% Date: 1/08/2019


% add paths examples
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\github\streamline-seqFISH\src\FindThreshold\', '-end');



%% Variables
% directories
experimentDir = 'L:\2020-01-26-E14L22Rik48hr-RNAFISH-IF';
experimentFiduciaryDir = 'L:\2020-01-26-E14L22Rik48hr-RNAFISH-IF';
experimentName = '2020-01-26-E14L22Rik48hr-RNAFISH-IF';
experimentLabel = '02-18-2019-bin2x2x1'; % Label for save folder in 'analysis'
githubDir = 'C:\github'; % set up path for github
dnafishExperimentLabel = '01-10-2019-beadalign'; % label from dnafish experiment
dnafishPath = fullfile(experimentFiduciaryDir, 'analysis', dnafishExperimentLabel); % experiment path from dna fish experiment to align to hyb 1 ch 1 offsets

% immuno variables
segoption = '2d'; % segmentation in '2d' or '3d'
immunoChArray = 1:2; % channels of immuno images
binSize = 2; % bin pixels for intensity average

% number variables
superres = 'radial3d'; % 'radial3d', 'gaussian'
typedots = 'log'; % use log filter to grab spots
posArray = 0:4; % array of positions
posArrayChaTform = 0:4; % array for fiduciary marker images
folderArray = 0:10; % folder array to align, process, decode
numCh = 2; % # of chanels to decode
channelAnalysis = 3; % which channel of ref points to use
numChFiduciary = 3; % # of channels for fiduciary markers
sizeZ = 25; % # of z-slices
numRefPointsCh3 = []; % number of reference points to get for alignment, 1000 for E14 rep2-2

% flag options
filtersigma = true; % filter using the sigma from the radial center algorithm

% alignment
longestEdges = 20; % number of edges
minEdgeMatch = 3;  % minimume number of edges to match
minDotMatch = 4;   % minimum spots to align

% threshold
numRefPoints = []; 
firstThreshold = [];
roimask = [];

imageJBackSubtract = false;



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
save(saveGlobalPath, 'experimentDir', 'experimentName', 'experimentLabel', 'superres', ...
    'posArray', 'folderArray', 'numCh', 'sizeZ', 'segoption', 'typedots', 'numChFiduciary', ...
    'numRefPointsCh3', 'posArrayChaTform', 'immunoChArray', 'dnafishExperimentLabel', 'dnafishPath', ...
    'githubDir', 'firstThreshold', 'roimask', 'binSize', ...
    'filtersigma', 'numRefPoints', 'longestEdges', 'minEdgeMatch', ...
    'minDotMatch');

       

%% get the reference threshold which is usually pos0
% if csv, confirm the column variable for the threshold
thresholdDir = fullfile(experimentDir, 'threshold');
threshold = 100;


%% Chromatic aberration and physical offsets
% load from the dnafish experiment
%chaAbbPath = getfile(fullfile(dnafishPath, 'points', 'variables'), 'refpoints-chromaticaberrations-initial-final', 'match');
%load(chaAbbPath, 'chaTformInitial', 'physicaltform', 'refPointsInitial', 'refInitialSaveName');


[chaTformInitial, ~, physicaltform, refPointsInitial, refInitialSaveName, ~] ...
    = debugbeadsonlyinitial(experimentFiduciaryDir, posArrayChaTform, numCh, sizeZ, superres, dnafishExperimentLabel);


for position = posArray
    
 
    refCh3SaveName = '';
    %% get the points using the autothreshold....and apply to all other hybs
    numRefPoints = 1000;
    [I, points, intensity, shadingcorr, adjustedThreshold] = outputprocessimagesimmunoUseDapiTform2(experimentDir, ...
        position, folderArray, numCh, typedots, superres, firstThreshold, numRefPoints, ...
        experimentName, experimentLabel, roimask, filtersigma, refCh3SaveName, ...
        githubDir, longestEdges, minEdgeMatch, minDotMatch, chaTformInitial, physicaltform, ...
        threshold, dnafishPath, imageJBackSubtract);
    
    

    % use the images to output the alignment
    switch segoption
        case '2d'
            L = [];
        case '3d' % need to debug and see if labeled cells are correct
            labelDir = fullfile(experimentDir, 'segmentation', ['Pos' num2str(position)]);
            labelCellPath = getfile(lableDir, '*.h5', 'match');
            L = geth5mask(labelCellPath);
    end
    
    avgintpercell(experimentDir, experimentName, experimentLabel, I(:,immunoChArray), ...
        L, position, folderArray, immunoChArray, segoption, binSize);

    
    experimentLabel = '02-18-2019-bin1x1x1';
    binSize = 1;
    avgintpercell(experimentDir, experimentName, experimentLabel, I(:,immunoChArray), ...
        L, position, folderArray, immunoChArray, segoption, binSize);
    

    %{
    experimentLabel = '02-03-2019-bin1x1x1-backsubtract';
    binSize = 1;
    uniqueString = 'imageTempProcess-3535fsfsg';
    for f = size(I,1)
        for ch = 1:numCh
            I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString,...
                experimentDir);
        end
    end
    avgintpercell(experimentDir, experimentName, experimentLabel, I(:,immunoChArray), ...
        L, position, folderArray, immunoChArray, segoption, binSize);
        %}
    
    experimentLabel = '02-18-2019-bin2x2x1';
    binSize = 2;
end