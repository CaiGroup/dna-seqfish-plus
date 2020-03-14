% Example Script to run the seqFISH pipleline on  immuno fluorescence ch1
% and 2, using ch3 as the aligning channel
% Tested on E14-rep2-2
%
% Requirement: need to run the dnafish alignment first
%
% Date: 1/08/2019


% add paths examples
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\github\streamline-seqFISH\src\FindThreshold\', '-end');



%% Variables
% directories
experimentDir = 'I:\2019-07-25-E14-DNA-seqFISH+rep2-2-DNAFISH-plate2';
experimentFiduciaryDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentLabel = '01-22-2019-test-fiduciary-align'; % Label for save folder in 'analysis'
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
folderArray = 0:20; % folder array to align, process, decode
numCh = 3; % # of chanels to decode
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
threshold = loadthreshold(thresholdDir, channelAnalysis); 


%% Chromatic aberration and physical offsets
% load from the dnafish experiment
%chaAbbPath = getfile(fullfile(dnafishPath, 'points', 'variables'), 'refpoints-chromaticaberrations-initial-final', 'match');
%load(chaAbbPath, 'chaTformInitial', 'physicaltform',jk 'refPointsInitial', 'refInitialSaveName');
[chaTformInitial, ~, physicaltform, refPointsInitial, refInitialSaveName, ~] ...
    = debugbeads(experimentFiduciaryDir, posArrayChaTform, numCh, sizeZ, superres, dnafishExperimentLabel);


for position = posArray
    

    %% make ch3 ref points into ch1 that will be used for decoding
    refPath = fullfile(experimentFiduciaryDir, 'analysis', dnafishExperimentLabel, 'points', 'pre_formated', refInitialSaveName{position+1});
    refCh3SaveName = printch3pointsch1(refPath, experimentFiduciaryDir, dnafishExperimentLabel, ...
        position);
    basePointDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points');
    preformatDir = fullfile(basePointDir, 'pre_formated');
    if exist(preformatDir, 'dir') ~= 7 
        mkdir(preformatDir);
    end
    positionsDir = fullfile(basePointDir, 'positions');
    if exist(positionsDir, 'dir') ~= 7
        mkdir(positionsDir);
    end
    refCh3Path = fullfile(experimentFiduciaryDir, 'analysis', dnafishExperimentLabel, 'points', 'pre_formated', refCh3SaveName);
    refCh3SavePath = fullfile(preformatDir, refCh3SaveName);
    copyfile(refCh3Path, refCh3SavePath);
    

    %% get the points using the autothreshold....and apply to all other hybs
    numRefPoints = size(refPointsInitial{position+1}(channelAnalysis).channels,1);
    if isempty(numRefPointsCh3)
        numRefPoints = numRefPoints * 30; % want 1000 for E14 data...need to test on brain
    else
        numRefPoints = numRefPointsCh3;
    end
    [I, points, intensity, shadingcorr, adjustedThreshold] = outputprocessimagesimmunoUseDapiTform(experimentDir, ...
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
   

    experimentLabel = '01-29-2020-1bin-backsubtract';
    binSize = 1;
    avgintpercell(experimentDir, experimentName, experimentLabel, I(:,immunoChArray), ...
        L, position, folderArray, immunoChArray, segoption, binSize);
end