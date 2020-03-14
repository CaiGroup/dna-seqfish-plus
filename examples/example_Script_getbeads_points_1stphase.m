% function to adjust the threshold for all positions by using a multiplier
% to compare the median intensity of the captured points - use E14 rep3-1
% as the reference

% get the images and points using the reference threshold and calculate the
% median reference intensity across all hybridizations....

% create a table and compare to the reference and adjust the multiplier to
% find a closer median intensity across all positions

addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
%% Variables
experimentLabel = '12-18-2019-pipelinetest';
superres = 'gaussian'; % 'radial3d', 'gaussian'
thresholdadjust = true; %false; % option to adjust threshold
posArray = 0;%[0 1 2 3]; %0:4;
folderArray = 0:79;
numCh = 2;
 % usually it is 3 but the adjusted thresholds Yodai provided only has 2 for E14-rep2-2
typedots = 'log';%'exons';% 'log'; % new log filter
experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';



%% get the reference threshold
chArray = 1:2;
%{
csvPath = getfile(fullfile(experimentDir, 'threshold'), ...
    '002-008-E14-rep2-2-DNAFISH-ch1-2-hyb1-80-threshold-logfilter-20191121.csv', 'match');
threshold = threshcsv2mat(csvPath, chArray);
%}
% load yodai's thresholds
csvPath = fullfile(experimentDir, 'threshold', '002-016-E14-rep2-2-DNAFISH-ch1-2-hyb1-80-threshold-logfilter-Pos0vs3-20191208.csv');
threshold = threshcsv2matyodai(csvPath, chArray);



%% get median intensity from reference points
numRounds = 5;
numChannels = 16;
refposition = 0;
chArray = 1:2;
sizeZ = 25;
%{
% fix teh pointsDir and points Name ....got median intensity manually first
% time
pointsDir = fullfile(experimentDir, 'analysis');
pointsName = ['points-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped-Pos0-2019-11-05'];
[medianIntensity, numRefPoints] = outputmedianintensity(pointsDir, pointsName, chArray, numRounds, numChannels);
%}
%% Roi Mask
roiPath = fullfile(experimentDir, 'segmentation', ['Pos' num2str(refposition)], 'RoiSet.zip');
if exist(roiPath, 'dir') ~= 7
    vertex = selfsegzip(roiPath);
    roimask = roi2mask(vertex, [2048 2048 sizeZ]);
end
pointsDir = fullfile(experimentDir, 'points');
pointsName = 'points-int-thresh-pos0-UpdatedThreshold-12-10-2019-MedianIntRatio-logdots-Pos3-UseThresholdPos0-debugadjustedThreshold.mat';
[numRefPoints, medRefIntensity] = loadrefpoints(pointsDir, pointsName, chArray, numRounds, numChannels, roimask);

%{
%% Get the Ref beads and chromatic aberrations from the initial and final beads
posArray = 0:4;
beadFolderName = 'initial_fiducial_markers';
[chaTform, chaTformAvgAll, chaTformAll] = outputallchromaticaberration(experimentDir, beadFolderName, posArray, numCh);
%}
%{
posArray = 0:4;
superbeadres = 'gaussian';
[chaTformInitial, chaTformFinal, chaTformAvgInitial, chaTformAvgFinal, physicalTform, refPoints] = debugbeads(experimentDir, posArray, numCh, sizeZ, superbeadres, saveEnding);
%}
%% load chatforms
load('I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points\pre_formated\refpoints-chromaticaberrations-initial-final.mat', 'chaTformInitial');
chaTform = chaTformInitial;

% load chaTorm 
%load('I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points\pre_formated\refpoints-chromaticaberrations.mat', 'chaTform');
% median intensity from updated log filter thresholds
%load('I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points\pre_formated\medianIntensity-E14-rep2-2-pos0-radial3d=log.mat', 'medianIntensity');
% get the number of points as well to use as a multiplier
refThreshold = threshold;




%% Loop through each position to get the raw points

for position = posArray
    
    %% get median intensity from older points
    %[threshold, medianIntensity] = getthreshmedianintfrompoints(experimentDir, ...
    %    position, numCh, numRounds, numChannels);
    numPointChannels = size(threshold,2);
    
    %% Roi Mask
    roiPath = fullfile(experimentDir, 'segmentation', ['Pos' num2str(position)], 'RoiSet.zip');
    if exist(roiPath, 'dir') ~= 7
        vertex = selfsegzip(roiPath);
        roimask = roi2mask(vertex, [2048 2048 sizeZ]);
    end
    
    %% Output points form the images
    [points, intensity, shadingcorr, adjustedThreshold] = outputprocessimages(experimentDir, ...
        position, folderArray, numCh, typedots, superres, medRefIntensity, numRefPoints, threshold, ...
        refposition, experimentLabel, roimask, thresholdadjust);
    
    %% Get the number of points and the median intensity to save using the roi masks as a filter
    [numPoints, medIntensity] = getnumrefpoints(points, intensity, chArray, numRounds, numChannels, roimask);
    % save the data
    saveDir = fullfile(experimentDir, 'points', 'variables');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    savePath = fullfile(saveDir, ['data-pos' num2str(position) '-' experimentName '-' experimentLabel '.mat']);
    save(savePath, 'points', 'intensity', 'shadingcorr', 'refThreshold', 'adjustedThreshold', ...
        'numRefPoints', 'numPoints', 'medIntensity', 'medRefIntensity', 'roimask');
  
    
    %% output the points and get them, also saves variables
    [points, intensity, pass] = outputrawpoints(points, intensity, adjustedThreshold, ...
        folderArray, chaTform, experimentDir, experimentName, position, ...
        numPointChannels, shadingcorr, experimentLabel);
    
end

