% function to adjust the threshold for all positions by using a multiplier
% to compare the median intensity of the captured points - use E14 rep3-1
% as the reference

% get the images and points using the reference threshold and calculate the
% median reference intensity across all hybridizations....

% create a table and compare to the reference and adjust the multiplier to
% find a closer median intensity across all positions
%{
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';



%% get the reference threshold
csvPath = getfile(fullfile(experimentDir, 'threshold'), ...
    '002-008-E14-rep2-2-DNAFISH-ch1-2-hyb1-80-threshold-logfilter-20191121.csv', 'match');
threshold = threshcsv2mat(csvPath, chArray);



%% get median intensity from reference points
numRounds = 5;
numChannels = 16;
refposition = 0;
chArray = 1:2;
medianIntensity = outputmedianintensity(experimentDir, refposition, chArray, numRounds, numChannels);



%% Get the Ref beads and chromatic aberrations from the initial and final beads
posArray = 0:4;
beadFolderName = 'initial_fiducial_markers';
[chaTform, chaTformAvgAll, chaTformAll] = outputallchromaticaberration(experimentDir, beadFolderName, posArray, numCh);
%}


%% Loop through each position to get the raw points
posArray = [0 2 3 4];
folderArray = 0:83;
numCh = 3;
numPointChannels = size(threshold,2); % usually it is 3 but the adjusted thresholds Yodai provided only has 2 for E14-rep2-2
typedots = 'log'; % new log filter
for position = posArray
    % Output processed images
    [I, dapiI, shadingcorr] = outputprocessimages(experimentDir, position, folderArray, numCh);
    sizeZ = size(I{1,1},3);
    
    % get adjusted threshold if not the reference position
    if position ~= refposition
        % get roi mask
        roiPath = fullfile(experimentDir, 'segmentation', ['Pos' num2str(position)], 'RoiSet.zip');
        if exist(roiPath, 'file') == 2
            vertex = selfsegzip(roiPath);
            roimask = roi2mask(vertex, [2048 2048 sizeZ]);
        else
            roimask = [];
        end
        % adjusted
        [adjustedThreshold, adjustedIntensity, adjustedPoints] = ...
            outputadjustedthreshold(I, medianIntensity, threshold, folderArray, ...
            numPointChannels, typedots, roimask);
    else
        % same threshold
        adjustedThreshold = threshold;
    end
    
    % output the points and get them, also saves variables
    [points, intensity, pass] = outputrawpoints(I, adjustedThreshold, folderArray, typedots, ...
        chaTform, experimentDir, experimentName, position, numPointChannels, shadingcorr);
    
end



    
    % get the offsets
    
    
    
    
    
%% decode the points
numRounds = 5;
numChannels = 16;
sqrtradius = 6;
segment = 'whole';
minseeds = 3;
alloweddiff = 2;
chArray = 1:2;
typedots = 'log';
superres = 'radial3d';
for postion = posArray
    outputdecodedpoints(experimentDir, experimentName, position, numRounds, ...
        numChannels, sqrtradius, segment, alloweddiff, typedots, superres, chArray, ...
        minseeds, chaTform, folderArray);
end

%% QC check for chromatic aberration corrections
% if the initial fiducial marker chaTforms will line up the final fiducial markers




