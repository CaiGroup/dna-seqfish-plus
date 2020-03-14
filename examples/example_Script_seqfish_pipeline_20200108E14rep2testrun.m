% Example Script to run the seqFISH pipleline
% Tested on E14-rep2-2
%
% To Do:
% 1. threshold function is specific to grab column 'Pos0'
% 2. add option to use 3d segmentation
% 3. update parameters for bead alignment
% 4. debug brain dataset
%
% Date: 12/19/2019


%% Variables
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');
experimentDir = 'E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentLabel = '2020-01-12-pipeline-test';
superres = 'radial3d'; % 'radial3d', 'gaussian'
thresholdadjust = true; %false; % option to adjust threshold
posArray = 0:4; %0:4;
posArrayAll = 0:4;
folderArray = 0:79;
numCh = 2;
chArray = 1:numCh;
sizeZ = 25;
numRounds = 5;
numChannels = 16;
sqrtradius = 6;
minseeds = 3;
alloweddiff = 2;
usechabboffsets = true;
usephysicaloffsets = true;
filtersigma = true;
longestEdges = 20; 
minEdgeMatch = 3; 
minDotMatch = 4;   
numRefPoints = [];
medRefIntensity = [];
refposition = 0;
segment = 'roi';


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
    'thresholdadjust', 'posArray', 'folderArray', 'numCh', 'sizeZ', 'numRounds', ...
    'numChannels', 'sqrtradius', 'alloweddiff', 'usechabboffsets', 'usephysicaloffsets', ...
    'filtersigma', 'numRefPoints', 'medRefIntensity', 'longestEdges', 'minEdgeMatch', ...
    'minDotMatch');
       


%% Step1: Chromatic aberration and physical offsets
[chaTform, ~, physicalTform, refPointsAlignedInitial, refInitialSaveName] ...
    = debugbeads(experimentDir, posArrayAll, numCh, sizeZ, superres, experimentLabel);

for position = posArray
    
    %% Step2: Process - grab the raw points from preprocessed images.
    [rawpoints, intensity, pointsCsvName, I] ...
        = seqfishprocess_initialpoints(experimentDir, experimentName, experimentLabel, ...
        position, folderArray, superres, chaTform, filtersigma);

    
    %% Step3: Compute tforms for all hybs with fiducial alignment and output aligned/corrected points and offsets
    [pointsch, offsets] = seqfishtformalign(experimentDir, experimentLabel, ...
        position, numRounds, numCh, folderArray, chaTform, physicalTform, ...
        usechabboffsets, usephysicaloffsets, refInitialSaveName{position+1}, pointsCsvName, longestEdges, minEdgeMatch, minDotMatch);
    
    
    %% Step4: Calculate the adjusted thresholding values relative to Pos0, using aligned points within ROIs...
    %         and then detect dots again with updated thresholding values,
    %         and then output unaligned points.
    
    % change the format of pointsch
    corrpoints = [];
    for ch = chArray
        for ro = 1:numRounds
            for bch = 1:numChannels
                corrpoints{(ro-1)*numChannels+bch,ch} = pointsch{ch,1}{ro,1}(bch).channels;
            end
        end
    end
    
    if position == refposition %% refposition doesn't need the thresholdinig correction.
    
        % may need to save variables to match the name same as step4,5 output.
        roiPath = fullfile(experimentDir, 'segmentation', ['Pos' num2str(refposition)], 'RoiSet.zip');
        if exist(roiPath, 'dir') ~= 7
            vertex = selfsegzip(roiPath);
            roimask = roi2mask(vertex, [2048 2048 sizeZ]);
        end
        [numPoints, medIntensity] = getnumrefpoints(corrpoints, intensity, chArray, numRounds, numChannels, roimask);
        medRefIntensity = medIntensity;
        numRefPoints = numPoints;
        finalpoints = pointsch;
        savePath = fullfile(saveGlobalDir, ['data-pos' num2str(position) '-' experimentName '-step5-finalpoints-' experimentLabel '.mat']);
        save(savePath, 'finalpoints');
    
    else %% Step4: Other positions  
        [rawpoints2, intensity2, numRefPoints, medRefIntensity] ...
            = seqfishprocess_alignedpoints(experimentDir, experimentName, experimentLabel, ...
            position, folderArray, superres, numChannels, ...
            numRounds, sizeZ, chaTform, numRefPoints, medRefIntensity, filtersigma, corrpoints, intensity, I, refposition);
        % check file saving. it seems same variables are saved many times.

    %% Step5: Align points using existing tforms computed at step2.
        for ch = chArray    
            pointsch2{ch} = alignpoints(folderArray, rawpoints2(:,ch), intensity2(:,ch), offsets{ch}, chaTform{ch}, numRounds, physicalTform{position+1,1}{ch,1}, usechabboffsets, usephysicaloffsets);
        end
        finalpoints = pointsch2;
        savePath = fullfile(saveGlobalDir, ['data-pos' num2str(position) '-' experimentName '-step5-finalpoints-' experimentLabel '.mat']);
        save(savePath, 'finalpoints');
       
    end
    %% Step6: decoding for the barcoding dataset.
    % check where intensity is passed to the final output. probably
    % pointsch
    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, finalpoints{channel}, segment, sqrtradius, alloweddiff, ...
        channel, minseeds, experimentLabel);
    end
end