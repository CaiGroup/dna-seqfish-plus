% main script to run DNA seqFISH+ analysis.
% ch1-2 and ch3 analysis are integrated.
% edit by Yodai Takei
% Last updated on 01/13/20.

%% Experiment dependent Variables
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');
experimentDir = 'E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH';
experimentLabel = '2020-01-14-pipeline-test-newalignment-ch1-3';
posArray = 0:4; % for decoding. needs Pos0 anyway. Don't start from Pos1 or later.
posArrayAll = 0:4;% for fiducial markers.
sizeZ = 25; % number of z slices. in the future, this should be automatically calculted.

% barcoding parameters. may need to be changed and compared.
sqrtradius = 6;
minseeds = 3;
alloweddiff = 2;

%% Other fixed variables for DNA seqFISH+. Some may still need to be changed to compare to.
numChAll = 3;
chArrayAll = 1:numChAll;
numCh = 2; % for barcoding analysis
chArray = 1:numCh;
folderArray = 0:79;
numRounds = 5; % barcoding rounds
numChannels = 16; % barcoding pseudo-channels.
superres = 'radial3d'; % 'radial3d', 'gaussian'
thresholdadjust = true; %false; % option to adjust threshold
usechabboffsets = true; % should be always true.
usephysicaloffsets = true; % default is true.
filtersigma = true;
longestEdges = 20; 
minEdgeMatch = 3; 
minDotMatch = 4;   
numRefPoints = [];
medRefIntensity = [];
refposition = 0;
segment = 'roi';
saveprocessedImages = true;
numCh3pointshyb = 60;

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
    'thresholdadjust', 'posArray', 'folderArray', 'numChAll','numCh', 'sizeZ', 'numRounds', ...
    'numChannels', 'sqrtradius', 'alloweddiff', 'usechabboffsets', 'usephysicaloffsets', ...
    'filtersigma', 'numRefPoints', 'medRefIntensity', 'longestEdges', 'minEdgeMatch', ...
    'minDotMatch');


%% Step1-6: for DNA seqFISH+ analysis.
%% Step1: Chromatic aberration and physical offsets
[chaTform, ~, physicalTform, refPointsAlignedInitial, refInitialSaveName] ...
    = debugbeads(experimentDir, posArrayAll, numChAll, sizeZ, superres, experimentLabel);


for position = posArray
    
    %% Step2: Process - grab the raw points from preprocessed images.
    [rawpoints, intensity, pointsCsvName, I] ...
        = seqfishprocess_initialpoints(experimentDir, experimentName, experimentLabel, ...
        position, folderArray, superres, chaTform, filtersigma, numChAll);

    
    %% Step3: Compute tforms for all hybs with fiducial alignment and output aligned/corrected points and offsets
    [pointsch, offsets, corrpoints] = seqfishtformalign(experimentDir, experimentLabel, ...
        position, numRounds, numChAll, folderArray, chaTform, physicalTform, ...
        usechabboffsets, usephysicaloffsets, refInitialSaveName{position+1}, pointsCsvName, longestEdges, minEdgeMatch, minDotMatch);
    % pointsch: already barcoded cell format
    % corrpoints: sequential cell format, cell(numFolders, numCh)
    % offsets: reference image/points are hyb1 ch1.
    
    %% Step3.2 Output images I if necessary.
    I3paint = seqfish_outputpreprocessedimages(experimentDir, experimentName, experimentLabel, I, offsets, chaTform, physicalTform{position+1}, position, usechabboffsets, usephysicaloffsets, saveprocessedImages);
    
    %% Step4 and step5: Calculate the adjusted thresholding values relative to Pos0, using aligned points within ROIs and then detect dots again with updated thresholding values, and then output unaligned points.
    if position == refposition %% refposition doesn't need any thresholdinig adjustment.
        roiPath = fullfile(experimentDir, 'segmentation', ['Pos' num2str(refposition)], 'RoiSet.zip');
        if exist(roiPath, 'dir') ~= 7
            vertex = selfsegzip(roiPath);
            roimask = roi2mask(vertex, [2048 2048 sizeZ]);
        end
        % calculate medIntensity with aligned points within ROIs.
        [numPoints, medIntensity] = getnumrefpoints(corrpoints, intensity, chArrayAll, numRounds, numChannels, roimask);
        medRefIntensity = medIntensity;
        numRefPoints = numPoints;
        finalpoints = pointsch;
        savePath = fullfile(saveGlobalDir, ['data-pos' num2str(position) '-' experimentName '-step5-finalpoints-' experimentLabel '.mat']);
        save(savePath, 'finalpoints');
    
        
    else %% Step4: Other positions
        % grab points again with updated thresholding values, using initial
        % points intensity info of the position vs Pos0.
        [rawpoints2, intensity2, numRefPoints, medRefIntensity] ...
            = seqfishprocess_alignedpoints(experimentDir, experimentName, experimentLabel, ...
            position, folderArray, superres, numChannels, ...
            numRounds, sizeZ, chaTform, numRefPoints, medRefIntensity, filtersigma, corrpoints, intensity, I, refposition, numChAll);

    %% Step5: Align points using existing tforms computed at step2, and other corrections if necessary.
        for ch = chArrayAll    
            pointsch2{ch} = alignpoints(folderArray, rawpoints2(:,ch), intensity2(:,ch), offsets{ch}, chaTform{ch}, numRounds, physicalTform{position+1,1}{ch,1}, usechabboffsets, usephysicaloffsets);
        end
        finalpoints = pointsch2;
        savePath = fullfile(saveGlobalDir, ['data-pos' num2str(position) '-' experimentName '-step5-finalpoints-' experimentLabel '.mat']);
        save(savePath, 'finalpoints');
       
    end
    
    %% Step6: decoding for the ch1,2 barcoding dataset.
    clearvars I
    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, finalpoints{channel}, segment, sqrtradius, alloweddiff, ...
        channel, minseeds, experimentLabel);
    end
    
    %% Step7: decoding for the ch3 barcoding dataset.
    if numChAll == 3
        pointsch3final = seqfishformatforch3(finalpoints{3}, folderArray, numRounds, numCh3pointshyb);
        ChromPaintIntensities_perpoints_v2(experimentDir, experimentName, experimentLabel, pointsch3final, I3paint, position);
        clearvars I3paint
    end
    
end