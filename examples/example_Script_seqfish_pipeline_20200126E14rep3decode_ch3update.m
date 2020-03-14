% main script to run DNA seqFISH+ analysis.
% Second time to run the ch3 with different thresholding values.
% threshold values were updated to the new csv file for ch3.
% pull offset, chaTform, physicalTform and I3paint from the previous analysis.
% edit by Yodai Takei
% Last updated on 01/31/20.

%% Experiment dependent Variables
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');
experimentDir = 'F:\Yodai\DNA+\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH';
experimentName = '2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH';
experimentLabel = '2020-01-29-full-decoding';
posArray = 0:5; % for decoding. needs Pos0 anyway. Don't start from Pos1 or later.
posArrayAll = 0:5;% for fiducial markers.
sizeZ = 25; % number of z slices. in the future, this should be automatically calculted.

% barcoding parameters. may need to be changed and compared.
sqrtradius = 3;
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
saveprocessedImages = false;
numCh3pointshyb = 60;

%% Step1-6: for DNA seqFISH+ analysis.
for position = posArray
    %% Step1: load previous variables.
    saveDecodeDir = fullfile(experimentDir, 'analysis', experimentLabel, 'decoding_variables');
    saveDecodePath = fullfile(saveDecodeDir,['decoding-variables-pos' num2str(position) '-' experimentName '.mat']);
    load(saveDecodePath, 'chaTform', 'physicalTform' ,'offsets');
    
    %% Step2: Process - grab the raw points from preprocessed images.
    [rawpoints, intensity, pointsCsvName, I] ...
        = seqfishprocess_initialpoints(experimentDir, experimentName, experimentLabel, ...
        position, folderArray, superres, chaTform, filtersigma, numChAll);

    %% Step3: Compute tforms for all hybs with fiducial alignment and output aligned/corrected points and offsets
    for ch = chArrayAll    
        pointsch{ch} = alignpoints(folderArray, rawpoints(:,ch), intensity(:,ch), offsets{ch}, chaTform{ch}, numRounds, physicalTform{position+1,1}{ch,1}, usechabboffsets, usephysicaloffsets);
    end
    
    numHybs = length(folderArray);
    corrpoints = cell(numHybs,numChAll);
    for ch = chArrayAll % actual imaging channels.
        for ro = 1:numRounds
            for bch = 1:numChannels
                corrpoints{(ro-1)*numChannels+bch,ch} = pointsch{1,ch}{ro,1}(bch).channels;
            end
        end
    end
    
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
       
    end
    
    %% Step6: decoding for the ch1,2 barcoding dataset.
    clearvars I
        % save all the variables for decoding to change the decoding paramters
    % at the next run if necessary.
    saveDecodeDir2 = fullfile(experimentDir, 'analysis', experimentLabel, 'decoding_variables_ch3update20200202');
    if exist(saveDecodeDir2, 'dir') ~= 7
        mkdir(saveDecodeDir2);
    end
    saveDecodePath = fullfile(saveDecodeDir2,['decoding-variables-pos' num2str(position) '-' experimentName '.mat']);
    save(saveDecodePath, 'experimentDir', 'experimentName', 'experimentLabel', 'superres', ...
    'thresholdadjust', 'posArray', 'position','folderArray', 'numChAll','numCh', 'chArray', 'sizeZ', 'numRounds', ...
    'numChannels', 'sqrtradius', 'alloweddiff', 'usechabboffsets', 'usephysicaloffsets', ...
    'filtersigma', 'numRefPoints', 'medRefIntensity', 'longestEdges', 'minEdgeMatch', ...
    'minDotMatch', 'finalpoints', 'segment', 'minseeds', 'numCh3pointshyb', ...
    'chaTform', 'physicalTform' ,'offsets');
    
    %% Step7: decoding for the ch3 barcoding dataset.
    if numChAll == 3
        pointsch3final = seqfishformatforch3(finalpoints{3}, folderArray, numRounds, numCh3pointshyb);
        ChromPaintIntensities_perpoints_v2(experimentDir, experimentName, experimentLabel, pointsch3final, I3paint, position);
        clearvars I3paint
    end
    
end