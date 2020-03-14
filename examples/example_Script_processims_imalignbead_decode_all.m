% script for processing yodai dna fish rep2 for each channel
%
% Dependences: 
% 1. Packages: preprocessing, processing, AlignImages
% 2. threshold: round by pseudochannel
% 3. barcodekey
% 4. (if necessary) roi segmentations
%
% Date: 8/14/2019



% Variables for processing
experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH';
fovArray = 3:-1:0;
folderArray = 0:80;
numHybs = length(folderArray) - 1; % last is a repeat
backgroundFolder = 'HybCycle_initial_background';
useBackground = false;
saveProcessIms = false;
backgroundSubtract = false;

% Variables for Decoding
numRounds = 5;
numChannels = 16;
superres = 'radial';
channelArrayBarcode = 1:2;
alloweddiff = 2;
typedots = 'exons';
sqrtradius = 6;
savePoints = true;

% Variables fo Bead Alignment
targetNumPoints = 400;
fixedThreshold = [500,450,500];%fixedThreshold = [];%get initial threshold from one position and check if other positions have the same number of points
%hybThreshold = [];
channelArrayBeads = 1:2;

for position = fovArray
    fprintf('Processing Position: %.0f\n\n', position);
    
    %% Process the images
    %[I, hybIms, tformDapi] = preprocessimages(experimentName, experimentDir, position, ...
    %    folderArray, useBackground, backgroundFolder, backgroundSubtract, saveProcessIms);
    
    % directories
    saveDir = fullfile(experimentDir, 'organizehybs', ['pos' num2str(position)]);
    
    %position 3 has processed images
    if position ~= 3
        % get initial dapi hyb
        refImPath = fullfile(experimentDir, 'HybCycle_0', ['MMStack_Pos' num2str(position) '.ome.tif']);
        [refIms, numRefDapi, numRefZSlice, ~, ~] = grabimseries(refImPath, position);
        dapiRef = refIms{numRefDapi};
        numCh = numRefDapi - 1;

        % get shading corrections
        backgroundFolderName = 'HybCycle_initial_background';
        [shadingcorr, backIms, backTform] = alignbackimsgetshading(dapiRef, experimentDir, ...
            backgroundFolderName, position);


        % background subtract for position 1
        fprintf('Loading Raw Images for position %.0f\n', position);
        dataPath = fullfile(saveDir, ['imagesHybDapi-pos' num2str(position) '-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH-2019-08-15.mat']);
        load(dataPath, 'hybIms', 'tform'); tformDapi = tform; % hybIms are aligned dapi images


        % Apply shading corrections
        I = useshadingsubtractback(experimentDir, hybIms, backIms, folderArray, numCh, shadingcorr);
        clearvars hybIms

       % Save Mat Files
        savePath = fullfile(saveDir, ['preProcessedData-pos' num2str(position) '-' experimentName '-2019-08-14.mat']);
        save(savePath, 'I', 'backIms', 'backTform', 'shadingcorr', 'position','-v7.3');%'dapiTform', 'shadingcorr', 'position','-v7.3');
    
    else
        savePath = fullfile(saveDir, ['preProcessedData-pos' num2str(position) '-' experimentName '-2019-08-14.mat']);
        load(savePath, 'I');
        dataPath = fullfile(saveDir, ['imagesHybDapi-pos' num2str(position) '-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH-2019-08-15.mat']);
        load(dataPath, 'tform'); tformDapi = tform;
    end
    
    
    %% Get all raw images
    rawHybIms = [];
    sizeC = [];
    for f = folderArray
        
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(experimentDir, ['HybCycle_' num2str(f)], imageName);
        [rawHybImsTemp, sizeC, ~, ~, ~] = grabimseries(imagePath, position);
        rawHybIms = cat(1, rawHybIms, rawHybImsTemp);
    end
    numHybChannels = sizeC - 1;
    rawHybIms = rawHybIms(:, 1:numHybChannels);
        
    
    
    %% Align by beads
    hybThreshold = [];
    [hybImsBeadAligned, hybThreshold, fixedThreshold, beadMatchPoints] = imalignbeads(experimentName, experimentDir, position, ...
        channelArrayBeads, I, rawHybIms, hybThreshold, fixedThreshold, tformDapi, targetNumPoints);


    % Process the images for each channel
    % Add Bead points and remove any of these points
    for ch = channelArrayBarcode
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, ...
            numfinalpoints, numpointspercell, seeds, points] = ...
            processimages(experimentDir, experimentName, position, ...
            numRounds, numChannels, hybImsBeadAligned(1:numHybs,ch), sqrtradius, typedots, ...
            superres, savePoints, alloweddiff, ch, beadMatchPoints);
    end
        %}
    
end