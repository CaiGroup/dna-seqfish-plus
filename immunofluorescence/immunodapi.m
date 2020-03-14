function [chaTform, initialThreshold] = immunodapi(experimentDir, experimentName, experimentLabel, ...
    beadsFolderName, posArray, sizeZ, numChannels, superres, binSize, searchradius, dnaFishDir, chaPosArray)

% used to grab the reference dapi, get offset from bead images from ch4 dapi
% to channel 1, apply offset to dapi, calculate immunofluorescence.
%
% date: 1/14/2019


% outputs the ref points for each channel as a matlab figure and save the
% points using the roi to filter points in cells.

    %% Paths
    addpath('C:\github\streamline-seqFISH\src\FindThreshold', '-end');
    
    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = datestr(dateStart, formatDate);
    
    
    
    % Set up variables
    manualThreshold = true; %false; %false; % use manual thrseshold or automatic threshold
    numCh = numChannels;
    chaTformArray = chaPosArray;
    chaTformInitialAll = cell(length(chaTformArray),1);
    refPointsAlignedInitial = chaTformInitialAll;
    intensityInitial = chaTformInitialAll;
    refPointsInitial = chaTformInitialAll;
    refInitialSaveName = chaTformInitialAll;
    chaTformAvgInitial = cell(numCh,1);
    chaTformInitial = cell(numCh,1);
    numFirstPoints = [];
    firstThreshold = [];
    L = [];
    segoption = '2d';
    folderArray = 0;

    dapiSaveDir = fullfile(experimentDir, 'analysis', experimentLabel);
    if exist(dapiSaveDir, 'dir') ~= 7
        mkdir(dapiSaveDir);
    end
    
    dapiRef = cell(length(posArray),1);
    for position = posArray
        %% Load the dapi ref images for each position in HybCycle0
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(dnaFishDir, ['HybCycle_' num2str(0)], imageName);
        [allIms, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
        dapiRef{position+1} = allIms{sizeC};
    end
    
    
    
    %% Get the chaTforms for each of the beads
    for position = chaTformArray
        roiPath = [];
        if exist(roiPath, 'file') == 2
            vertex = selfsegzip(roiPath);
            roimask = roi2mask(vertex, [2048 2048 sizeZ]);
        else 
            roimask = [];
        end
        
        [chaTformInitialAll{position+1}, refPointsAlignedInitial{position+1}, ...
            intensityInitial{position+1}, firstThreshold, numFirstPoints] = ...
            getrefbeadsfilter(experimentDir, position, beadsFolderName, roimask, ...
            superres, manualThreshold, experimentLabel, firstThreshold, numFirstPoints, ...
            numCh, searchradius, false); % this is the chromatic aberration corrections for one position

        for ch = 1:numCh
            if position == 0
                chaTformAvgInitial{ch} = zeros(4,4,length(chaTformArray));
            end
            % each channel has a tform....take average
            chaTformAvgInitial{ch}(:,:,position+1) = chaTformInitialAll{position+1}{ch}.T;
        end
    end
    for ch = 1:numCh
        chaTformInitial{ch} = affine3d(mean(chaTformAvgInitial{ch},3));
    end
    initialThreshold = firstThreshold;
    dapiRefAligned = cell(1, length(posArray));
    for position = posArray
        dapiRefAligned{position+1} = imwarp(dapiRef{position+1}, chaTformInitial{numCh}, 'OutputView', imref3d(size(dapiRef{position+1})));

        % output the avg intensity csv
        dapiTemp = cell(1, numCh);
        dapiTemp{numCh} = dapiRefAligned{position+1};
        avgintpercell(experimentDir, experimentName, experimentLabel, dapiTemp, ...
            L, position, folderArray, numCh, segoption, binSize);
        
        experimentLabelPrev = experimentLabel;
        experimentLabel = '02-06-2019-bin2x2x1';
        binSize = 2;
        avgintpercell(experimentDir, experimentName, experimentLabel, dapiTemp, ...
            L, position, folderArray, numCh, segoption, binSize);
        
        experimentLabel = experimentLabelPrev;
        binSize = 1;
    end
    
    
    %% Output raw reference fiduciary markers
    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'pre_formated');
    if exist(saveDir, 'dir') ~= 7 
        mkdir(saveDir);
    end
    
    
    
    %% Physical Tform after chromatic aberration offsets applied to reference points
    refSaveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'beads-tforms');
    if exist(refSaveDir, 'dir') ~= 7
        mkdir(refSaveDir);
    end
    save(fullfile(refSaveDir, ['beads-chromaticaberrations-dapiref' endingDateString '.mat']), ...
        'refPointsInitial', 'intensityInitial', 'chaTformInitial', 'chaTformInitialAll', 'chaTformArray', ...
        'refInitialSaveName');
    
    chaTform = chaTformInitial;

end