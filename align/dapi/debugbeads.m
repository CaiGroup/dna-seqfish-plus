function [chaTformInitial, chaTformFinal, physicaltform, refPointsAlignedInitial, ...
    refInitialSaveName] = debugbeads(experimentDir, posArray, numChannels, sizeZ, ...
    superres, experimentLabel)
% outputs the ref points for each channel as a matlab figure and save the
% points using the roi to filter points in cells.

    %% Paths
    addpath('C:\github\streamline-seqFISH\src\FindThreshold', '-end');
    
    
    
    %% Get the Ref beads and chromatic aberrations from the initial and final beads
    initFolderName = 'initial_fiducial_markers';
    finalFolderName = 'final_fiducial_markers';
    
    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = datestr(dateStart, formatDate);
    
    
    
    % Set up variables
    manualThreshold = true; %false; %false; % use manual thrseshold or automatic threshold
    numCh = numChannels;
    chaTformArray = posArray;
    chaTformInitialAll = cell(length(chaTformArray),1);
    chaTformFinalAll = chaTformInitialAll;
    refPointsAlignedInitial = chaTformInitialAll;
    refPointsAlignedFinal = chaTformInitialAll;
    intensityInitial = chaTformInitialAll;
    intensityFinal = chaTformInitialAll;
    refPointsInitial = chaTformInitialAll;
    refInitialSaveName = chaTformInitialAll;
    refPointsFinal = chaTformInitialAll;
    chaTformAvgInitial = cell(numCh,1);
    chaTformAvgFinal = cell(numCh, 1);
    chaTformInitial = cell(numCh,1);
    chaTformFinal = cell(numCh, 1);
    numFirstPoints = [];
    firstThreshold = [];

    
    
    %% Get the chaTforms for each of the beads
    for position = chaTformArray
        roiPath = fullfile(experimentDir, 'segmentation', ['Pos' num2str(position)], 'RoiSet.zip');
        if exist(roiPath, 'file') == 2
            vertex = selfsegzip(roiPath);
            roimask = roi2mask(vertex, [2048 2048 sizeZ]);
        end
        
        [chaTformInitialAll{position+1}, refPointsAlignedInitial{position+1}, ...
            intensityInitial{position+1}, firstThreshold, numFirstPoints] = ...
            getrefbeadsfilter(experimentDir, position, initFolderName, roimask, ...
            superres, manualThreshold, experimentLabel, firstThreshold, numFirstPoints); % this is the chromatic aberration corrections for one position
        [chaTformFinalAll{position+1}, refPointsAlignedFinal{position+1}, ...
            intensityFinal{position+1}, firstThreshold, numFirstPoints] = ...
            getrefbeadsfilter(experimentDir, position, finalFolderName, roimask, ...
            superres, manualThreshold, experimentLabel, firstThreshold, numFirstPoints);
        for ch = 1:numCh
            if position == 0
                chaTformAvgInitial{ch} = zeros(4,4,length(chaTformArray));
                chaTformAvgFinal{ch} = zeros(4,4,length(chaTformArray));
            end
            % each channel has a tform....take average
            chaTformAvgInitial{ch}(:,:,position+1) = chaTformInitialAll{position+1}{ch}.T;
            chaTformAvgFinal{ch}(:,:,position+1) = chaTformFinalAll{position+1}{ch}.T;
        end
    end
    for ch = 1:numCh
        chaTformInitial{ch} = affine3d(mean(chaTformAvgInitial{ch},3));
        chaTformFinal{ch} = affine3d(mean(chaTformAvgFinal{ch},3));
    end
    
    
    
    %% Use chaTform on reference beads
    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'pre_formated');
    if exist(saveDir, 'dir') ~= 7 
        mkdir(saveDir);
    end
    for position = chaTformArray
        endingString = [initFolderName '-' experimentLabel];
        [refPointsInitial{position+1}, refInitialSaveName{position+1}] = outputrefpoints(refPointsAlignedInitial{position+1}, chaTformInitial, intensityInitial{position+1}, position, endingString, saveDir, numCh);
        endingString = [finalFolderName '-' experimentLabel];
        [refPointsFinal{position+1}, ~] = outputrefpoints(refPointsAlignedFinal{position+1}, chaTformFinal, intensityFinal{position+1}, position, endingString, saveDir, numCh);
    end
    
    
    
    %% Physical Tform after chromatic aberration offsets applied to reference points
    physicaltform = getphysicaloffsets(posArray, refPointsInitial, ...
        intensityInitial, initFolderName, numCh, chaTformInitial);
    refSaveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'variables');
    if exist(refSaveDir, 'dir') ~= 7
        mkdir(refSaveDir);
    end
    save(fullfile(refSaveDir, ['refpoints-chromaticaberrations-initial-final-' endingDateString '.mat']), ...
        'refPointsInitial', 'refPointsFinal', 'physicaltform', 'intensityInitial', 'intensityFinal', ...
        'chaTformInitial', 'chaTformFinal', 'chaTformInitialAll', 'chaTformFinalAll', 'chaTformArray', ...
        'refInitialSaveName');

end