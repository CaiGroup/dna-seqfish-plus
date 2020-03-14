function [chaTform, chaTformAvgAll, chaTformAll] = outputallchromaticaberration(experimentDir, beadFolderName, posArray, numChannels)
% outputs the positions for all beads across all positions


    %% Get the Ref beads and chromatic aberrations from the initial and final beads
    %beadFolderName = 'initial_fiducial_markers';
    numCh = numChannels;
    chaTformArray = posArray;
    chaTformAll = cell(length(chaTformArray),1);
    refPointsAlignedAll = chaTformAll;
    intensity = chaTformAll;
    refPoints = chaTformAll;
    chaTformAvgAll = cell(numCh,1);
    chaTform = cell(numCh,1);
    % have pos ch tforms
    for position = chaTformArray
        [chaTformAll{position+1}, refPointsAlignedAll{position+1}, intensity{position+1}] = getrefbeads(experimentDir, position, beadFolderName); % this is the chromatic aberration corrections for one position
        for ch = 1:numCh
            if position == 0
                chaTformAvgAll{ch} = zeros(4,4,length(chaTformArray));
            end
            % each channel has a tform....take average
            chaTformAvgAll{ch}(:,:,position+1) = chaTformAll{position+1}{ch}.T;
        end
    end
    for ch = 1:numCh
        chaTform{ch} = affine3d(mean(chaTformAvgAll{ch},3));
    end
    %% Use chaTform on reference beads
    saveDir = fullfile(experimentDir, 'points', 'pre_formated');
    if exist(saveDir, 'dir') ~= 7 
        mkdir(saveDir);
    end
    for position = chaTformArray
        % use the refPointsAlignedAll
        endingString = ['radial3d-' beadfoldername '-raw-intensity'];
        refPoints{position+1} = outputrefpoints(refPointsAlignedAll{position+1}, chaTform, intensity{position+1}, position, endingString, saveDir, numCh);
    end
    save(fullfile(saveDir, 'refpoints-chromaticaberrations.mat'), 'refPoints', 'intensity', 'chaTform');

end