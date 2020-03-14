function [] = printinitfinalbeaddapi(experimentDir, posArray)
% print the dapi bead images to view if there is any expansion.

initName = 'initial_fiducial_markers';
finalName = 'final_fiducial_markers';
initDir = fullfile(experimentDir, initName);
finalDir = fullfile(experimentDir, finalName);

    for position = posArray
        % get images
        initPath = fullfile(initDir, ['MMStack_Pos' num2str(position) '.ome.tif']);
        [initI, isizeC, isizeZ, ~, ~] = grabimseries(initPath, position);
        finalPath = fullfile(finalDir, ['MMStack_Pos' num2str(position) '.ome.tif']);
        [finalI, fsizeC, fsizeZ, ~, ~] = grabimseries(finalPath, position);

        % find tform
        tform = grabtform(finalI{fsizeC}, initI{isizeC});
        finalDapi = applydapitform(finalI, tform);

        % print together
        images = cell(2,1);
        images{1} = initI{isizeC};
        images{2} = finalDapi{fsizeC};

        % save the image
        saveDir = initDir;
        startString = 'checkBeadAlignment';
        endingString = '2019-12-05';
        savefolchimage(position, images, saveDir, startString, endingString);


    end




end
