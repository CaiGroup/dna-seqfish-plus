function [statsFilt, boundaries, shape2d] = getsegmentationavgintdapi(L, rawImage, position, saveDir, experimentName)
% function gets the cellinfo using the multicut and simple segmentations
%
% L = cellnum; % or labeled image
%
% Date: 11/1/2019
% Author: Nico Pierson
% Email: nicogpt@caltech.edu

    %% Assign to the cell matrix using the segmentation
    stats = regionprops3(L,'Volume','Centroid','BoundingBox','VoxelList');
    % filter out the cells that don't have an area of 5000 or more
    %minCellVolume = 10000; % About 70x70 px for a 2d slice * 20 = 10,000
    %maxCellVolume = 150000;% 250*250*25 = 1,500,000
    %keepIndex = find(stats.Volume > minCellVolume & stats.Volume < maxCellVolume);
    statsFilt = stats;%stats(keepIndex, :);

    %% Get the Boundaries of the cell
    fprintf('Start Getting Boundaries of the Cells, Position %.0f...\n', position);
    numCells = size(statsFilt, 1);
    statsVoxelList = statsFilt.VoxelList;
    numZSlice = size(rawImage, 3);
    saveVideoPath = fullfile(saveDir, ['dapi-seg-' experimentName '-pos' num2str(position)]);
    [boundaries, shape2d] = makevideovoxellistalphashapeavgint(statsVoxelList, position, numZSlice, saveVideoPath, rawImage);


    numCellFields = 3;
    cellinfo = cell(numCellFields + 1, size(statsFilt,1) + 1);
    cellinfo{2,1} = 'Volume';
    cellinfo{3,1} = 'Centroid';
    cellinfo{4,1} = 'BoundingBox';
    for col = 1:size(statsFilt,1)
        cellinfo{1, col+1} = ['cell ' num2str(col)];
    end
    statsFinal = statsFilt;
    statsFinal.VoxelList = [];
    insertCellInfo = table2cell(statsFinal)';
    cellinfo(2:end,2:end) = insertCellInfo(1:numCellFields,:);

    %% save spots that will be used on the cluster
    saveFileName = ['boundaries3d-Data-pos' num2str(position)];
    saveDataHPC = fullfile(saveDir, saveFileName);
    save(saveDataHPC, 'numCells', 'cellinfo', 'statsFilt', 'boundaries', 'shape2d');
            
end