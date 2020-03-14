function [statsFilt, boundaries, shape2d] = getsegmentation(segPath, multiName, simpleName, position, savePath, experimentName, stainName)
% function gets the cellinfo using the multicut and simple segmentations
%
% Date: 7/26/2019
% Author: Nico Pierson
% Email: nicogpt@caltech.edu

    %% Assign to the cell matrix using the segmentation
    multiSegPath = fullfile(segPath, multiName);
    simpleSegPath = fullfile(segPath, simpleName);
    simpleMask = geth5mask(simpleSegPath);
    multicutMask = geth5mask(multiSegPath);
    segmentIm = subtractmulticut(multicutMask, simpleMask);
    

    stats = regionprops3(segmentIm,'Volume','Centroid','BoundingBox','VoxelList');
    % filter out the cells that don't have an area of 5000 or more
    %minCellVolume = 100000; % About 70x70 px for a 2d slice * 20 = 100,000
    %maxCellVolume = 1500000;% 250*250*25 = 1,500,000
    %keepIndex = find(stats.Volume > minCellVolume & stats.Volume < maxCellVolume);
    statsFilt = stats;%stats(keepIndex, :);

    %% Get the Boundaries of the cell
    fprintf('Start Getting Boundaries of the Cells, Position %.0f...\n', position);
    numCells = size(statsFilt, 1);
    statsVoxelList = statsFilt.VoxelList;
    numZSlice = size(segmentIm, 3);
    saveVideoPath = fullfile(savePath, [stainName '-seg-' experimentName '-pos' num2str(position)]);
    %%%%% Change function makevideovoxellistalphashape to accept different
    %%%%% stains as the baseimage
    [boundaries, shape2d] = makevideovoxellistalphashape(statsVoxelList, position, numZSlice, saveVideoPath, segPath);


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
    saveDataHPC = fullfile(savePath, saveFileName);
    save(saveDataHPC, 'numCells', 'statsVoxelList', 'cellinfo', 'statsFilt', 'boundaries', 'shape2d');
            
end