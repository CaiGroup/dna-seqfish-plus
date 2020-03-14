function [boundaries, shape2d] = visualize3dsegmentation(saveDir, label, position, rawImage, cellnum)
% makes avi video of the labeled cells
%
% Variables:
% 'cellnum' is a labeled matrix
% 'label' is the name of the file to be saved
% 'saveDir' is the directory to save the video
% 'rawImage' is the raw image to overlay the segmentation boundaries


    % add path to bfmatlab package to load images
    addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');

    % variables
    numZSlices = size(rawImage, 3);
    stats = regionprops3(cellnum, 'VoxelList'); % get cells from labeled image
    voxellist = stats.VoxelList; % get the voxel list

    % save the video
    saveVideoFilePath = fullfile(saveDir, [label '-pos' num2str(position)]);
    [boundaries, shape2d] = makevideovoxellistalphashapeavgint(voxellist, position, numZSlices, saveVideoFilePath, rawImage);
end