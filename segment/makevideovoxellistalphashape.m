function [boundaries, shape2d] = makevideovoxellistalphashape(voxellist, pos, numZSlices, saveVideoFilePath, rawImagePath)
% function that saves the boundaries of a cell using the voxellist
%
%     delete(gcp('nocreate')); % delete the current number of workers
%     parpool('local',10); % to set up the numbwer of workers
%
% Date: June 3, 2019

    % load the rawImages
    numCells = length(voxellist);
    nisslePath = fullfile(rawImagePath, ['pos' num2str(pos) 'nissle.h5']);
    %nisslePath = fullfile(rawImagePath, ['nisslepos' num2str(pos) '.h5']);
    nissleImage = geth5masknissle(nisslePath);
    warning('off','all');

    boundaries = cell(numZSlices,1);
    shape2d = cell(numZSlices, 1);
    parpool('local',4);
    parfor z = 1:numZSlices
        fprintf('Geting Vertices: %.0f\n', z);
        boundaries{z} = cell(numCells, 1);
        shape2d{z} = cell(numCells, 1);
        for c = 1:numCells
            %shp = alphaShape(voxellist{c});
            %shp = alphaShape(voxellist{c}, 1, 'RegionThreshold', 1, 'HoleThreshold', 80);
            zInd = find(voxellist{c}(:,3) == z);
            % assign zInd to each vertices
            drawVertices = voxellist{c}(zInd, :);
            shape2d{z}{c} = alphaShape(drawVertices(:,1:2));
            if shape2d{z}{c}.numRegions > 0
                [~, vertices] = boundaryFacets(shape2d{z}{c},1); % get region 1
            else
                vertices = [];
            end
            boundaries{z}{c} = vertices;
        end
    end
    
    upperImageThreshold = 1800; % get the original image: % 8000 for Control1 %1800 for Control2 %40000 for aggression and mating
    makevideofromboundaries(boundaries, nissleImage, numZSlices, saveVideoFilePath, upperImageThreshold);
    
end