function [boundaries, shape2d] = makevideovoxellistalphashapeavgint(voxellist, pos, numZSlices, saveVideoFilePath, rawImage)
% function that saves the boundaries of a cell using the voxellist
%
%     delete(gcp('nocreate')); % delete the current number of workers
%     parpool('local',10); % to set up the numbwer of workers
%
% Date: June 3, 2019

    % load the rawImages
    numCells = length(voxellist);
    warning('off','all');

    boundaries = cell(numZSlices,1);
    shape2d = cell(numZSlices, 1);
    delete(gcp('nocreate'));
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
    
    upperImageThreshold = 2000; % get the original image: % 8000 for Control1 %1800 for Control2 %40000 for aggression and mating
    makevideofromboundaries(boundaries, rawImage, numZSlices, saveVideoFilePath, upperImageThreshold);
    
end