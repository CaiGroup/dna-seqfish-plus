function [finalPosList, numCells] = assignspots2cells(roiPathName, spots)
% function assigns spots to cells if there is a set of ROIs for the cells
%
% Input: path to roi zip file, spots data 1xbarcode cell with the positions
%
% Date: 7/26/2019

    %% Variables
    numGenes = length(spots);
    noroi = false;
    if exist(roiPathName, 'file')   
        vertex = selfsegzip(roiPathName);
        numCells = length(vertex);
    else
        numCells = 1;
        noroi = true;
    end
    
    finalPosList = cell(numGenes, numCells + 1);
    if noroi
        finalPosList(:,2) = spots;
    else
        for i = 1:numCells
            for j = 1:numGenes
            if ~isempty(spots{j})
                    include = inpolygon(spots{j}(:,1),spots{j}(:,2),vertex(i).x,vertex(i).y);
                    getspots = spots{j}(include,:);
                    finalPosList{j, i+1} = getspots;
                    spots{j}(include,:) = []; % remove the previous spots
                end
            end
        end
    end

end