function pointsincells = filterpointspercell(roiPathName, points)
% function assigns spots to cells if there is a set of ROIs for the cells
%
% Input: path to roi zip file, spots data 1xbarcode cell with the positions
%
% Date: 7/26/2019

    %% Variables
    numRounds = length(points);
    numChannels = length(points{1});
    vertex = selfsegzip(roiPathName);
    numCells = length(vertex);
    pointsincells = cell(1, numRounds);
    

    for j = 1:numRounds
        pointsincells{j} = struct('channels', cell(1, numChannels));
        pointsincells{j} = struct('intensity', cell(1, numChannels));
        pointsincells{j} = struct('scaledIntensity', cell(1, numChannels));
        for k = 1:numChannels
            if ~isempty(points{j}(k).channels)
                for i = 1:numCells
                    include = inpolygon(points{j}(k).channels(:,1),points{j}(k).channels(:,2),vertex(i).x,vertex(i).y);
                    if i == 1
                        includeall = include;
                    end
                    includeall = or(include, includeall);
                end
                pointsincells{j}(k).channels = points{j}(k).channels(includeall,:);
                pointsincells{j}(k).intensity = points{j}(k).intensity(includeall,:);
                pointsincells{j}(k).scaledIntensity = points{j}(k).scaledIntensity(includeall,:);

            end
        end
    end

end