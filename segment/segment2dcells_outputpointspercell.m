function [pointspercell, numCells] = segment2dcells_outputpointspercell(experimentDir, experimentLabel, chArray, position, points)

segPath = fullfile(experimentDir, 'segmentation');  
[pointspercell, numCells] = segment2dcells(segPath, chArray, position, points);
saveDir = fullfile(experimentDir, 'analysis',  experimentLabel, 'pointsList');
if exist(saveDir, 'dir') ~= 7
    mkdir(saveDir);
end
savePointsFilePath = fullfile(saveDir, ['pointsList-pos' num2str(position) '-' experimentLabel '.csv']); % Main data to save
outputpointspercell(pointspercell, savePointsFilePath, position);

end

function [pointspercell, numCells] = segment2dcells(segPath, channelArray, position, points)
    roiSegPath = fullfile(segPath, ['Pos' num2str(position)], 'RoiSet.zip');
    vertex = selfsegzip(roiSegPath);
    numCells = numel(vertex);
    numChannels = length(channelArray);
    pointspercell = cell(numCells, 1);
    for c = 1:numCells
        pointspercell{c} = points;
        numHybs = size(points, 1);
        for f = 1:numHybs
            for ch = 1:numChannels
                pointsSelection = points{f}(ch).channels;
                ind = inpolygon(pointsSelection(:,1), pointsSelection(:,2), vertex(c).x, vertex(c).y);
                pointspercell{c}{f}(ch).channels = pointspercell{c}{f}(ch).channels(ind,:);
                pointspercell{c}{f}(ch).intensity = pointspercell{c}{f}(ch).intensity(ind,:);
                %pointspercell{c}{f}(ch).scaledIntensity = pointspercell{c}{f}(ch).scaledIntensity(ind,:);
            end
        end
        
    end
end

function [] = outputpointspercell(pointspercell, savePointsFilePath, position)
    % variables and directory
    numCells = size(pointspercell,1);
    numHybs = size(pointspercell{1},1);
    numChannels = size(pointspercell{1,1}{1,1},1);

    % set up csv file
    fileID = fopen(savePointsFilePath,'w');
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                    'fov', 'cellID', 'hybID', 'chID', 'x', 'y', 'z', 'int');
    %fprintf(fileID,'\n'); % first new line doesn't create a new line
    
    for c = 1:numCells
        for f = 1:numHybs
            for ch = 1:numChannels
                numPoints = size(pointspercell{c}{f}(ch).channels,1);
                if numPoints>0
                    for p = 1:numPoints
                        pointSelection = pointspercell{c}{f}(ch);
                        x = pointSelection.channels(p,1);
                        y = pointSelection.channels(p,2);
                        z = pointSelection.channels(p,3);
                        int = pointSelection.intensity(p);
                        fprintf(fileID,'%d,%d,%d,%d,%.3f,%.3f,%.3f,%.1f\n', ...
                                                position, c, f, ch, x, y, z, int);
                    end
                end
            end
        end
    end

end
