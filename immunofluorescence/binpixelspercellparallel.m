function [sumCellsHybs] = binpixelspercellparallel(binSize, stats, shape2d, image, savePath, endingString, fov, folderArray, channel, segoption)
% use 2d boundaries per z-slice
%
% use the bonding box to bin the image by 2x2 boxes
% for each z slice; for each bin; check if bin is in alphsape; calculate
% average intensity and get centroid
%
% segoption for '3d' or '2d'
% what is the option for '2d' same as vertex(each cell).x or y
% 
% Date: 8/26/2019
% Author: Nico Pierson

    %{
    % set up csv file
    listSavePath = fullfile(savePath, ['outputData-ImmunoFluorescence-Hyb' num2str(hybcycle) '-Pos' num2str(fov) '-Ch' num2str(channel) endingString '.csv']);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s\n', ...
                    'fov', 'hybID', 'cellID', 'avgint', 'xIdx', 'yIdx', 'zIdx');
    fprintf(fileID,'\n'); % first new line doesn't create a new line
%}


    pixelsPerBin = binSize * binSize;


    switch segoption
        case '3d'
            numCells = size(stats, 1);
            sumCellsHybs = cell(numCells, 1); % save sum of each cell
            parfor cellidx = 1:numCells
                [sumCellsHybs{cellidx}] = binpixels3d(shape2d, stats, fov, savePath, channel, cellidx, endingString, ...
                    image, binSize, folderArray, pixelsPerBin);
            end
    case '2d'
        numCells = numel(shape2d);
        sumCellsHybs = cell(numCells,1); % save sum of each cell
        parfor cellidx = 1:numCells
            [sumCellsHybs{cellidx}] = binpixels2d(shape2d, fov, savePath, channel, cellidx, endingString, ...
                image, binSize, folderArray, pixelsPerBin);
        end
    end
    
end
