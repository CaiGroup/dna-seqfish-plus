function [sumCellsHybs] = binpixelspercell(binSize, stats, shape2d, image, savePath, endingString, fov, folderArray, channel, segoption)
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
            numHybCycles = length(folderArray);
            sumCellsHybs = cell(numCells, numHybCycles); % save sum of each cell
            for c = 1:numCells
                
                % set up csv file
                listSavePath = fullfile(savePath, ['outputData-ImmunoFluorescence-Ch' num2str(channel) '-Pos' num2str(fov) '-Cell' num2str(c) endingString]);
                fileID = fopen(listSavePath,'w');
                fprintf(fileID,'%s,%s,%s,%s,%s\n', ...
                    'hybID','avgint', 'xIdx', 'yIdx', 'zIdx');

               % fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                    %'fov', 'chID', 'cellID', 'hybID', 'avgint', 'xIdx', 'yIdx', 'zIdx');
                
                boundbox = stats.BoundingBox(c,:);
                % get the z
                z1Dim = ceil(boundbox(3)); % get lower but change
                z2Dim = z1Dim + boundbox(6) - 1;
                if z1Dim == 0 % change zslice to 1
                    z1Dim = 1;
                end
                indZ = z1Dim:1:z2Dim; % increment one for each zslice
                
                % Get the start indices for the bounding box: take into account if odd
                x1Dim = ceil(boundbox(1)); % get top number
                x2Dim = x1Dim + boundbox(4) - 1;
                y1Dim = ceil(boundbox(2));
                y2Dim = y1Dim + boundbox(5) - 1;
                indX = x1Dim:x2Dim;
                indY = y1Dim:y2Dim;
                hybIter = 1;
                for hybcycle = folderArray % need to subtract by 1 for original folder
                    rawimage = image{hybIter}; % assign to each folder
                    for z = indZ
                    
                        if shape2d{z}{c}.Alpha ~= Inf
                            
                            boundI = rawimage(indY, indX, z);
                            binI = binavg(boundI, binSize);  % binning image by 2 in x and y % check if pixel is in bounding box....x1 correlates to 1-2: 2 to 3-4; 3 to 5-6
                            [sizeY, sizeX, sizeZ] = size(binI); % how to use 2 x 2 matrix and check if in bounding box
                            % set up base indices
                            baseInd = [binSize, binSize];
                            IND = 1:pixelsPerBin;
                            [baseY,baseX] = ind2sub(baseInd,IND);

                            for i = 1:sizeX
                                for j = 1:sizeY % left most point for each bin
                                    % original indices to check if in shape
                                    xInc = (i - 1) * binSize;
                                    yInc = (j - 1) * binSize;
                                    indShapeX = (baseX + xInc + x1Dim - 1)'; % add from original x1Dim and y1Dim
                                    indShapeY = (baseY + yInc + y1Dim - 1)';
                                    avgInt = binI(j,i);
                                    %indShapeZ = (ones(1,length(indShapeX))*z)';


                                    % can you check if in shape for array
                                    
                                    %[c hybcycle z i j]
                                    %if c == 1 && hybcycle == 1 && z == 2 && i == 16 && j == 20
                                    %    disp('debug');
                                    %end
                                    if any(inShape(shape2d{z}{c}, indShapeX, indShapeY)) && avgInt ~= 0

                                        avgX = mean(indShapeX);
                                        avgY = mean(indShapeY);

                                        % print each line
                                        %origIndices = [indShapeX, indShapeY, indShapeZ];
                                        %origIndicesString = mat2str(origIndices);

                                        fprintf(fileID,'%.0f,%.3f,%.2f,%.2f,%d\n', ...
                                            hybcycle, avgInt, avgX, avgY, z);
                                        if ~isempty(sumCellsHybs{c, hybIter})
                                            sumCellsHybs{c, hybIter} = sumCellsHybs{c} + avgInt * pixelsPerBin;
                                        else
                                            sumCellsHybs{c, hybIter} = avgInt * pixelsPerBin;
                                        end
                                        %fprintf(fileID,'%d,%d,%d,%d,%.3f,%.2f,%.2f,%d\n', ...
                                        %    fov, channel, c, hybcycle, avgInt, avgX, avgY, z);
                                    end

                                end 
                            end
                        end
                    end
                    hybIter = hybIter + 1;
                end
                fclose(fileID);
            end
    case '2d'
        numCells = numel(shape2d);
        numHybCycles = length(folderArray);
        sumCellsHybs = cell(numCells, numHybCycles); % save sum of each cell
        for c = 1:numCells
            
            % set up csv file
            listSavePath = fullfile(savePath, ['outputData-ImmunoFluorescence-Ch' num2str(channel) '-Pos' num2str(fov) '-Cell' num2str(c) endingString]);
            fileID = fopen(listSavePath,'w');
            fprintf(fileID,'%s,%s,%s,%s,%s\n', ...
                'hybID','avgint', 'xIdx', 'yIdx', 'zIdx');
           % fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                %'fov', 'chID', 'cellID', 'hybID', 'avgint', 'xIdx', 'yIdx', 'zIdx');
            
            
            % get indices of the polygon shape
            [xlim,ylim] = boundingbox(polyshape(shape2d(c).x,shape2d(c).y));
            indZ = 1:size(image{1},3);
            % Get the start indices for the bounding box: take into account if odd
            x1Dim = xlim(1); % get top number
            x2Dim = xlim(2);
            y1Dim = ylim(1);
            y2Dim = ylim(2);
            indX = x1Dim:x2Dim;
            indY = y1Dim:y2Dim;
            hybIter = 1;
            for hybcycle = folderArray % need to subtract by 1 for original folder
                rawimage = image{hybcycle}; % assign to each folder
                for z = indZ

                
                %for hybcycle = folderArray % need to subtract by 1 for original folder
                    
                    boundI = rawimage(indY, indX, z);
                    binI = binavg(boundI, binSize);  % binning image by 2 in x and y % check if pixel is in bounding box....x1 correlates to 1-2: 2 to 3-4; 3 to 5-6
                    [sizeY, sizeX, sizeZ] = size(binI); % how to use 2 x 2 matrix and check if in bounding box
                    % set up base indices
                    baseInd = [binSize, binSize];
                    IND = 1:pixelsPerBin;
                    [baseY,baseX] = ind2sub(baseInd,IND);

                    for i = 1:sizeX
                        for j = 1:sizeY % left most point for each bin
                            % original indices to check if in shape
                            xInc = (i - 1) * binSize;
                            yInc = (j - 1) * binSize;
                            indShapeX = (baseX + xInc + x1Dim - 1)'; % add from original x1Dim and y1Dim
                            indShapeY = (baseY + yInc + y1Dim - 1)';
                            avgInt = binI(j,i);
                            %indShapeZ = (ones(1,length(indShapeX))*z)';

                            % can you check if in shape for array
                            if any(inpolygon(indShapeX, indShapeY, shape2d(c).x, shape2d(c).y))  && avgInt ~= 0
                                avgX = mean(indShapeX);
                                avgY = mean(indShapeY);

                                % print each line
                                %origIndices = [indShapeX, indShapeY, indShapeZ];
                                %origIndicesString = mat2str(origIndices);
                                fprintf(fileID,'%.0f,%.3f,%.2f,%.2f,%d\n', ...
                                                hybcycle, avgInt, avgX, avgY, z);
                                if ~isempty(sumCellsHybs{c, hybIter})
                                    sumCellsHybs{c, hybIter} = sumCellsHybs{c} + avgInt * pixelsPerBin;
                                else
                                    sumCellsHybs{c, hybIter} = avgInt * pixelsPerBin;
                                end
                                            %fprintf(fileID,'%d,%d,%d,%d,%.3f,%.2f,%.2f,%d\n', ...
                                            %    fov, channel, c, hybcycle, avgInt, avgX, avgY, z);

                            end
                        end
                    end  
                end
                hybIter = hybIter + 1;
            end
        end
        fclose(fileID);
    end
    %fclose(fileID);
end
