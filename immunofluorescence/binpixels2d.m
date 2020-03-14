function [sumCellsHybs] = binpixels2d(shape2d, fov, savePath, channel, cellidx, endingString, ...
    image, binSize, folderArray, pixelsPerBin)

    % set up csv file
    listSavePath = fullfile(savePath, ['outputData-ImmunoFluorescence-Ch' num2str(channel) '-Pos' num2str(fov) '-Cell' num2str(cellidx) endingString]);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID,'%s,%s,%s,%s,%s\n', ...
        'hybID','avgint', 'xIdx', 'yIdx', 'zIdx');
    % fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s\n', ...
        %'fov', 'chID', 'cellID', 'hybID', 'avgint', 'xIdx', 'yIdx', 'zIdx');


    % get indices of the polygon shape
    [xlim,ylim] = boundingbox(polyshape(shape2d(cellidx).x,shape2d(cellidx).y));
    indZ = 1:size(image{1},3);
    % Get the start indices for the bounding box: take into account if odd
    x1Dim = xlim(1); % get top number
    x2Dim = xlim(2);
    y1Dim = ylim(1);
    y2Dim = ylim(2);
    indX = x1Dim:x2Dim;
    indY = y1Dim:y2Dim;
    if any(indX < 1)
        indX(indX < 1) = [];
    elseif any(indX > 2048)
        indX(indX > 2048) = [];
    end
    if any(indY < 1)
        indY(indY < 1) = [];
    elseif any(indY > 2048)
        indY(indY > 2048) = [];
    end
    hybIter = 1;
    sumCellsHybs = cell(1, length(folderArray));
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
                    if any(inpolygon(indShapeX, indShapeY, shape2d(cellidx).x, shape2d(cellidx).y))  && avgInt ~= 0
                        avgX = mean(indShapeX);
                        avgY = mean(indShapeY);

                        % print each line
                        %origIndices = [indShapeX, indShapeY, indShapeZ];
                        %origIndicesString = mat2str(origIndices);
                        fprintf(fileID,'%.0f,%.3f,%.2f,%.2f,%d\n', ...
                                        hybcycle, avgInt, avgX, avgY, z);
                        if ~isempty(sumCellsHybs{1, hybIter})
                            sumCellsHybs{hybIter} = sumCellsHybs{hybIter} + avgInt * pixelsPerBin;
                        else
                            sumCellsHybs{hybIter} = avgInt * pixelsPerBin;
                        end
                                    %fprintf(fileID,'%d,%d,%d,%d,%.3f,%.2f,%.2f,%d\n', ...
                                    %    fov, channel, c, hybcycle, avgInt, avgX, avgY, z);

                    end
                end
            end  
        end
        hybIter = hybIter + 1;
    end
    fclose(fileID);
end