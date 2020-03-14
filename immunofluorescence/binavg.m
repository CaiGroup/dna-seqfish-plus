function binI = binavg(I, binSize)
% averages the image depending on the bin size, 2 = bin images by 2 or in
% half.
%
% Assume I is rectangular

    i2=1;
    j2=1;
    binStep = binSize - 1;
    [xDim, yDim, zDim] = size(I);
    binI = zeros(ceil(xDim/binSize), ceil(yDim/binSize), zDim);
    for i = 1:binSize:xDim
        for j = 1:binSize:yDim
            for k = 1:zDim
                % check if indies are out of bounds
                xlimit = getlastdim(i, xDim, binStep);
                ylimit = getlastdim(j, yDim, binStep);
                sumXInd = i:xlimit;
                sumYInd = j:ylimit;
                numBinPixels = length(sumXInd) * length(sumYInd);
                binI(i2, j2, k) = sum(sum(I(sumXInd, sumYInd, k))) / numBinPixels;
            end
            j2 = j2 + 1;
        end
        i2 = i2 + 1;
        j2 = 1;
    end
    

end

function limit = getlastdim(init, dimSize, binStep)
    stepInd = init + binStep;
    if stepInd <= dimSize
        limit = stepInd;
    else
        % find the limit
        for i = stepInd:-1:init
            startLimit = i - 1;
            if startLimit <= dimSize
                limit = startLimit;
                break;
            end
        end
    end
end