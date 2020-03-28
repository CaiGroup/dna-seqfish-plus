function I = getchim(dimensionOrder, chArray, sizeC, sizeZ, numImages, result)
% gets the image from the specific channel from the cell array in result of
% grabimseries.m
%

    %% variables
    I = cell(1, length(chArray));
    
    
    %% get the image from result
    index = 1;
    for i = chArray
        if strcmp(dimensionOrder, 'XYCZT')
            indices = i:sizeC:numImages;
        elseif strcmp(dimensionOrder, 'XYZCT')
            startIndex = sizeZ * (i - 1) + 1;
            endIndex = startIndex + sizeZ - 1;
            indices = startIndex:endIndex;
        end
        for j = indices
            I{index} = cat(3,I{index}, result{1,1}{j,1});
        end
        index = index + 1;
    end

end