function I = applydapitform(image, tform)    
% applys tform to image in cell array
%
% Return image with applied transformation
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 2/26/2019
% Modified:

    % Variables
    numChannels = length(image);
    numDim = length(tform.T);
    imageIsCell = false;
    if numDim ~= 3 && numDim ~= 4
        error 'tform does not have the right number of dimensions';
    end
    
    % Check if image is cell with difference channels
    if iscell(image)
        imageIsCell = true;
        I = cell(1, numChannels);
    end

    if numDim == 4
        if imageIsCell
            for i = 1:numChannels
                I{i} = imwarp(image{i}, tform, 'OutputView', imref3d(size(image{i})));
            end
        else
            I = imwarp(image, tform, 'OutputView', imref3d(size(image)));
        end
    elseif numDim == 3
        if imageIsCell
            numZ = size(image{1},3);
            for i = 1:numChannels
                for z = 1:numZ
                    I{i}(:,:,z) = imwarp(image{i}(:,:,z), tform, 'OutputView', imref2d(size(image{i}(:,:,z))));
                end
            end
        else
            numZ = size(image,3);
            for z = 1:numZ
                I(:,:,z) = imwarp(image(:,:,z), tform, 'OutputView', imref2d(size(image(:,:,z))));
            end
        end
    end
    
    
end