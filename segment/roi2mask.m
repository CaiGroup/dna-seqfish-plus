function BW = roi2mask(vertex, imageSize)
% Makes a binary mask from vertices of 2d ROIs
% 
% Author: Nico Pierson
% Date: October 11, 2018
% imageSize is usually [2048 2048 1] for a single .tif image

    % Set size cells
    numCells = length(vertex);
    BW = zeros(imageSize);
    x = imageSize(1);
    y = imageSize(2);
    z = imageSize(3);

    for i = 1:numCells

        pm = poly2mask(vertex(i).x,vertex(i).y,x,y);
        pm3 = repmat(pm,[1 1 z]);

        % Add polymask to binaryMask
        BW = or(BW, pm3); 
    end

end