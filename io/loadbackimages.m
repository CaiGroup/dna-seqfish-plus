function [backIms, backTform] = loadbackimages(backImPath, position, refImage)

    [backIms, numDapi, numZSlice, ~, ~] = grabimseries(backImPath, position);
    % dapi is the 4 in the cell
    numCh = numDapi - 1;

    if numZSlice < 16
        % divide image into 4 pieces from 1 image if zslices < 16
        [backImDivide, ~] = imdivideby4(backIms{numDapi});
        [backImRefDivide, numZSliceDivide] = imdivideby4(refImage);
    else
        backImDivide = backIms{numDapi};
        backImRefDivide = refImage;
        numZSliceDivide = 15; % arbitrary
    end

    initialRadius = 0.0625; %0.0625 for 3d is default
    numIterations = 100; % 100 is default
    backTform = grabtform(backImDivide, backImRefDivide, initialRadius, numIterations);
    % remove the z transformation
    if numZSliceDivide >= 16
        backTform.T(4,3) = 0;
    end
    % apply the tform
    for ch = 1:numCh
        if length(backTform.T) == 3
            outputView = imref2d(size(backIms{ch}));
        elseif length(backTform.T) == 4
            outputView = imref3d(size(backIms{ch}));
        end
        backIms{ch} = imwarp(backIms{ch}, backTform, 'OutputView', outputView);
    end

end