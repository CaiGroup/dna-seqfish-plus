function [backIms, numDapi, maxZ, backTform] = alignbackimages(dapiRef, backImPath, position)
% function aligns the background images

    [backIms, numDapi, maxZ, ~, ~] = grabimseries(backImPath, position);
    numCh = numDapi - 1;
    initialRadius = 0.0625; %0.0625 for 3d is default
    numIterations = 100; % 100 is default
    backTform = grabtform(backIms{numDapi}, dapiRef, initialRadius, numIterations);

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