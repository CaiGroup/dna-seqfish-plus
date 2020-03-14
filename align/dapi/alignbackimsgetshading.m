function [shadingcorr, backIms, backTform] = alignbackimsgetshading(dapiRef, experimentDir, ...
    backgroundFolderName, position)
% aligns the background images to the first hyb dapi and gets the uneven
% illumination shading corrections
%
% Date: 8/16/2019


    % Get background images
    backImBasePath = fullfile(experimentDir, backgroundFolderName);
    backImPath = fullfile(backImBasePath, ['MMStack_Pos' num2str(position) '.ome.tif']);
    if exist(backImPath, 'file') == 2
        [backIms, numDapi, numZSlice, ~, ~] = grabimseries(backImPath, position);
        % dapi is the 4 in the cell
        numCh = numDapi - 1;

        if numZSlice < 16
            % divide image into 4 pieces from 1 image if zslices < 16
            [backImDivide, ~] = imdivideby4(backIms{numDapi});
            [backImRefDivide, numZSliceDivide] = imdivideby4(dapiRef);
        else
            backImDivide = backIms{numDapi};
            backImRefDivide = dapiRef;
            numZSliceDivide = 15; % arbitrary
        end

        initialRadius = 0.0625; %0.0625 for 3d is default
        numIterations = 100; % 100 is default
        lastwarn('');
        backTform = grabtform(backImDivide, backImRefDivide, initialRadius, numIterations);
        [warnMsg, warnId] = lastwarn;
        if isempty(warnMsg) % check if there is a warning for image registration
            for ch = 1:numCh
                if length(backTform.T) == 3
                    outputView = imref2d(size(backIms{ch}));
                elseif length(backTform.T) == 4
                    outputView = imref3d(size(backIms{ch}));
                end
                backIms{ch} = imwarp(backIms{ch}, backTform, 'OutputView', outputView);
            end
        else
            warning('skipped dapi alignment of background image to first dapi');
        end

        % Get Shading Corrections
        shadingcorr = shadingcorrection(backIms(1:numCh));
    else
        error('background image: %s\nNot Found', backImPath);
    end


end