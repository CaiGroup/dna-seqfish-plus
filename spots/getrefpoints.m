function [dots, intensity] = getrefpoints(images, threshold, typeDetectDots, savepoints, savepointsDir)
% get the points in each hyb as a reference

    % variables
    numCh = size(images, 2);
    numHybs= size(images, 1);
    dots = cell(numHybs, numCh);
    intensity = dots;
    % set up directories
    if exist(savepointsDir, 'dir') ~= 7
        mkdir(savepointsDir);
    end
    for ch = 1:numCh
        % set up directories
        savepointsChDir = fullfile(savepointsDir, ['ch' num2str(ch)]);
        if exist(savepointsChDir, 'dir') ~= 7
            mkdir(savepointsChDir);
        end
        for i = 1:numHybs
            [dots{i, ch}, intensity{i, ch}, ~, ~] = detectdotsv2(images{i, ch}, threshold{i}(ch), typeDetectDots, savepoints, i, savepointsChDir);
        end
    end


end