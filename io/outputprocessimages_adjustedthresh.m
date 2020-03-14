function [points, intensity, adjustedThreshold] = outputprocessimages_adjustedthresh(experimentDir, ...
    position, folderArray, numCh, typedots, superres, medRefIntensity, medIntensity, threshold, refposition, ...
    experimentLabel, filtersigma, I)
% output the processeed images depending on the position so the points can
% be extracted for the bead alignment

    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end

    saveFigDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', ['points_debug_pos' num2str(position)]);
    if exist(saveFigDir, 'dir') ~= 7
        mkdir(saveFigDir);
    end
    
    
    
    %% Declare Variables
    points = cell(length(folderArray),numCh);
    intensity = cell(length(folderArray), numCh);
    intcheck = cell(length(folderArray), numCh);
    sigma = cell(length(folderArray), numCh);
    xyPixSize = 1;
    zPixSize = 1;
    adjustedThreshold = ones(length(folderArray),numCh) * 999999;
            


    %% load pos1-4 images, grab points using same treshold vales
    % load threshold
    %load('I:\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\threshold\thresholdAllCh-pos0-numHybCycles85-numCh-3-E14-DNA-seqFISH+rep3-1-DNAFISH.mat');
    % get images and grab points

 
    for f = 1:length(folderArray)
        for ch = 1:numCh
            if ~isempty(medIntensity)&&~isempty(medRefIntensity)
                adjustedThreshold(f,ch) = threshold(f,ch) * (double(medIntensity{f,ch})/double(medRefIntensity{f,ch}));
            else
                adjustedThreshold(f,ch) = threshold(f,ch);
            end
            %% grab the points
            saveFigPath = fullfile(saveFigDir, ['figure-Pos' num2str(position) '-Folder' num2str(f) '-Ch' num2str(ch) '.fig']);
            [pointsTemp, intcheck{f,ch}, ~, ~] = detectdotsv2(I{f, ch}, adjustedThreshold(f,ch), typedots, true, saveFigPath, 1);
            switch superres
                case 'radial3d'
                    [points{f,ch}, intensity{f,ch}, sigma{f,ch}] = SuperResPoints(pointsTemp,I{f,ch},xyPixSize,zPixSize);
                case 'gaussian'
                    [points{f,ch}, intensity{f,ch}] = getgaussian(pointsTemp, I{f,ch});
                    [~, ~, sigma{f,ch}] = SuperResPoints(pointsTemp,I{f,ch},xyPixSize,zPixSize, false);
            end

            % filter using sigma
            if filtersigma && ~isempty(points{f,ch})
                removeidx = [];
                for i = 1:size(points{f,ch},1)
                    sigcheck = sigma{f,ch}(i);
                    if sigcheck < 0.58 % just to remove shot noise.
                        removeidx = cat(1, removeidx, i);
                    end
                end
                points{f,ch}(removeidx, :) = [];
                intensity{f,ch}(removeidx, :) = [];
                sigma{f,ch}(removeidx, :) = [];
            end

            %}
            I{f,ch} = [];

        end
    end

    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'points');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end
    if position ~= refposition
        savePath = fullfile(saveDirPath, ['points-int-thresh-pos' num2str(position) '-' experimentLabel '.mat']);
    else
        savePath = fullfile(saveDirPath, ['points-int-thresh-pos' num2str(position) '-' experimentLabel '.mat']);
    end
    save(savePath, 'points', 'intensity', 'adjustedThreshold', 'sigma');
    %toc

end