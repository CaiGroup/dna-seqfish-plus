function [points, pointsch, pointschpercell, intensity] ...
    = seqfishprocess_preprocessedimages(experimentDir, experimentName, experimentLabel, ...
    position, folderArray, superres, chaTform, usechabboffsets, physicalTforms, usephysicaloffsets, filtersigma, I,...
    segPath)


    %% Paths
    addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    

    %% Variables
    typedots = 'log';%'exons';% 'log'; % new log filter
    numCh = 2;
    refposition = 0;
    numRounds = 1; % non-barcoding analysis.

    %% get the reference threshold
    chArray = 1:numCh;
    % load Pos0 most recent thresholds
    %csvPath = fullfile(experimentDir, 'threshold', '002-016-E14-rep2-2-DNAFISH-ch1-2-hyb1-80-threshold-logfilter-Pos0vs3-20191208.csv');
    csvPath = getfile(fullfile(experimentDir, 'threshold'), 'threshold', 'match');
    fileExt = csvPath(end-2:end);
    if strcmp(fileExt, 'mat')
        % just load the data
        t = load(csvPath);
        threshold = t.threshold;
    elseif strcmp(fileExt, 'csv')
        threshold = threshcsv2matyodai(csvPath, chArray);
    end
    %% Output points form the images "I"
    [points, intensity, adjustedThreshold] = outputprocessimages_adjustedthresh(experimentDir, ...
    position, folderArray, numCh, typedots, superres, [], [], threshold, refposition, ...
    experimentLabel, filtersigma, I);

    % save the data
    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'variables');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    savePath = fullfile(saveDir, ['data-pos' num2str(position) '-' experimentName '-rawpoints-' experimentLabel '.mat']);
    save(savePath, 'points', 'intensity',  'adjustedThreshold');
    
    pointsch = cell(size(folderArray,2), 1);
    for i = 1:size(folderArray,2)
        pointsch{i} = struct('channels', cell(numCh,1), 'intensity', cell(numCh, 1));
    end
    
    for ch = chArray 
        for folder = folderArray
            idx = folder + 1;
            pointsTemp = points(:,ch);
            pointsTemp = pointsTemp{idx};
            
            if usechabboffsets
                pointsch{idx}(ch).channels = transformPointsForward(chaTform,  pointsTemp);
            else
                pointsch{idx}(ch).channels = pointsTemp;
            end
            if usephysicaloffsets
                pointsch{idx}(ch).channels = transformPointsForward(physicalTforms,  pointsch{idx}(ch).channels);
            end
            pointsch{idx}(ch).intensity = intensity{idx,ch};

        end
    end
    
    if usechabboffsets && usephysicaloffsets
        savePath = fullfile(saveDir, ['data-pos' num2str(position) '-points-chphycorr-' experimentLabel '.mat']);
        save(savePath, 'points', 'pointsch', 'intensity',  'adjustedThreshold');
    elseif usechabboffsets
        savePath = fullfile(saveDir, ['data-pos' num2str(position)  '-points-chcorr-' experimentLabel '.mat']);
        save(savePath, 'points', 'pointsch', 'intensity',  'adjustedThreshold');
    elseif usephysicaloffsets
        savePath = fullfile(saveDir, ['data-pos' num2str(position) '-points-phycorr-' experimentLabel '.mat']);
        save(savePath, 'points', 'pointsch', 'intensity',  'adjustedThreshold');  
    end
    
    %% find dots within ROIs. assuming 2D ROIs for now.
    if ~isempty(segPath)
        % calculate points within ROIs.
        [pointschpercell, numCells] = segment2dcells(segPath, chArray, position, pointsch);
        
        % Output the points per cell
        saveDir2 = fullfile(experimentDir, 'analysis',  experimentLabel, 'pointsList');
        if exist(saveDir2, 'dir') ~= 7
            mkdir(saveDir2);
        end
        savePointsFilePath = fullfile(saveDir2, ['pointsList-pos' num2str(position) '-' experimentLabel '.csv']); % Main data to save
        outputpointspercell(pointschpercell, savePointsFilePath, position)

        % save points and data
        save(savePath, 'points', 'pointsch', 'intensity', 'adjustedThreshold', 'numCells', 'pointschpercell', 'chArray', 'folderArray', 'position');
    end
    
end

function [pointspercell, numCells] = segment2dcells(segPath, channelArray, position, points)
    roiSegPath = fullfile(segPath, ['Pos' num2str(position)], 'RoiSet.zip');
    vertex = selfsegzip(roiSegPath);
    numCells = numel(vertex);
    numChannels = length(channelArray);
    pointspercell = cell(numCells, 1);
    for c = 1:numCells
        pointspercell{c} = points;
        numHybs = size(points, 1);
        for f = 1:numHybs
            for ch = 1:numChannels
                pointsSelection = points{f}(ch).channels;
                ind = inpolygon(pointsSelection(:,1), pointsSelection(:,2), vertex(c).x, vertex(c).y);
                pointspercell{c}{f}(ch).channels = pointspercell{c}{f}(ch).channels(ind,:);
                pointspercell{c}{f}(ch).intensity = pointspercell{c}{f}(ch).intensity(ind,:);
                %pointspercell{c}{f}(ch).scaledIntensity = pointspercell{c}{f}(ch).scaledIntensity(ind,:);
            end
        end
        
    end
end

function [] = outputpointspercell(pointspercell, savePointsFilePath, position)
    % variables and directory
    numCells = size(pointspercell,1);
    numHybs = size(pointspercell{1},1);
    numChannels = size(pointspercell{1,1}{1,1},1);

    % set up csv file
    fileID = fopen(savePointsFilePath,'w');
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                    'fov', 'cellID', 'hybID', 'chID', 'x', 'y', 'z', 'int');
    %fprintf(fileID,'\n'); % first new line doesn't create a new line
    
    for c = 1:numCells
        for f = 1:numHybs
            for ch = 1:numChannels
                numPoints = size(pointspercell{c}{f}(ch).channels,1);
                if numPoints>0
                    for p = 1:numPoints
                        pointSelection = pointspercell{c}{f}(ch);
                        x = pointSelection.channels(p,1);
                        y = pointSelection.channels(p,2);
                        z = pointSelection.channels(p,3);
                        int = pointSelection.intensity(p);
                        fprintf(fileID,'%d,%d,%d,%d,%.3f,%.3f,%.3f,%.1f\n', ...
                                                position, c, f, ch, x, y, z, int);
                    end
                end
            end
        end
    end

end

