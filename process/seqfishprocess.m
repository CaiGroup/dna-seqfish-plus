function [rawpoints, intensity, pointsCsvName, numRefPoints, medRefIntensity] ...
    = seqfishprocess(experimentDir, experimentName, experimentLabel, ...
    position, folderArray, thresholdadjust, superres, numChannels, ...
    numRounds, sizeZ, chaTform, numRefPoints, medRefIntensity, filtersigma)


    %% Paths
    addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    

    %% Variables
    refposition = 0;
    typedots = 'log';%'exons';% 'log'; % new log filter
    numCh = 2;



    %% get the reference threshold
    chArray = 1:numCh;
    % load Pos0 most recent thresholds
    %csvPath = fullfile(experimentDir, 'threshold', '002-016-E14-rep2-2-DNAFISH-ch1-2-hyb1-80-threshold-logfilter-Pos0vs3-20191208.csv');
    csvPath = getfile(fullfile(experimentDir, 'threshold'), 'threshold', 'match');
    threshold = threshcsv2matyodai(csvPath, chArray);


    
    %% Roi Mask
    roiPath = fullfile(experimentDir, 'segmentation', ['Pos' num2str(refposition)], 'RoiSet.zip');
    if exist(roiPath, 'dir') ~= 7
        vertex = selfsegzip(roiPath);
        roimask = roi2mask(vertex, [2048 2048 sizeZ]);
    end
    numPointChannels = size(threshold,2);


    %% Output points form the images
    [points, intensity, shadingcorr, adjustedThreshold] = outputprocessimages(experimentDir, ...
        position, folderArray, numCh, typedots, superres, medRefIntensity, numRefPoints, threshold, ...
        refposition, experimentLabel, roimask, thresholdadjust, filtersigma);

    %% Get the number of points and the median intensity to save using the roi masks as a filter
    if isempty(numRefPoints) && isempty(medRefIntensity)
        [numPoints, medIntensity] = getnumrefpoints(points, intensity, chArray, numRounds, numChannels, roimask);
        if position == refposition
            medRefIntensity = medIntensity;
            numRefPoints = numPoints;
        end
    end

    % save the data
    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'variables');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    savePath = fullfile(saveDir, ['data-pos' num2str(position) '-' experimentName '-' experimentLabel '.mat']);
    save(savePath, 'points', 'intensity', 'shadingcorr', 'adjustedThreshold', ...
        'numRefPoints', 'medRefIntensity', 'roimask');


    %% output the points and get them, also saves variables
    [rawpoints, intensity, pointsCsvName] = outputrawpoints(points, intensity, adjustedThreshold, ...
        folderArray, chaTform, experimentDir, experimentName, position, ...
        numPointChannels, shadingcorr, experimentLabel);

end

