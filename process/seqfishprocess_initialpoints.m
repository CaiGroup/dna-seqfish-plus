function [rawpoints, intensity, pointsCsvName, I] ...
    = seqfishprocess_initialpoints(experimentDir, experimentName, experimentLabel, ...
    position, folderArray, superres, chaTform, filtersigma, numCh)

    
    %% Paths
    addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    

    %% Variables
    refposition = 0;
    typedots = 'log';%'exons';% 'log'; % new log filter
    %numCh = 2;
    thresholdadjust = false;


    %% get the reference threshold
    chArray = 1:numCh;
    % load Pos0 most recent thresholds
    %csvPath = fullfile(experimentDir, 'threshold', '002-016-E14-rep2-2-DNAFISH-ch1-2-hyb1-80-threshold-logfilter-Pos0vs3-20191208.csv');
    csvPath = getfile(fullfile(experimentDir, 'threshold'), 'threshold', 'match');
    threshold = threshcsv2matyodai(csvPath, chArray);

    %% Output points form the images
    numPointChannels = size(threshold,2);
    [points, intensity, shadingcorr, adjustedThreshold, I] = outputprocessimages(experimentDir, ...
        position, folderArray, numCh, typedots, superres, [], [], threshold, ...
        refposition, experimentLabel, [], thresholdadjust, filtersigma);



    %% save the data
    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'variables');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    savePath = fullfile(saveDir, ['data-pos' num2str(position) '-' experimentName '-step2-' experimentLabel '.mat']);
    save(savePath, 'points', 'intensity', 'shadingcorr', 'adjustedThreshold');


    %% output the points and get them, also saves variables
    [rawpoints, intensity, pointsCsvName] = outputrawpoints(points, intensity, adjustedThreshold, ...
        folderArray, chaTform, experimentDir, experimentName, position, ...
        numPointChannels, shadingcorr, experimentLabel);

end

