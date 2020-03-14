function [] = outputdecodedpoints(experimentDir, experimentName, position, numRounds, ...
    numChannels, sqrtradius, segment, alloweddiff, typedots, superres, chArray, ...
    minseeds, chaTform, folderArray)
% function to decode points for each channel

    %% Decode from points
    % Set up Variables
    %{
    numRounds = 5;
    numChannels = 16;
    sqrtradius = 6;
    segment = 'whole';
    experimentName = '2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH'; 
    minseeds = 3;
    alloweddiff = 2;
    chArray = 1:2;
    typedots = 'log';
    superres = 'radial3d';
    %}


    pointsDir = fullfile(experimentDir, 'points');
    [pointsch, offsets] = alignpointswrapper(chArray, position, pointsDir, chaTform, numRounds, folderArray, experimentName); % need to adjust this to fix the pipeline

    % save data/ change the saved data
    save([experimentDir '\points\pointsAlignedBeadsJonathan-pos' num2str(position) '-chs1_2.mat'], 'pointsch', 'offsets', 'experimentDir', 'position', 'numRounds', ...
        'segment', 'sqrtradius', 'minseeds', 'alloweddiff', 'experimentName', 'chArray', 'chaTform');




    %% Decode the points


    % To load the processed images
    %dataDir = fullfile(experimentDir, 'processedimages\pos0\preProcessedData-pos0-2019-09-09-brain-rep2-2-DNAFISH-2019-09-28.mat');
    %load(dataDir, 'I');


    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, pointsch{channel}, segment, sqrtradius, typedots, superres, alloweddiff, ...
        channel, minseeds);
    end

end