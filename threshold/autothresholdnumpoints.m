function thresholdHold = autothresholdnumpoints(I, numTargetPoints, numIters, numWorkers, channelArray, typedots)
% gets the threshold based on the number of points
%
% Date: 8/14/2019

    %% Initialize Variables
    numHybCycles = size(I, 1);
    thresholdHold = cell(numHybCycles, 1);
    numChannels = length(channelArray);
    threshold = cell(numHybCycles, numChannels);
    switch typedots
        case 'exons'
            x1 = 500;
            x2 = 100000;
        case 'introns'
            x1 = 1; % start at 100
            x2 = 7500;
    end
    
    %% Find the Minimum Threshold
    if numWorkers > 1 % how to check for number of workers on a computer
        delete(gcp('nocreate'));
        parpool('local', numWorkers);
        parfor h = 1:numHybCycles
            if isempty(thresholdHold{h})
                minThreshold = cell(1, numChannels);
                for ch = 1:numChannels
                    minThreshold{ch} = findminthreshold(x1, x2, numIters, numTargetPoints, I{h,ch}, typedots);
                end
                thresholdHold{h} = cell2mat(minThreshold);
            end
        end
    else
        for h = 1:numHybCycles
            if isempty(thresholdHold{h})
                minThreshold = cell(1, numChannels);
                for ch = 1:numChannels
                    minThreshold{ch} = findminthreshold(x1, x2, numIters, numTargetPoints, I{h,ch}, typedots);
                end
                thresholdHold{h} = cell2mat(minThreshold);
            end
        end
    end



    %% Reorganize threshold
    %{
    for i = 1:numHybCycles
        for ch = 1:numChannels
            threshold{i,ch} = thresholdHold{i}{ch};
        end
    end
    %}
    
end