function [threshold, medianIntensity] = getthreshmedianintfrompoints(experimentDir, ...
    position, numCh, numRounds, numChannels)
% get the threshold from each channel and median intensity from the points
% using E14rep2-2 older data

    %% Get Median Intensity
    %experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
    %position = 0;
    %numCh = 2;
    %numRounds = 5;
    %numChannels = 16;
    numHybs = numRounds * numChannels;
    medianIntensity = cell(numHybs, numCh);
    for ch = 1:numCh
        % get the points
        pointsPath = getfile(fullfile(experimentDir, ['analysis\exons-spot-detector\2error-sqrt6-ch' num2str(ch)]), ['points-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped-Pos' num2str(position) '-2019-11-06.mat'], 'strict');
        load(pointsPath, 'points');
        for r = 1:numRounds
            for c = 1:numChannels
                hybidx = (r-1)*numChannels + c;
                medianIntensity{hybidx, ch} = median(points{r}(c).intensity);
            end
        end
    end

    %% load the thresholds of pos0 for both channels
    threshold = ones(numHybs, numCh) * 999999;
    t1  = load('I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\threshold\ch1\rep2-2-channel1-thresholds.mat', 'threshold');
    threshold(:,1) = ch2wholethreshold(t1.threshold);
    t2 = load('I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\threshold\ch2\rep2-2-channel2-thresholds.mat', 'threshold');
    threshold(:,2) = ch2wholethreshold(t2.threshold);

end



