function threshold = orgthrehsold2round(thresh, numRounds, numChannels)
% organize a hyb by 1 threshold to round by channel

    threshold = ones(numRounds, numChannels) * 999999;
    for i = 1:numRounds
        for j = 1:numChannels
            idx = (i-1)*numChannels+j;
            threshold(i,j) = thresh(idx);
        end
    end

end