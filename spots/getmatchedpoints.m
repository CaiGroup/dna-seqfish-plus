function points = getmatchedpoints(allPoints)
% function gets all similar points across all channels using reference points
%
% output: points structure with field 'channel' for every matching point
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 4/1/2019

    numChannels = length(allPoints);
    points = struct('channels', cell(1,numChannels), 'intensity', cell(1,numChannels));
    
    for k = 1:length(allPoints{1}.channels)
        index = zeros(1, numChannels);
        index(1) = k; % initialize index to row in first channel
        for channel = 2:numChannels
            [~,~,indextemp] = intersect(allPoints{1}.channels(k,:),allPoints{channel}.ref(:,:),'rows');
            if ~isempty(indextemp)
                index(channel) = indextemp;
            end
        end
        if any(~any(index, 1)) % test if any indices are zero
            % do nothing
        else
            for channel = 1:numChannels
                if channel ~= 1
                    points(channel).channels = cat(1, points(channel).channels, allPoints{channel}.points(index(channel),:));
                    points(channel).intensity = cat(1, points(channel).intensity, allPoints{channel}.intmatch(index(channel),:));
                else
                    points(channel).channels = cat(1, points(channel).channels, allPoints{channel}.channels(k,:));
                    points(channel).intensity = cat(1, points(channel).intensity, allPoints{channel}.intensity(k,:));
                end
            end
        end

    end
end

