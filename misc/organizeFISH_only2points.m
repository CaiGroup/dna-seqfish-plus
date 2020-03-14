function [points] = organizeFISH_only2points(FISH_only)
% organizes points FISH_only from Linus' seqFISH+ experiment into format of
% points

    numRounds = length(FISH_only);
    numChannels = max(FISH_only{1}(:,3));
    points = cell(1, numRounds);
    for i = 1:numRounds
        points{i} = struct('channels', cell(1, numChannels), 'intensity', cell(1, numChannels), 'scaledIntensity', cell(1, numChannels));
        for j = 1:numChannels
            tempPoints = FISH_only{i};
            points{i}(j).channels = tempPoints(tempPoints(:,3) == j,:);
            % set zslice equal to 1
            points{i}(j).channels(:,3) = 1; % z is 1 for 2d images
            points{i}(j).intensity = ones(size(points{i}(j).channels, 1), 1);
            points{i}(j).scaledIntensity = ones(size(points{i}(j).channels, 1), 1);
        end
    end

end

%save('data-Linus_cleared_091218_brain-FP-Test-sqrt3', 'alloweddiff', 'barcodekey', 'numChannels', 'numRounds', 'points', 'position', 'projectName', 'sqrtradius', 'projectDir');