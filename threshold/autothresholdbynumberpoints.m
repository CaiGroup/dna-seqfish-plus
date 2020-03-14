function [points, intensity, threshold, numpointError] = autothresholdbynumberpoints(image, refNumPoints, startThreshold, typedots, varargin)%, superres)
% get the same number of points for the image



    %% Set up optional Parameters
    argsLimit = 1;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:adjustthreshold:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end   
    % set defaults for optional inputs
    optargs = {[]};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [roimask] = optargs{:};



% find the range-how...1st try gives low or high...move 20% each time, if
% low set as new low if not move 20% again from the threshold

    %% initial number of points
    errorThreshold = 35;
    threshold = [];
    [points, intensity, ~, ~] = detectdotsv2(image, startThreshold, typedots, false, '', 1);
    removeInd = [];
    if ~isempty(roimask)
        numPoints = size(points,1);
        for i = 1:numPoints
            if ~roimask(points(i,2), points(i,1), points(i,3))
                removeInd = cat(1, removeInd, i);
            end
        end
        points(removeInd,:) = [];
        intensity(removeInd) = [];
    end
    % cases for different types of super resolution
    numNewPoints = size(points, 1);
    numpointError = refNumPoints - numNewPoints;
    foundRange = 1;
    if abs(numpointError) > errorThreshold
        if numpointError > 0
            high = startThreshold;
            newThreshold = high;
            lowRange = false;
        elseif numpointError < 0
            low = startThreshold;
            newThreshold = low;
            lowRange = true;
        end
    else
        threshold = startThreshold;
        foundRange = 0;
    end

    % Find low and high range
    while(foundRange)
        numpointError = refNumPoints - numNewPoints;
        %check if difference is low enough
        if abs(numpointError) > errorThreshold

            if lowRange % started with low range
                if numpointError < 0 % change to new low range
                    low = newThreshold;
                    newThreshold = newThreshold *1.2; % increase by 10%
                else
                    high = newThreshold;
                    break;
                end
            else
                % started with high range
                if numpointError > 0
                    high = newThreshold;
                    newThreshold = newThreshold *0.8; % decrease by 10%
                else
                    low = newThreshold;
                    break;
                end
            end

            [points, intensity, ~, ~] = detectdotsv2(image, newThreshold, typedots, false, '', 1);
            % check if points are in roi mask
            removeInd = [];
            if ~isempty(roimask)
                numPoints = size(points,1);
                for i = 1:numPoints
                    if ~roimask(points(i,2), points(i,1), points(i,3))
                        removeInd = cat(1, removeInd, i);
                    end
                end
                points(removeInd,:) = [];
                intensity(removeInd) = [];
            end
            
            % cases for different types of super resolution
            %[points, ~] = SuperResPoints(pointsTemp,image,xyPixSize,zPixSize);
            numNewPoints = size(points, 1);
        else
            threshold = newThreshold;
            break;
        end

    end


    % binary search for right threshold
    if isempty(threshold)
        if low == high
            disp('debug here');
        end
        [threshold, numpointError, points, intensity] = binarysearchthreshold(low, high, refNumPoints, image, typedots, roimask);
    end


end

function [threshold, pointError, points, intensity] = binarysearchthreshold(low, high, refNumPoints, image, typedots, roimask)
    mid = (low + high) / 2;
    threshold = mid;
    [points, intensity, ~, ~] = detectdotsv2(image, threshold, typedots, false, '', 1);
    
    % get points only in roimask
    errorThreshold = 35;
    removeInd = [];
    if ~isempty(roimask)
        numPoints = size(points,1);
        for i = 1:numPoints
            if ~roimask(points(i,2), points(i,1), points(i,3))
                removeInd = cat(1, removeInd, i);
            end
        end
        points(removeInd,:) = [];
        intensity(removeInd) = [];
    end
            
    numNewPoints = size(points, 1);
    pointError = refNumPoints - numNewPoints;
    if abs(pointError) < errorThreshold
        return;
    elseif pointError < 0
        high = mid;
        binarysearchthreshold(low, high, refNumPoints, image, typedots, roimask);
    elseif pointError > 0
        low = mid;
        binarysearchthreshold(low, high, refNumPoints, image, typedots, roimask);
    end
end

