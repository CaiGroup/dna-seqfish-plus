function [points, intensity, threshold, medianError] = level2samethreshold(image, refInt, startThreshold, typedots, varargin)%, superres)
% function to adjust the threshold that has the closest median pixel
% intensity of retrieved points
%
% Date: 11/15/2019

%[pointsTemp, ~, ~, ~] = detectdotsv2(I{f, ch}, threshold(f,ch)*0.8, 'log', false, '', 1);
%[points{f,ch}, intensity{f,ch}] = SuperResPoints(pointsTemp,I{f,ch},xyPixSize,zPixSize);
%medInt = median(intensity{f,ch});

% Test case: E14-rep3-1-DNAFISH
% 3 images - speed sequential: 1min parfor: 1.38min
% 9 images - speed sequential: 2.26min; parfor(3): 3.2min; parfor(9): 2min



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

    %% initial median 
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
    %xyPixSize = 1;
    %zPixSize = 1;
    %[points, ~] = SuperResPoints(pointsTemp,image,xyPixSize,zPixSize);
    medInt = median(intensity);
    medianError = refInt - medInt;
    foundRange = 1;
    if abs(medianError) > 100
        if medianError < 0
            high = startThreshold;
            newThreshold = high;
            lowRange = false;
        elseif medianError > 0
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
        medianError = refInt - medInt;
        %check if difference is low enough
        if abs(medianError) > 100

            if lowRange % started with low range
                if medianError > 0 % change to new low range
                    low = newThreshold;
                    newThreshold = newThreshold *1.3; % increase by 10%
                else
                    high = newThreshold;
                    break;
                end
            else
                % started with high range
                if medianError < 0
                    high = newThreshold;
                    newThreshold = newThreshold *0.7; % decrease by 10%
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
            medInt = median(intensity);
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
        [threshold, medianError, points, intensity] = binarysearchthreshold(low, high, refInt, image, typedots);
    end


end

function [threshold, medianError, points, intensity] = binarysearchthreshold(low, high, refInt, image, typedots)
    mid = (low + high) / 2;
    threshold = mid;
    [points, intensity, ~, ~] = detectdotsv2(image, threshold, typedots, false, '', 1);
    % cases for different types of super resolution
    %[points, intensity] = SuperResPoints(pointsTemp,image,1,1);
    medInt = median(intensity);
    medianError = refInt - medInt;
    if abs(medianError) < 100
        return;
    elseif medianError < 0
        high = mid;
        binarysearchthreshold(low, high, refInt, image, typedots);
    elseif medianError > 0
        low = mid;
        binarysearchthreshold(low, high, refInt, image, typedots);
    end
end

