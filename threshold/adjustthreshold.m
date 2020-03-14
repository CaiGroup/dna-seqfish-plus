function [points, intensity, threshold, medianError] = adjustthreshold(image, ...
    refInt, startThreshold, typedots, varargin)%, superres)
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
    medianError = [];
    [points, intensity, ~, ~] = detectdotsv2(image, startThreshold, typedots, false, '', 1);
    % only get the median intensity of the points inside the rois
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
    medInt = median(intensity);
    threshold = round((medInt*startThreshold)/refInt);
    
end
