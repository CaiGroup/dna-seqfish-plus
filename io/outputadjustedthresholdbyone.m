function [adjustedThreshold, adjustedIntensity, adjustedPoints] = ...
    outputadjustedthresholdbyone(image, medianIntensity, numRefPoints, startThreshold, ...
    typedots, varargin)
% output adjusted threshold


    %% Set up optional Parameters
    argsLimit = 1;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:outputadjustedthreshold:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end   
    % set defaults for optional inputs
    optargs = {[]};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [roimask] = optargs{:};
    
    
    
    %% initial median - use the number of points ratio as the mutliplier
    %refInt = medianIntensity;
    [adjustedPoints, adjustedIntensity, ~, ~] = detectdotsv2(image, startThreshold, typedots, false, '', 1);
    % only get the median intensity of the points inside the rois
    removeInd = [];
    if ~isempty(roimask)
        numPoints = size(adjustedPoints,1);
        for i = 1:numPoints
            if ~roimask(adjustedPoints(i,2), adjustedPoints(i,1), adjustedPoints(i,3))
                removeInd = cat(1, removeInd, i);
            end
        end
        adjustedPoints(removeInd,:) = [];
        adjustedIntensity(removeInd) = [];
    end
    
    %adjustedThreshold = round((medInt*startThreshold)/refInt);
    numAdjustedPoints = size(adjustedPoints,1);
    adjustedThreshold = startThreshold * (double(median(adjustedIntensity))/double(medianIntensity));
   

end