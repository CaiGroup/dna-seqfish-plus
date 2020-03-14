function [adjustedThreshold, adjustedIntensity, adjustedPoints] = ...
    outputadjustedthreshold(I, medianIntensity, threshold, folderArray, ...
    numCh, typedots, varargin)
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
    
    
    
    %% get the adjusted threshold for the log filter for pos 0
    adjustedPoints = cell(length(folderArray),numCh);
    adjustedIntensity = adjustedPoints;
    medianError = adjustedPoints;
    adjustedThreshold = ones(length(folderArray),numCh) * 999999;
    for ch = 1:numCh
        parfor folder = 1:length(folderArray)
            [adjustedPoints{folder,ch}, adjustedIntensity{folder,ch}, adjustedThreshold(folder,ch), ...
                medianError{folder,ch}] = adjustthreshold(I{folder,ch}, medianIntensity{folder,ch}, ...
                threshold(folder,ch), typedots, roimask);
        end
    end
end