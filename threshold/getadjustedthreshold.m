function adjustedThreshold = getadjustedthreshold(chArray, numRounds, numChannels, ...
    pointsDir, pointsName, I, typedots, varargin)
% get the adjusted threshold for a different position by using the manually
% thresholded points from the reference position
%
% make function to get the adjusted threshold an

%% Set up optional Parameters
    argsLimit = 1;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:getadjustedthreshold:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end   
    % set defaults for optional inputs
    optargs = {[]};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [roimask] = optargs{:};

%% work on points and intensity and getting the output
    % how to enter the refmedian
    %chArray = 1:2;

    %numRounds = 5;
    %numChannels = 16;
    %position = 0;
    %pointsDir = fullfile(experimentDir, 'analysis', ['2error-sqrt6-ch' num2str(c)]);
    %pointsName = ['points-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped-Pos' num2str(position) '-2019-11-06.mat'];
    numHybs = numRounds * numChannels;
    numCh = length(chArray);
    [refPoints, refIntensity, threshold] = mat2hybxcell(pointsDir, pointsName, chArray, numRounds, numChannels);
    % need function to get the points and the intensity for the firs position
    % that is manually thresholded


    % get the median ref intensity
    medianIntensity = refIntensity;
    numRefPoints = cell(numHybs, numCh);
    for i = 1:size(medianIntensity,1)
        for ch = chArray
            medianIntensity{i,ch} = median(medianIntensity{i,ch});
            numRefPoints{i,ch} = size(refPoints{i,ch},1);
        end
    end
                
                
    %% get the points
    points = cell(numHybs,numCh);
    intensity = points;
    adjustedThreshold = ones(numHybs,numCh) * 999999;
    for ch = 1:numCh
        parfor f = 1:numHybs
            [points{f,ch}, intensity{f,ch}, adjustedThreshold(f,ch)] = ...
                adjustthreshold(I{f,ch}, medianIntensity{f,ch}, numRefPoints{f,ch}, ...
                threshold(f,ch), typedots, roimask);
        end
    end
    
end