function [points, threshold, matchFinal] = matchpointsims(image, varargin)
%function [points, tformCh] = matchpointsallch(imagePath, compareString, numOfChannels)
% create function to compare points acrosss all channels
%
% Output: returns structure match and mismatch with nested matrix of 'ref' and
% 'points'
%
% Update: 
% 1. 8/12/2019 - cleaned function; added optional parameters for threshold,
% back radius and sliding parabolar in imagej. Return threshold that can be
% used for other channels.
%
% To do: 
% 1. Make min threshold function quicker
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 4/3/2019
% Updated: 8/12/2019



    %% Set up optional Parameters
    argsLimit = 11;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:matchpointsims:TooManyInputs', ...
            'requires at most 11 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 0
        if ~isnumeric(varargin{1}) && ~isscalar(varargin{1}) && ~iscell(varargin{1}) && ~ismat(varargin{1})
            error('src:matchpointsims:WrongInput', ...
                'matchpointsims threshold requires type int or cell or mat');
        end
    end
    if numvarargs > 1
        if ~ischar(varargin{2}) 
            error('src:matchpointsims:WrongInput', ...
                'matchpointsims var typedots requires type string');
        elseif ~strcmp(varargin{2}, 'exons') && ~strcmp(varargin{2}, 'introns') && ~strcmp(varargin{2}, 'exons2d') && ~strcmp(varargin{2}, 'log')
            error('src:matchpointsims:WrongInput', ...
                'matchpointsims var typedots requires type string: "exons" or "introns" or "exons2d"');
        end
    end
    if numvarargs > 2
        if varargin{3} ~= 1 && varargin{3} ~= 0
            error('src:matchpointsims:WrongInput', ...
                'matchpointsims process requires type boolean');
        end
    end
    if numvarargs > 3
        if ~isnumeric(varargin{4}) && ~isscalar(varargin{4})
            error('src:matchpointsims:WrongInput', ...
                'pmatchpointsims var backradius requires type int');
        end
    end
    if numvarargs > 4
        if varargin{5} ~= 1 && varargin{5} ~= 0
            error('src:matchpointsims:WrongInput', ...
                'matchpointsims sliding requires type boolean');
        end
    end
    if numvarargs > 5
        if ~isnumeric(varargin{6}) && ~isscalar(varargin{6})
            error('src:matchpointsims:WrongInput', ...
                'pmatchpointsims var numRefPoints requires type int');
        end
    end
    if numvarargs > 6
        if ~isnumeric(varargin{7})
            error('src:matchpointsims:WrongInput', ...
                'pmatchpointsims var searchradius requires type int');
        end
    end
    if numvarargs > 7
        if varargin{8} ~= 1 && varargin{8} ~= 0
            error('src:matchpointsims:WrongInput', ...
                'matchpointsims debug requires type boolean');
        end
    end
    if numvarargs > 8
        if ~ischar(varargin{9}) 
            error('src:matchpointsims:WrongInput', ...
                'matchpointsims var saveFigPath requires type string');
        end
    end
    % set defaults for optional inputs
    optargs = {[], 'exons', true, 3, false, [], sqrt(3), false, pwd, 'gaussian', false};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [threshold, typedots, process, backradius, sliding, target, searchradius, debug, saveFigPath, superres, filtersigma] = optargs{:};

    
    
    %% Declare Variables
    %fijiDirectory = checkfijipath();
    numOfChannels = length(image);
    matchFinal = cell(1, numOfChannels);
    mismatchFinal = cell(1, numOfChannels);
    dotsMoving = cell(1, numOfChannels);
    intensityCh = cell(1, numOfChannels);
    
    
    
    %% Background Subtract Images
    if process
        I = cell(1, numOfChannels);
        uniqueString = 'tempImage-aioeshfoi43534';
        despeckle = true;
        for c = 1:numOfChannels
            I{c} = imagejbackgroundsubtraction(image{c}, ...
                uniqueString, saveFigPath, backradius, sliding, despeckle);
        end
    else
        I = image;
    end
    
    options.processimages = false; % do it before
    options.typedetectdots = typedots; % 'exons' or 'introns'
    options.backradius = backradius;
    options.sliding = sliding;
    
    if isempty(threshold) && isempty(target)     %threshold = 475; % for yodai rep2; 1200 for sliding parabolar
        threshold = ones(1, numOfChannels) * 99999;
        for i = 1:numOfChannels
            threshold(i) = getthreshold(I{i}, options);
        end
    end
    
    
    
    
    %% Get the Reference points with threshold - make function to get doubled the amount of reference points
    % Get the minimum Threshold
    if isempty(threshold)
        threshold = ones(1, numOfChannels) * 99999;
        numIters = 7;
        switch typedots
            case 'exons'
                x1 = 10000;
                x2 = 100000;
            case 'introns'
                x1 = 100; % start at 100
                x2 = 7500;
        end
        for i = 1:numOfChannels
            minThreshold = findminthreshold(x1, x2, numIters, target, I{i}, typedots);
            threshold(i) = minThreshold;
        end
        fprintf('Threshold:, %.2f, %.2f, %.2f\n', threshold);
    end
    
    % check if cell
    if iscell(threshold)
        threshold = cell2mat(threshold);
    end
    % Get points for first channel
    intensity = cell(1, numOfChannels);
    pointsRef = cell(1, numOfChannels);
        filteroption = false;
    [dots, ~, ~, ~] = detectdotsv2(I{1}, threshold(1), typedots);
    switch superres
        case 'gaussian'
            [dotsRef, intensityRef]= getgaussian(dots, image{1});
            [~, ~, sigmaRef] = SuperResPoints(dots,image{1},1,1, filteroption);
            matchFinal{1}.sigma = sigmaRef;
        case 'radial3d'
            [dotsRef, intensityRef, sigmaRef] = SuperResPoints(dots,image{1},1,1);
            matchFinal{1}.sigma = sigmaRef;
            
        case 'radial'
            dotsRef = getradialcenter(dots, image{1});
    end
    matchFinal{1}.channels = dotsRef;
    matchFinal{1}.intensity = intensityRef;
    intensityCh{1} = intensityRef;
    dotsMoving{1} = dotsRef;

    
    
    %% Match points across Channels

    for j = 2:numOfChannels
        [dotsTemp, ~, ~, ~] = detectdotsv2(I{j}, threshold(j), typedots);
        %[dotsMoving{j}, intensityCh{j}] = getgaussian(dotsTemp, image{j});
        %[dotsMoving{j}, intensityCh{j}, sigma] = SuperResPoints(dotsMoving{j}, image{j},1,1); 
        %[~, ~, sigma] = SuperResPoints(dotsTemp, image{j},1,1, filteroption); 
        [dotsMoving{j}, intensityCh{j}, sigma] = SuperResPoints(dotsTemp, image{j},1,1);

        [matchFinal{j}, mismatchFinal{j}, ~, ~] = matchpoints(dotsRef, dotsMoving{j}, searchradius, intensityRef, intensityCh{j}, sigmaRef, sigma);
    end
    
    for ch = 1:numOfChannels
        removeidx = [];
        if ch == 1
            %adjPoints = matches(ch).channels;
            numPoints = size(matchFinal{ch}.channels,1);
            int = matchFinal{ch}.intensity;
        else
            %adjPoints = matches(ch).points;
            numPoints = size(matchFinal{ch}.points,1);
            int = matchFinal{ch}.intmatch;
        end
        sigma = matchFinal{ch}.sigma;

        if filtersigma
            for i = 1:numPoints

                %if sigma(i) < 0.58 || sigma(i) > 0.72 || int(i) < 500 %E14
                if sigma(i) < 0.58 || int(i) < 100 % changed for 48hr
                    removeidx = cat(1, removeidx, i);
                end

            end
            if ch == 1
                matchFinal{ch}.channels(removeidx,:) = [];
                matchFinal{ch}.intensity(removeidx,:) = [];
                matchFinal{ch}.sigma(removeidx,:) = [];
            else
                matchFinal{ch}.points(removeidx,:) = [];
                matchFinal{ch}.intmatch(removeidx,:) = [];
                matchFinal{ch}.sigma(removeidx,:) = [];
                matchFinal{ch}.intref(removeidx,:) = [];
                matchFinal{ch}.ref(removeidx,:) = [];
                matchFinal{ch}.distance(removeidx,:) = [];
            end
        end
    end
    
    %% Organize Points that Match to the same row
    points = getmatchedpoints(matchFinal);
    
    
    
    %% Show figure to check
    if debug
        for i = 1:numOfChannels
            f = figure();
            % View image with dots: can use logFish or fish image
            imshow(max(image{i}, [], 3), 'DisplayRange', [min(min(max(image{i},[],3))) mean(mean(max(image{i}, [], 3))) + 1000], 'InitialMagnification', 'fit');
            hold on;
            scatter(points(i).channels(:, 1), points(i).channels(:, 2), 100, 'mx', 'LineWidth', 1);
            scatter(dotsMoving{i}(:, 1), dotsMoving{i}(:, 2), 100, 'go', 'LineWidth', 1);
            hold off;
            savefig(f, saveFigPath);
            close all;
        end
    end
end

