function [points] = processpoints(experimentDir, experimentName, position, numRounds, numChannels, ...
    sqrtradius, alloweddiff, processedImages, varargin)
% processes the images by getting the points 
% 
% Requirements: 
% 1. processed images (image x 1 cell array)
% 2. thresholds (round x channel matrix) - in 'threshold' folder in project
% directory
%
% Inputs: project directory path, project name will be used to save data
% and directories, position (aka field of view), number of barcode rounds,
% number of channels, radius for point search radius, allowed number of
% barcodes needed, processed images.
%
% Optional Parameters: 1. type of spot detection 2. type of super resolve 
% for dot location 3. boolean for saving points.
%
% Default: 'exons' and 'none'
%
% To do: 
% 1. Clean up code - streamline and clear code
%
% Date: 8/3/2019
% Author: Nico Pierson
 
    %% Set up optional Parameters
    argsLimit = 5;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:processpoints:TooManyInputs', ...
            'requires at most 5 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 0
        if ~ischar(varargin{1}) 
            error('src:processpoints:WrongInput', ...
                'processpoints var typedots requires type string');
        elseif ~strcmp(varargin{1}, 'exons') && ~strcmp(varargin{1}, 'introns') ...
                && ~strcmp(varargin{1}, 'exons2d') && ~strcmp(varargin{1}, 'log') ...
                && ~strcmp(varargin{1}, 'log2d')
            error('src:processpoints:WrongInput', ...
                'processpoints var typedots requires type string: "exons" or "introns" or "exons2d"');
        end
    end
    if numvarargs > 1
        if ~ischar(varargin{2}) 
            error('src:processpoints:WrongInput', ...
                'processpoints superresolve requires type string');
        elseif ~strcmp(varargin{2}, 'none') && ~strcmp(varargin{2}, 'gaussian') && ~strcmp(varargin{2}, 'radial') && ~strcmp(varargin{2}, 'radial3d')
            error('src:processpoints:WrongInput', ...
                'processpoints var superresolve requires type string: "none" or "gaussian" or "radial" or "compressed"');
        end
    end
    if numvarargs > 2
        if ~isnumeric(varargin{3}) && ~isscalar(varargin{3})
            error('src:processpoints:WrongInput', ...
                'processpoints var channel requires type int');
        end
    end
    % set defaults for optional inputs
    optargs = {'exons', 'none', [], 1, []};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [typedots, superresolve, channel, multiplier, adjustedthreshold] = optargs{:};


    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    dateSaveString = datestr(dateStart, formatDate);



    %% initialize the variables
    % Directories and Paths
    if isempty(channel)
        explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius)];
    else
        explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius) '-ch' num2str(channel)];
    end
    % Directories for saving the csv files for each position
    analysisFolder = 'analysis';
    savePath = fullfile(experimentDir, analysisFolder, explabel);
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    saveFileNameEnd = [experimentName '-Pos' num2str(position) '-' dateSaveString];
    saveFileName = [savePath filesep 'points-' saveFileNameEnd '.mat']; % Main data to save



    %% Load threshold
    thresholdFolder = 'threshold';
    if isempty(channel)
        thresholdPath = getfile(fullfile(experimentDir, thresholdFolder), 'threshold');
    else
        channelFolder = ['ch' num2str(channel)];
        thresholdPath = getfile(fullfile(experimentDir,thresholdFolder, channelFolder), 'threshold');
    end
    t = load(thresholdPath, 'threshold');
    threshold = t.threshold;



    %% Retrieve Points, Organize in Barcode Rounds and Colocalize the points
    pointsDir = fullfile(savePath, 'points-error');
    if ~exist(pointsDir, 'dir')
        mkdir(pointsDir);
    end
    pointsSaveDir = fullfile(pointsDir, ['PointsCheck-' experimentName '-Pos' num2str(position)]);
    points = detectdotsperbarcodev3(processedImages, threshold, numRounds, ...
        numChannels, typedots, superresolve, pointsSaveDir, multiplier, adjustedthreshold);

    %% Save points
    if ~isempty(pointsDir)
        save(saveFileName, 'points', 'threshold');
    end
    
end

