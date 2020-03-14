function points = detectdotsperhybcycle(experimentDir, imageStruct, folderArray, channel, varargin)
% process the points for each hyb round for sequential experiments
%
% Requirements: 
% 1. images in cell hybs x 1 cell array
% 2. threshold is a m x n matrix where m is the number of barcode
% rounds and n is the number of pseudochannel rounds
%
% inputs: processed image, threshold (integer), 
% optional inputs: flags for viewing the dot images, save the image, and
% path for saving.
%
% outputs: points is structure of 'channels' with x, y, and z
% location of a point, 'intensity' for raw intensity, and 'scaledIntensity'
% for scaled intensity using the mean
%
% Author: Nico Pierson
% Date: 8/27/2019
% Email: nicogpt@caltech.edu


    %% Set up optional Parameters
    argsLimit = 3;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:detectdotsperhybcycle:TooManyInputs', ...
            'requires at most 3 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs > 0
        if ~ischar(varargin{1}) 
            error('src:detectdotsperhybcycle:WrongInput', ...
                'detectdotsperhybcycle var detectTypePoints requires type string');
        elseif ~strcmp(varargin{1}, 'exons') && ~strcmp(varargin{1}, 'introns') && ~strcmp(varargin{1}, 'exons2d')
            error('src:detectdotsperhybcycle:WrongInput', ...
                'detectdotsperhybcycle var refinePointLocation requires type string: "exons" or "introns" or "exons2d"');
        end
    end
    if numvarargs > 1
        if ~ischar(varargin{2}) 
            error('src:detectdotsperhybcycle:WrongInput', ...
                'detectdotsperhybcycle var refinePointLocation requires type string');
        elseif ~strcmp(varargin{2}, 'none') && ~strcmp(varargin{2}, 'gaussian') && ~strcmp(varargin{2}, 'radial') && ~strcmp(varargin{2}, 'compressed')
            error('src:detectbarcodepoints:WrongInput', ...
                'detectdotsperhybcycle var refinePointLocation requires type string: "none" or "gaussian" or "radial" or "compressed"');
        end
    end
    if numvarargs >= argsLimit
        if ~ischar(varargin{3}) 
            error('src:detectdotsperhybcycle:WrongInput', ...
                'processimages pointsSavePath requires type string');
        end
    end

    % set defaults for optional inputs
    optargs = {'exons', 'none', []};
    
    % assign defaults
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [detectTypePoints, refinePointLocation, pointsSavePath] = optargs{:};
    
    

    %% Initialize Variables
    saveFig = true;
    numHybs = length(folderArray);
    points = cell(numHybs, 1);
    
    
    %% Find the points using the Threshold Values
    for f = 1:numHybs % starts at 0
        hybIndex = f + 1;
        fprintf('Finding Dots in HybCycle: %.0f channel', folderArray(f));
        points{f} = struct('channels', cell(numHybs, 1));
        points{f} = struct('intensity', cell(numHybs, 1));
        points{f} = struct('scaledIntensity', cell(numHybs, 1));
        for ch = channel
            fprintf('%.0f ', ch);
            
            
            %% Get threshold per channel
            thresholdFolder = 'threshold';
            channelFolder = ['ch' num2str(ch)];
            thresholdPath = fullfile(experimentDir, thresholdFolder, channelFolder);
            if exist(thresholdPath, 'dir') ~= 0
                thresholdPathFind = getfile(thresholdPath, 'threshold');
                t = load(thresholdPathFind, 'threshold');
                threshold = t.threshold;
            else
                error('detectdotsperhybcycle:Threshold not found at: %s', thresholdPath);
            end
            
            
            
            figSavePath = [pointsSavePath '-Hyb' num2str(folderArray(f)) 'Ch' num2str(ch) '.fig'];
            [points{f}(ch).channels, points{f}(ch).intensity, ~, ~] = detectdotsv2(imageStruct{hybIndex, ch}, ...
                threshold(f), detectTypePoints, saveFig, figSavePath);
        
            %% Create options for different kind of dot localization - none, gaussian, radialcenter, compression
            switch refinePointLocation
                case 'none'
                    % do nothing
                case 'gaussian'
                    points{f}(ch).channels = getgaussian(points{f}(ch).channels, imageStruct{hybIndex, ch});
                case 'radial'
                    points{f}(ch).channels = getradialcenter(points{f}(ch).channels, imageStruct{hybIndex, ch});
                case 'compressed'
                    % need to add this part
                otherwise
                    error 'wrong optional input';
            end
            
            %% get scaled intensity - Not sure if needed
            points{f}(ch).scaledIntensity = double(points{f}(ch).intensity)/mean(points{f}(ch).intensity);

        
        end
        fprintf('\n');
    end
    
    % make the video of the points
    if ~saveFig % turn off for now: currently debugging
        displayRange = [0 300]; % how to change this; need to update
        % need to update for making a video using only hyb images
        makevideobarcode(imageStruct, points, numRounds, numChannels, threshold, displayRange, pointsSavePath);
    end
    
end


