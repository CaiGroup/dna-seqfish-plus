function points = detectdotsperbarcodev3(imageStruct, threshold, numRounds, numChannels, varargin)
% process the points for each barcode round, and pseudo channel round for
% barcode experiments.
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
% outputs: points is structure of 'round' and 'channels' with x, y, and z
% location of a point
%
% Update:
% 1. Made options to gaussian fit points, use radial center or to
% use compressed sensing.
% 2. Added option to get point locations
%
% Author: Nico Pierson
% Date: 4/16/2019
% Email: nicogpt@caltech.edu
% Updated: 8/4/2019

    %% Set up optional Parameters
    argsLimit = 5;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:detectdotsperbarcodev3:TooManyInputs', ...
            'requires at most 5 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs > 0
        if ~ischar(varargin{1}) 
            error('src:detectdotsperbarcodev2:WrongInput', ...
                'detectdotsperbarcodev3 var detectTypePoints requires type string');
        elseif ~strcmp(varargin{1}, 'exons') && ~strcmp(varargin{1}, 'introns') ...
                && ~strcmp(varargin{1}, 'exons2d') && ~strcmp(varargin{1}, 'log') ...
                && ~strcmp(varargin{1}, 'log2d')
            error('src:detectdotsperbarcodev2:WrongInput', ...
                'detectdotsperbarcodev3 var refinePointLocation requires type string: "exons" or "introns" or "exons2d"');
        end
    end
    if numvarargs > 1
        if ~ischar(varargin{2}) 
            error('src:detectdotsperbarcodev2:WrongInput', ...
                'detectdotsperbarcodev3 var refinePointLocation requires type string');
        elseif ~strcmp(varargin{2}, 'none') && ~strcmp(varargin{2}, 'gaussian') && ~strcmp(varargin{2}, 'radial') && ~strcmp(varargin{2}, 'radial3d')
            error('src:detectbarcodepoints:WrongInput', ...
                'detectdotsperbarcodev3 var refinePointLocation requires type string: "none" or "gaussian" or "radial" or "compressed"');
        end
    end
    if numvarargs > 2
        if ~ischar(varargin{3}) 
            error('src:detectdotsperbarcodev3:WrongInput', ...
                'processimages pointsSavePath requires type string');
        end
    end

    % set defaults for optional inputs
    optargs = {'exons', 'none', [], 1, []};
    
    % assign defaults
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [detectTypePoints, refinePointLocation, pointsSavePath, multiplier, adjustedthreshold] = optargs{:};
    
    

    %% Initialize Variables
    saveFig = true;
    
    barcodeArray = 1:numRounds;
    channelArray = 1:numChannels;
    points = cell(1, numRounds);
    
    
    %% Find the points using the Threshold Values
    for barcode = barcodeArray
        fprintf('Finding Dots in Barcode Round: %.0f channel ', barcode);
        points{barcode} = struct('channels', cell(numChannels, 1), 'intensity', cell(numChannels, 1));
        %points{barcode} = struct('scaledIntensity', cell(numChannels, 1));
        for ch = channelArray
            fprintf('%.0f ', ch);
            hybIndex = (barcode - 1) * numChannels + ch; 
            figSavePath = [pointsSavePath '-Round' num2str(barcode) '-Hyb' num2str(ch) '.fig'];
            
            % need the adjusted threshold if it is used
            if ~isempty(adjustedthreshold)
                threshold = adjustedthreshold;
            end
            % check if threshold is in the right format
            if size(threshold,1) > numChannels
                threshold = orgthrehsold2round(threshold, numRounds, numChannels);
            end
            
            [points{barcode}(ch).channels, points{barcode}(ch).intensity, ~, ~] = detectdotsv2(imageStruct{hybIndex}, ...
                threshold(barcode, ch), detectTypePoints, saveFig, figSavePath, multiplier);
        
            %% Create options for different kind of dot localization - none, gaussian, radialcenter, compression
            switch refinePointLocation
                case 'none'
                    % do nothing
                case 'gaussian'
                    points{barcode}(ch).channels = getgaussian(points{barcode}(ch).channels, imageStruct{hybIndex});
                case 'radial'
                    [points{barcode}(ch).channels, points{barcode}(ch).intensity] = getradialcenter(points{barcode}(ch).channels, imageStruct{hybIndex});
                case 'radial3d'
                    [points{barcode}(ch).channels, points{barcode}(ch).intensity] = SuperResPoints(points{barcode}(ch).channels, imageStruct{hybIndex},1,1);
                otherwise
                    error 'wrong optional input';
            end
            
            %% get scaled intensity - Not sure if needed
            %points{barcode}(ch).scaledIntensity = double(points{barcode}(ch).intensity)/mean(points{barcode}(ch).intensity);

        
        end
        fprintf('\n');
    end
    
    % make the video of the points
    saveVideo = false;
    if saveVideo % turn off for now: currently debugging
        displayRange = [0 300]; % how to change this; need to update
        makevideobarcode(imageStruct, points, numRounds, numChannels, threshold, displayRange, pointsSavePath);
    end
    
end

