function [finalPosList, dotlocations, numpointconsensus, numdotlocations, ...
        numfinalpoints, numpointspercell, seeds, points] = ...
        processimages(experimentDir, experimentName, position, ...
        numRounds, numChannels, processedImages, segment, varargin)
% processes the images by getting the points, decoding the points, and
% finding the false positive rate
% 
% Requirements: 
% 1. processed images (image x 1 cell array)
% 2. thresholds (round x channel matrix) - in 'threshold' folder in project
% directory
% 3. barcode key in excel - in 'barcodekey' folder in project directory
% 4. RoiSet.zip for 2d segmentation - in project directory
%
% Inputs: project directory path, project name will be used to save data
% and directories, position (aka field of view), number of barcode rounds,
% number of channels, directory with full codebook, and processed Images.
%
% Optional Parameters: 1. type of spot detection 2. type of super resolve 
% for dot location 3. path for saving points 4. allowed diff 5. square
% root radius
%
% Default: 'exons' ['exons', 'introns', 'exons2d'], 'none' ['none',
% 'gaussian', 'radialcenter'], 1, 6
%
% To do: 
% 1. Clean up code - streamline and clear code
% 2. Add segmentation for labeled images from ilastik
% 3. 3d functions for segmentation
% 4. Test for parallelization
% 5. Update Readme.txt and make it clear and concise
% 6. Add an example
% 7. preprocessing package: Signal Retention for the repeat hyb
%                           Add the pixel shift for alignmnet
%
% Update:
% 1. 8/14/2019 - added optional parameters for processing by channel from
% experiments formatted similarly to seqFISH+
%
% Date: 8/3/2019
% Author: Nico Pierson

    %% Set up optional Parameters
    argsLimit = 11;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:processimages:TooManyInputs', ...
            'requires at most 11 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 0
        if ~isnumeric(varargin{1}) && ~isscalar(varargin{1})
            error('src:processimages:WrongInput', ...
                'processmages sqrtradius requires type int');
        end
    end
    if numvarargs > 1
        if ~ischar(varargin{2}) 
            error('src:processimages:WrongInput', ...
                'processimages var typedots requires type string');
        elseif ~strcmp(varargin{2}, 'exons') && ~strcmp(varargin{2}, 'introns') ...
                && ~strcmp(varargin{2}, 'log') && ~strcmp(varargin{2}, 'log2d')
            error('src:processimages:WrongInput', ...
                'processimages var typedots requires type string: "exons" or "introns" or "log" or "log2d"');
        end
    end
    if numvarargs > 2
        if ~ischar(varargin{3}) 
            error('src:processimages:WrongInput', ...
                'processimages superresolve requires type string');
        elseif ~strcmp(varargin{3}, 'none') && ~strcmp(varargin{3}, 'gaussian') && ~strcmp(varargin{3}, 'radial') && ~strcmp(varargin{3}, 'radial3d')
            error('src:processimages:WrongInput', ...
                'processimages var superresolve requires type string: "none" or "gaussian" or "radial" or "radial3d"');
        end
    end
    if numvarargs > 3
        if ~isnumeric(varargin{4}) && ~isscalar(varargin{5})
            error('src:processimages:WrongInput', ...
                'processimages var alloweddiff requires type int');
        end
    end
    if numvarargs > 4
        if ~isnumeric(varargin{5}) && ~isscalar(varargin{6})
            error('src:processimages:WrongInput', ...
                'processimages var channel requires type int');
        end
    end
    if numvarargs > 6
        if ~ismatrix(varargin{7})
            error('src:processimages:WrongInput', ...
                'processimages var removePoints requires type matrix');
        end
    end
    %if numvarargs > 7
    %    if isempty(varargin{8})
    %        error('src:processimages:WrongInput', ...
    %            'processimages var dotlocations is required');
    %    end
    %end
    if numvarargs > 8
        if ~ismatrix(varargin{9})
            error('src:processimages:WrongInput', ...
                'processimages var removeInd requires type matrix');
        end
    end
    % set defaults for optional inputs
    optargs = {6, 'exons', 'none', 1, [], [], [], [], [], 1, []};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [sqrtradius, typedots, superresolve, alloweddiff, channel, minseeds, ...
        removePoints, dotlocations, removeInd, multiplier, adjustedthreshold] = optargs{:}; % beadPoints used to remove bead points

    
    
    %% Read the Barcode excel file
    barcodeFolder = 'barcodekey';
    if isempty(channel)
        barcodekeyPath = getfile(fullfile(experimentDir,barcodeFolder), 'barcode');
    else
        channelFolder = ['ch' num2str(channel)];
        barcodekeyPath = getfile(fullfile(experimentDir,barcodeFolder, channelFolder), 'barcode');
    end
    barcodekey = readbarcode(barcodekeyPath, 'header');

 
    
    if isempty(removePoints) && isempty(dotlocations)
        %% Process Points
        % Need to add returning variable for number of points in each round
        [points] = processpoints(experimentDir, experimentName, position, numRounds, ...
            numChannels, sqrtradius, alloweddiff, processedImages, typedots, ...
            superresolve, channel, multiplier, adjustedthreshold);
    else
        %% Remove Points - make function to remove points within a certain radius
        % need the previous points and the dotlocations, with indices to remove
        % certain points
        points = removepoints(removePoints, dotlocations, removeInd);
    end


    %% Decode Points for each ROI or labeled cell
    % Get the path for segmentation - need to make useful for 3d
    segmentFolder = 'segmentation';
    segmentPath = fullfile(experimentDir, segmentFolder, ['Pos' num2str(position)], 'RoiSet.zip');
    %segment = 'roi'; % can use ['roi', '3d', 'whole']
    saveoption = 'list'; % can use ['matrix', 'list'] % does both, decide what to do with the options
    if isempty(channel)
        explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius)];
    else
        explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius) '-ch' num2str(channel)];
    end
    [finalPosList, dotlocations, numpointconsensus, numdotlocations, ...
        numfinalpoints,  numpointspercell, seeds] = processallcells(experimentDir, experimentName, ...
        points, position, numRounds, numChannels, barcodekey, segmentPath, ...
        segment, sqrtradius, alloweddiff, saveoption, explabel, minseeds);
    
    
    
    %{
    %% False Positive Rate
    % Get the full barcode key
    barcodekeyfull = readbarcode(fullcodebookDir, 'only data');
    if ~isempty(barcodekeyfull)
        [finalPosListFP, dotlocationsFP, numtotalFP, seedsFP, barcodekeyfull] ...
            = falsepositive(finalPosList, dotlocations, points, avgpointsperround,...
            barcodekeyfull, savePath, projectName, numRounds, numChannels, ...
            alloweddiff, sqrtradius);
    end
    %}

end