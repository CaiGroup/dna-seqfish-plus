function [finalPosList, dotlocations, numpointconsensus, numdotlocations, ...
        numfinalpoints, numpointspercell, seeds, points] = ...
        processimagespoints(experimentDir, experimentName, position, ...
        numRounds, numChannels, points, segment, varargin)
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
    argsLimit = 9;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:processimages:TooManyInputs', ...
            'requires at most 9 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 0
        if ~isnumeric(varargin{1}) && ~isscalar(varargin{1})
            error('src:processimages:WrongInput', ...
                'processmages sqrtradius requires type int');
        end
    end
    if numvarargs > 1
        if ~isnumeric(varargin{2}) && ~isscalar(varargin{2})
            error('src:processimages:WrongInput', ...
                'processimages var alloweddiff requires type int');
        end
    end
    if numvarargs > 2
        if ~isnumeric(varargin{3}) && ~isscalar(varargin{3})
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
    if numvarargs > 7
        if isempty(varargin{8})
            error('src:processimages:WrongInput', ...
                'processimages var dotlocations is required');
        end
    end
    if numvarargs >= argsLimit
        if ~ismatrix(varargin{9})
            error('src:processimages:WrongInput', ...
                'processimages var removeInd requires type matrix');
        end
    end
    % set defaults for optional inputs
    optargs = {6, 1, [], [], '', [], [], []};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [sqrtradius, alloweddiff, channel, minseeds, experimentLabel, removePoints, dotlocations, removeInd] = optargs{:}; % beadPoints used to remove bead points

    
    
    %% Read the Barcode excel file
    barcodeFolder = 'barcodekey';
    if isempty(channel)
        barcodekeyPath = getfile(fullfile(experimentDir,barcodeFolder), 'barcode');
    else
        channelFolder = ['ch' num2str(channel)];
        barcodekeyPath = getfile(fullfile(experimentDir,barcodeFolder, channelFolder), 'barcode');
    end
    barcodekey = readbarcode(barcodekeyPath, 'header');

 


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
        segment, sqrtradius, alloweddiff, saveoption, explabel, minseeds, experimentLabel);
    
    
    
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