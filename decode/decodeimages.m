function [finalPosList, PosList, dotlocations, numpointconsensus, numdotlocations, ...
    numfinalpoints, seeds] = decodeimages(points, barcodekey, numRounds, numChannels, varargin)
% function decodes the points using the barcodekey for each set of ROI
% polygon (cell, region), or each 3d polygon.
% 
% Requirements: 
% 1. barcode key - in 'barcodekey' folder in project directory
% 2. roi.zip - in project directory
%
% Inputs: project directory path, project name will be used to save data
% and directories, position (aka field of view), 
%
% To do: 
% 1. add rois for 3d labeled images
% 2. Optional arguments for radius and alloweddiff
% 3. Add other csv matrices for saving: seeds, PosList, etc.
% 4. Generate the save path elsewhere
% 5. Create functions for efficiency
%
% Optional Parameters: 1. allowed difference or number of minimum barcodes
% to be dropped 2. square root radius for decoding points
% points.
%
% Default: 1 and 6
%
% Date: 8/3/2019
% Author: Nico Pierson
 
    %% Set up optional Parameters
    argsLimit = 3;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:decodeimages:TooManyInputs', ...
            'requires at most 2 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 0
        if ~ischar(varargin{1}) && ~isscalar(varargin{1})
            error('src:decodeimages:WrongInput', ...
                'decodeimages var typedots requires type int');
        end
    end
    if numvarargs > 1
        if ~isnumeric(varargin{2}) && ~isscalar(varargin{2})
            error('src:decodeimages:WrongInput', ...
                'decodeimages superresolve requires type int');
        end
    end
    if numvarargs > 2
        if ~isnumeric(varargin{3}) && ~isscalar(varargin{3})
            error('src:decodeimages:WrongInput', ...
                'decodeimages minseeds requires type int');
        end
    end
    % set defaults for optional inputs
    optargs = {1, 6, []};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [alloweddiff, sqrtradius, minseeds] = optargs{:};



    
    %% initialize the variables
    if isempty(minseeds)
        minseeds = numRounds - 1;
    end

    
    

    %% Find the barcodes
    barcodekeycode = barcodekey.barcode;
    [dotlocations,seeds] = BarcodeFinder(numChannels, points, numRounds, barcodekeycode,sqrtradius,alloweddiff);


    

    %% Filter the Based on the number of Seeds 
    [finalPosList, PosList, dotlocations, numdotlocations, numpointconsensus, numfinalpoints] = filterseedsv2(seeds, dotlocations, minseeds);

    
end

