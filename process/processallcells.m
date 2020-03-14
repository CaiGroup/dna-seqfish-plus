function [finalPosList, dotlocations, numpointconsensus, numdotlocations, ...
        numfinalpoints, numpointspercell, seeds] = processallcells(projectDir, projectName, ...
        points, position, numRounds, numChannels, barcodekey, segmentPath, varargin)
% processes each cell, decoding the points, and returning the seeds,
% dotlocations, numtotaldotsconsensus, averagepointsperround
% 
% Requirements: 
% 1. barcode key in excel - in 'barcodekey' folder in project directory
% 2. segmentPath for roi or 3d labeled images
%
% Inputs: project directory path, project name will be used to save data
% and directories, position (aka field of view), number of barcode rounds,
% number of channels, directory with full codebook, and path to segment
% file.
%
% Optional Parameters: 1. type of segmentation
%
% Default: 'roi' ['roi', '3d']
%
% To do: 
% 1. Clean up code - streamline and clear code
% 2. Add segmentation for labeled images from ilastik (3d)
% 3. add print functions
% 4. To parallelize: use two functions, one for segmenting, and one for
% decoding per cell
%
% Date: 8/8/2019
% Author: Nico Pierson

    %% Set up optional Parameters
    argsLimit = 7;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:processallcells:TooManyInputs', ...
            'requires at most 4 optional inputs');
    end
    % Error for type of arguments
    if numvarargs == 1 
        if ~ischar(varargin{1}) 
            error('src:processallcells:WrongInput', ...
                'processallcells var segment requires type string');
        elseif ~strcmp(varargin{1}, 'roi') && ~strcmp(varargin{1}, '3d') && ~strcmp(varargin{1}, '2d')
            error('myfun:processallcells:WrongInput', ...
                'processallcells var segment requires type string: "roi" or "3d" or "2d"');
        end
    end
    if numvarargs > 1
        if ~isnumeric(varargin{2}) && ~isscalar(varargin{2})
            error('src:processimages:WrongInput', ...
                'processallcells sqrtradius requires type int');
        end
    end
    if numvarargs > 2
        if ~ischar(varargin{3}) && ~isscalar(varargin{3})
            error('src:processimages:WrongInput', ...
                'processimages var alloweddiff requires type int');
        end
    end
    if numvarargs > 3
        if ~ischar(varargin{4}) 
            error('src:processallcells:WrongInput', ...
                'processallcells var saveoption requires type string');
        elseif ~strcmp(varargin{4}, 'matrix') && ~strcmp(varargin{4}, 'list')
            error('myfun:processallcells:WrongInput', ...
                'processallcells var saveoption requires type string: "matrix" or "list"');
        end
    end
    if numvarargs > 4
        if ~ischar(varargin{5}) 
            error('src:processallcells:WrongInput', ...
                'processallcells var explabel requires type string');
        end
    end
    % set defaults for optional inputs
    optargs = {'roi', 6, 1, 'matrix', [], [], ''};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [segment, sqrtradius, alloweddiff, saveoption, explabel, minseeds, experimentLabel] = optargs{:};

    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    dateSaveString = datestr(dateStart, formatDate);

    
    
    %% initialize the variables
    % Directories and Paths
    if isempty(explabel)
        explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius)];
    end
    % Directories for saving the csv files for each position
    analysisFolder = 'analysis';
    %posFolder = ['pos' num2str(position)];
    savePath = fullfile(projectDir, analysisFolder, experimentLabel, explabel);
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    projectSaveName = [projectName '-Pos' num2str(position) '-' explabel '-' dateSaveString];

    
    
    %% Segment the cells
    % returns pointsxcell in cell by 1 cell array
    if exist(segmentPath, 'file')
        switch segment
            case 'roi'
                % use ROIs from ImageJ to segment cells
                [pointsxcell, numCells, numpointspercell] = segmentpoints2cells(segmentPath, points, segment);
            case '3d'
                % use labeled images from ilastik
                [pointsxcell, numCells, numpointspercell] = segmentpoints2cells(segmentPath, points, segment);
            case 'whole'
                [pointsxcell, numCells, numpointspercell] = segmentpoints2cells(segmentPath, points, segment);
            
            otherwise
                error 'myfun:processallcells:WrongInput var segment invalid argument';
        end
    else
        switch segment
            case '2d'
                [pointsxcell, numCells, numpointspercell] = segmentpoints2cells(segmentPath, points, segment);
            case 'whole'
                % process the whole field of view
                [pointsxcell, numCells, numpointspercell] = segmentpoints2cells(segmentPath, points, segment);
            otherwise
                segment = 'whole';
                % process the whole field of view
                [pointsxcell, numCells, numpointspercell] = segmentpoints2cells(segmentPath, points, segment);
        end
    end



    %% Decode Points for each ROI or labeled cell
    posList = cell(numCells, 1);
    finalPosList = posList;
    dotlocations = posList;
    numpointconsensus = posList;
    numdotlocations = posList;
    seeds = posList;
    numfinalpoints = posList;
    fprintf('Cell:');
    for c = 1:numCells
        fprintf(' %.0f', c);
        [finalPosList{c}, posList{c}, dotlocations{c}, numpointconsensus{c},  numdotlocations{c}, ...
            numfinalpoints{c}, seeds{c}] = decodeimages(pointsxcell{c}, ...
            barcodekey, numRounds, numChannels, alloweddiff, sqrtradius, minseeds);
    end
    fprintf('\n');
    
    
    %% Organize the cell arrays to output data
    outputdata(savePath, projectSaveName, barcodekey, finalPosList, posList, dotlocations, ...
        numpointconsensus,  numdotlocations, numpointspercell, numfinalpoints, seeds, saveoption);

end