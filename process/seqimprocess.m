function [pointspercell, points] = seqimprocess(experimentDir, experimentName, ...
    processedImages, position, folderArray, channelArray, segoption, varargin)
% processes the images by getting the points, assigning to cells, then
% outputing to csv files.
% 
% Requirements: 
% 1. processed images (image x 1 cell array)
% 2. thresholds (round x channel matrix) - in 'threshold' folder in project
% directory
% 3. barcode key in excel - in 'barcodekey' folder in project directory
% 4. RoiSet.zip for 2d segmentation - in project directory
%
% Inputs: project directory path, project name will be used to save data
% and directories, position (aka field of view), folder array, channel 
% array, and processed images.
%
% Assumes threshold is in 'experimentDir\threshold\ch1\*.mat' directory for
% channel 1 threshold as a hybcycle by 1 matrix.
%
% Optional Parameters: 
% 1. square root radiustype of spot detection 
% 2. type of super resolve for dot location 
% 3. path for saving points 
%
% Default: 
% 1. square root radius for colocalization
% 2. 'exons' ['exons', 'introns', 'exons2d']
% 3. 'none' ['none', 'gaussian', 'radialcenter']
%
% To do: 
% 1. segoption '3d' for 3d segmentation
%
% Date: 8/27/2019
% Author: Nico Pierson

    %% Set up optional Parameters
    argsLimit = 3;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:processimages:TooManyInputs', ...
            'requires at most 3 optional inputs');
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
        elseif ~strcmp(varargin{2}, 'exons') && ~strcmp(varargin{2}, 'introns') && ~strcmp(varargin{2}, 'exons2d')
            error('src:processimages:WrongInput', ...
                'processimages var typedots requires type string: "exons" or "introns" or "exons2d"');
        end
    end
    if numvarargs > 2
        if ~ischar(varargin{3}) 
            error('src:processimages:WrongInput', ...
                'processimages superresolve requires type string');
        elseif ~strcmp(varargin{3}, 'none') && ~strcmp(varargin{3}, 'gaussian') && ~strcmp(varargin{3}, 'radial') && ~strcmp(varargin{3}, 'compressed')
            error('src:processimages:WrongInput', ...
                'processimages var superresolve requires type string: "none" or "gaussian" or "radial" or "compressed"');
        end
    end


    % set defaults for optional inputs
    optargs = {6, 'exons', 'none'};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [sqrtradius, typedots, superresolve] = optargs{:}; % beadPoints used to remove bead points

    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = datestr(dateStart, formatDate);
    
    
    %% Set up Directories and Save names
    explabel = ['sqrt' num2str(sqrtradius)];
    analysisFolder = 'analysis';
    savePath = fullfile(experimentDir, analysisFolder, explabel);
    savePointsCheck = fullfile(savePath, 'pointscheck');
    if ~exist(savePointsCheck, 'dir')
        mkdir(savePointsCheck);
    end
    saveFileNameEnd = [experimentName '-Pos' num2str(position) '-' dateSaveString];
    saveFileName = fullfile(savePath, ['pointsData-' saveFileNameEnd '.mat']); % Main data to save
    pointsSavePath = fullfile(savePointsCheck, ['PointsCheck-' experimentName '-Pos' num2str(position)]);


    %% Retrieve Points, Organize in Barcode Rounds and Colocalize the points
    % processedImages should be a hyb by channel cell array for the number of
    % folders
    points = detectdotsperhybcycle(experimentDir, processedImages, folderArray, channelArray, typedots, superresolve, pointsSavePath);
    
     


    %% Segment points based on cells and output
    segPath = fullfile(experimentDir, 'segmentation', ['ch' num2str(channelArray)]);
    [pointspercell, numCells] = segment2dcells(segPath, channelArray, position, points);
    
    
    
    %% Output the points per cell
    savePointsFilePath = fullfile(savePath, ['pointsList-' saveFileNameEnd '.csv']); % Main data to save
    outputpointspercell(pointspercell, savePointsFilePath, position)

    % save points and data
    save(saveFileName, 'points', 'threshold', 'numCells', 'pointspercell', 'channelArray', 'folderArray', 'position');
end

function [] = outputpointspercell(pointspercell, savePath, position)
    % variables and directory
    numCells = size(pointspercell,1);
    numHybs = size(pointspercell{1},1);
    numChannels = numel(pointspercell{1});

    % set up csv file
    fileID = fopen(savePath,'w');
    fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                    'fov', 'cellID', 'hybID', 'chID', 'x', 'y', 'z', 'int');
    fprintf(fileID,'\n'); % first new line doesn't create a new line
    
    for c = 1:numCells
        for f = 1:numHybs
            for ch = 1:numChannels
                numPoints = size(pointspercell{c}{f}(ch).channels,1);
                for p = 1:numPoints
                    pointSelection = pointspercell{c}{f};
                    x = pointSelection(ch).channels(p,1);
                    y = pointSelection(ch).channels(p,2);
                    z = pointSelection(ch).channels(p,3);
                    int = pointSelection(ch).intensity(p);
                    fprintf(fileID,'%d,%d,%d,%d,%.3f,%.3f,%.3f,%.1f\n', ...
                                            position, c, f, ch, x, y, z, int);
                end
            end
        end
    end

end

function [pointspercell, numCells] = segment2dcells(segPath, channelArray, position, points)
    roiSegPath = fullfile(segPath, ['Pos' num2str(position)], 'RoiSet.zip');
    vertex = selfsegzip(roiSegPath);
    numCells = numel(vertex);
    numChannels = length(channelArray);
    pointspercell = cell(numCells, 1);
    for c = 1:numCells
        pointspercell{c} = points;
        numHybs = size(points, 1);
        for f = 1:numHybs
            for ch = 1:numChannels
                pointsSelection = points(ch).channels;
                ind = inpolygon(pointsSelection(:,1), pointsSelection(:,2), vertex(c).x, vertex(c).y);
                pointspercell{c}{f}(ch).channels = pointspercell{c}{f}(ch).channels(ind,:);
                pointspercell{c}{f}(ch).intensity = pointspercell{c}{f}(ch).intensity(ind,:);
                pointspercell{c}{f}(ch).scaledIntensity = pointspercell{c}{f}(ch).scaledIntensity(ind,:);
            end
        end
        
    end
end