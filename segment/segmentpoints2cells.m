function [pointsxcell, numCells, pointsPerCell] = segmentpoints2cells(segmentPath, pointsRef, varargin)
% create a function to divide the points by cell using the RoiSet.zip
% segmentation from Fiji
%
% To do: create a function for 3d cell segmentation, or clean up old one
%
% Example: roiPathName =
% 'I:\FalsePositiveTestBrainLinus-2019-08-06\RoiSet.zip';
%
% Bug Fix: 
% 1. 8/10/2019- points deleting points used and next cell will have not points
%
% Implementation: use the parfor to make it faster
% Date: 8/6/2019
% Author: Nico Pierson

    %% Set up optional Parameters
    argsLimit = 1;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:segmentpoints2cells:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    % Error for type of arguments
    if numvarargs == 1
        if ~ischar(varargin{1}) 
            error('src:segmentpoints2cells:WrongInput', ...
                'segmentpoints2cells var segment requires type string');
        elseif ~strcmp(varargin{1}, 'roi') && ~strcmp(varargin{1}, '3d') && ~strcmp(varargin{1}, 'whole') ...
                && ~strcmp(varargin{1}, '2d')
            error('myfun:segmentpoints2cells:WrongInput', ...
                'segmentpoints2cells var segment requires type string: "roi" or "3d" or "whole"');
        end
    end
    
    % set defaults for optional inputs
    optargs = {'roi'};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [segment] = optargs{:};
    
    sizeZ = max(pointsRef{1}(1).channels(:,3));
    switch segment
        case 'roi'
            
            vertex = selfsegzip(segmentPath);
            numCells = length(vertex);
        case '3d'
            % add code here
        case 'whole'
            numCells = 1;
        case '2d'
            numCells = sizeZ;
        otherwise
            error 'myfun:segmentpoints2cells:WrongInput var segment invalid argument';
    end
    
    
    
    numRounds = length(pointsRef);
    numChannels = length(pointsRef{1});
    pointsxcell = cell(numCells, 1);
    pointsPerCell = pointsxcell;
    points = pointsRef; % use if need to not reset points list
    
    for i = 1:numCells
        %points = pointsRef; % use if want to reset points list
        pointsTemp = cell(1, numRounds);
        pointsPerCell{i} = 0;
        for r = 1:numRounds
            pointsTemp{r} = struct('channels', cell(1, numChannels), 'intensity', cell(1, numChannels));
            %pointsTemp{r} = struct('scaledIntensity', cell(1, numChannels));
            for c = 1:numChannels
                switch segment
                    case 'roi'
                        polygonIn = inpolygon(points{r}(c).channels(:,1), points{r}(c).channels(:,2), vertex(i).x, vertex(i).y);
                        pointsTemp{r}(c).channels = points{r}(c).channels(polygonIn,:);
                        pointsTemp{r}(c).intensity = points{r}(c).intensity(polygonIn,:);
                        %pointsTemp{r}(c).scaledIntensity = points{r}(c).scaledIntensity(polygonIn,:);
                        points{r}(c).channels(polygonIn, :) = []; % remove used points
                        points{r}(c).intensity(polygonIn, :) = []; % remove used points
                        %points{r}(c).scaledIntensity(polygonIn, :) = []; % remove used points
                        pointsPerCell{i} = pointsPerCell{i} + size(pointsTemp{r}(c).channels, 1);
                    case '3d'
                        %add code
                  
                    case 'whole'
                        % points remain the same
                        pointsTemp{r}(c).channels = points{r}(c).channels;
                        pointsTemp{r}(c).intensity = points{r}(c).intensity;
                        %pointsTemp{r}(c).scaledIntensity = points{r}(c).scaledIntensity;
                        pointsPerCell{i} = pointsPerCell{i} + size(pointsTemp{r}(c).channels, 1);
                    case '2d' % separate points based on zslice to analyze separately
                        idx = find(points{r}(c).channels(:,3) == i);
                        pointsTemp{r}(c).channels = points{r}(c).channels(idx,:);
                        pointsTemp{r}(c).intensity = points{r}(c).intensity(idx);
                        pointsPerCell{i} = pointsPerCell{i} + size(pointsTemp{r}(c).channels, 1);
                        
                end
                % add for filler
                %pointsTemp{r}(c).intensity = zeros(length(pointsTemp{r}(c).channels),1);
            end
        end
        
        pointsxcell{i} = pointsTemp;
    end



end