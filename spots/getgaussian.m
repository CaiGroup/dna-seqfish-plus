function [gaussPoints, gaussInt] = getgaussian(points, image, varargin)
% getgaussian gets points and returns gaussian fit points
%
% Author: Nico Pierson
% Date: 4/4/2019
% Email: nicogpt@caltech.edu
% Modified:
    
    %% Set up optional Parameters
    numvarargs = length(varargin);
    if numvarargs > 1
        error('myfuns:gaussPoints:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {3}; % default of using 7 x 7 pixel grid for gaussian function
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Place optional args in memorable variable names
    [pixelRadius] = optargs{:};

    %% Adjust Bead Coordinates in Dots Image with 2D Gaussian Fit
    % Initialize Variables
    options = optimset('Display','off');
    numberOfCoordinates = size(points,1);
    numberOfGaussianParameters = 7;
    gaussData = zeros(numberOfCoordinates, numberOfGaussianParameters);
    zSliceData = zeros(numberOfCoordinates, 1);
    
    for dotID = 1:numberOfCoordinates % For each dot location 
        zSlice = points(dotID, 3);
        dataRangeX = round((-pixelRadius:pixelRadius) + points(dotID,1)); % parenthesis are necessary
        dataRangeY = round((-pixelRadius:pixelRadius) + points(dotID,2)); % Need to round before inserting into image
        data = double(image(dataRangeY, dataRangeX, zSlice)); % select region in a 3x3 pixel radius % error NOT IN the IMAGE
        x0 = [min(min(data)) max(max(data))-min(min(data)) pixelRadius+1 pixelRadius+1 1 1 45]; % x0  = [min max-min(range) 4 4 1 1 45] is the starting point
        dataSize = size(data); % 7 x 7 array
        f = @(x)gauss2dfrobenium(x,dataSize,data); % Cost function for 2D Gaussian image
        % find minimum frobenium norm of the 2d Gaussian adjusted pixels
        % and return the values for min max-min x0 y0 stdx stdy and theta
        [minCoords, fva] = fminsearch(f,x0,options); 
        gaussData(dotID,:) = minCoords; % Assign value of gaussian fitted coordinates
        % Assign zSlice
        zSliceData(dotID) = zSlice;

    end
    %{
    %% filter out hot pixels...by accepting only points with a stdy of 1 or greater
    removeInd = find(gaussData(:,6) < 1);
    gaussData(removeInd,:) = [];
    points(removeInd,:) = [];
    zSliceData(removeInd) = [];
    %}
    % Take the 3rd and 4th columns 
    gaussPoints = gaussData(:,3:4); 
    % Adjust the gaussian coordinates by adding to original coordinates and
    % subtracting pixel width
    gaussPoints(:,1:2) = points(:,1:2) + gaussData(:,[3 4]) - [pixelRadius, pixelRadius]-1; % keep gaussData(:,[3 4]) for x and y
    % Add zSlices to col 3
    gaussPoints(:,3) = zSliceData;
    gaussInt = gaussData(:,1) + gaussData(:,2);
    %{
    %% add intensity and filter
    
    removeInt = find(gaussInt < 800);
    gaussPoints(removeInt,:) = [];
    gaussInt(removeInt) = [];
%}
end
