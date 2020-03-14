function [I1, I2, D1, D2] = grabdapi(image1, image2, varargin)
% grabdapi returns the image1, image2 with the hybs and dapi as separate
% files, assuming dapi is the last channel in both images.
%
% Inputs: image path and part of the unique string of image
%
% Optional inputs: index of channel (starting at 1; can use vector [1
% 3] to grab channels 1 and 3) or leave blank to get all channels, zslice
% index can be a vector [3:8] to grab zslices from 3 to 8 or [4 6 9] to
% grab only 4, 6, and 9 zslices.
%
%
% To Do List:
% 1. NEED DOCUMENTATION AND EXAMPLES.
% 2. Debug options if dapi channel is included.
%
% Optional Parameters: imageOrg, strCompare, zSliceIndex
% imageOrg is an option to organize images from two options: 'yodia intron
% rep4' from xyczt order, and 'linus 10k' from xyzct order
%
%
% zSliceIndex: is an option to get one specific zSlice in a stack.
%
% Outputs: transformed images
%
% Author: Nico Pierson
% Date: 2/26/2019
% Modified: 

    %% Set up optional Parameters for z-slice index
    numvarargs = length(varargin);
    if numvarargs > 4
        error('myfun:alignim:TooManyInputs', ...
            'requires at most 3 optional inputs');
    end

    % Error for type of arguments: zslice index, default are all the
    % zslices, while inputing a vector only gets the zslices wanted
    if numvarargs >= 3
        if ischar(varargin{3})
            if ~strcmp(varargin{2}, 'all')
                error('myfun:alignim:WrongInput', ...
                    'z-slice index requires type string: ''all'' or an index of zslices, ex: 1 or 3:8');
            end
        elseif ~isnumeric(varargin{3}) && ~(mod(varargin{3}, 1) == 0)
            error('myfun:alignim:WrongInput', ...
                'z-slice index is not numeric or an integer: requires type int');
        end
    end

    %% Check if bfmatlab is in the path
    try
        
        % Check if Fiji.app is in the MATLAB path
        foundbfPath = false;
        bfmatlabFolder = 'bfmatlab';
        pathCell = regexp(path, pathsep, 'split');
        for i = 1:size(pathCell, 2)
            pathCell2 = regexp(pathCell{i}, filesep, 'split');
            if strcmp(pathCell2{end}, bfmatlabFolder)
                foundbfPath = true;
                fprintf('bfmatlab directory already in path\n');
                break;
            end
        end
       
        % if not on path, add to MATLAB path
        if ~foundbfPath
            bfmatlabDirectory = getdirectory(bfmatlabFolder);
            addpath(bfmatlabDirectory);
        end
    catch
        error(['bfmatlab package not found or unable to add to path' ...
            'Download package at https://downloads.openmicroscopy.org/bio-formats/4.4.9/']);
        % later add download function:
        % https://downloads.openmicroscopy.org/bio-formats/4.4.9/artifacts/bfmatlab.zip
        % or /loci_tools.jar  file
    end

    
    
    %% Check file existance
    if ~exist(image1, 'file') && ~exist(image2, 'file')
        error 'image variables do not exist';
    end
    
    % Use bfopen to get images with OMEmetadata
    data1 = bfopen(image1);
    data2 = bfopen(image2);
    % get metadata for fileorder
    omeMeta1 = data1{1, 4};
    omeMeta2 = data2{1, 4};
    dimensionOrder1 = convertCharsToStrings(omeMeta1.getPixelsDimensionOrder(0).getValue().toCharArray);
    dimensionOrder2 = convertCharsToStrings(omeMeta2.getPixelsDimensionOrder(0).getValue().toCharArray);
    % number of channels
    numChannels1 = omeMeta1.getChannelCount(0);
    numChannels2 = omeMeta1.getChannelCount(0);
    fprintf('Image1 has %.0f channels\n', numChannels1);
    fprintf('Image2 has %.0f channels\n', numChannels2);
    % number of z-slices
    zSlices1 = omeMeta1.getPixelsSizeZ(0).getValue();
    zSlices2 = omeMeta2.getPixelsSizeZ(0).getValue();
    fprintf('Image1 has %.0f zslices\n', zSlices1);
    fprintf('Image2 has %.0f zslices\n', zSlices2);
    % number of series
    seriesCount = size(data1, 1);
    series1 = data1{1, 1};
    series2 = data2{1, 1};
    metadataList = data1{1, 2};
    % number of planes in series1
    series1_planeCount = size(series1, 1);
    
    
    
    %% Set Defaults for optional inputs
    optargs = {numChannels1, 1:(numChannels1-1), 1:zSlices1, false};
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of dapi channel is the last one, default of index
    % channel are to align all the channels
    [dapiChannel, indexOfChannel, zSliceIndices, viewImage] = optargs{:};

    
    %% Set flag for all channels or just 1
    if ~isnumeric(indexOfChannel) && ~(mod(indexOfChannel, 1) == 0) 
        error('myfun:alignim:WrongInput', ...
            'index of channel is not numeric or an integer: requires type int');
    elseif isscalar(indexOfChannel) && ~(mod(indexOfChannel, 1) == 0)
        % if single channel and integer
    end
    
    
    
    %% Find the directory to grab tif images
    %warning ('off','all');
    % the channel attribute gives the channels of the tiff stack to be
    % loaded; for 2 z-slice 4 channel image ex: index 1 = [1 5]; index 2 =
    % [2 6], index 3 = [3 7] and index 4 = [4 8]
    
    indicesOfDirectory1 = []; % initialize
    indicesOfDapi = [];
    fprintf('Grabbing Indices for Channels');
    fprintf('% .0f', indexOfChannel);
    fprintf('\n');
    fprintf('Grabbing Indices for zSlices');
    fprintf('% .0f', zSliceIndices);
    fprintf('\n');
    switch dimensionOrder1 % Switch to find the Directory of the Tif images based on organization of images
        case 'XYCZT' % if order is channels then zslices, ex [4 8] are dapi in 4 channels and 2 zslices
            fprintf('Image order: XYCZT...\n');
            finalIndex = numChannels1 * zSlices1; % only get the channels for the images
            directoryRange = 1:numChannels1:finalIndex;
            for channel = indexOfChannel % add indices to each channel
                channelIndices = directoryRange + (channel-1); % zslices skip per channel
                % grab only zSlices needed
                channelIndices = channelIndices(zSliceIndices);
                % add to main indices
                indicesOfDirectory1 = cat(2, indicesOfDirectory1, channelIndices); 
            end
            
            % Get the dapi indices
            for channel = dapiChannel % add indices to each channel
                channelIndices = directoryRange + (channel-1); % zslices skip per channel
                % grab only zSlices needed
                channelIndices = channelIndices(zSliceIndices);
                % add to main indices
                indicesOfDapi = cat(2, indicesOfDapi, channelIndices); 
            end
            
            

        case 'XYZCT' % if zslices then channels, ex [7 8] are dapi in 4 channels and 2 zslices
            fprintf('Image order: XYZCT...\n');
            for channel = indexOfChannel
                channelIndices = (1:zSlices1) + (channel-1) * zSlices1; % all zslices together
                % grab only zSlices needed
                channelIndices = channelIndices(zSliceIndices);
                % add to main indices
                indicesOfDirectory1 = cat(2, indicesOfDirectory1, channelIndices);   
            end
            
            % Get the dapi Indices
            for channel = dapiChannel
                channelIndices = (1:zSlices1) + (channel-1) * zSlices1; % all zslices together
                % grab only zSlices needed
                channelIndices = channelIndices(zSliceIndices);
                % add to main indices
                indicesOfDapi = cat(2, indicesOfDapi, channelIndices);   
            end
            
        otherwise
             error 'Invalid dimension order for loaded image';
    end
   
    
    
    %% Grab the Images from the Tif Stack
    % Get all images in field of view range, otherwise get only images in
    % field of range.
    I1 = [];
    I2 = [];
    D1 = [];
    D2 = [];
    % Need to get all the images and the images for dapi...in the end the
    % images will be transformed and use imerror to get the colocalization
    % of the dots
    for index = indicesOfDirectory1 % indices of Directory should be the same for the first and second image
        % Get images from tiff stack
        L1 = series1{index, 1};
        L2 = series2{index, 1};
        % position is not used if there is a string to match to
        if index > 1
            I1 = cat(3, I1, L1); % Catenate L to stack of I
            I2 = cat(3, I2, L2);
        else
            I1 = L1;
            I2 = L2;
        end
    end
    
    for dapi = indicesOfDapi % indices of Directory should be the same for the first and second image
        % Get images from tiff stack
        LD1 = series1{dapi, 1};
        LD2 = series2{dapi, 1};
        % position is not used if there is a string to match to
        if dapi > 1
            D1 = cat(3, D1, LD1); % Catenate L to stack of I
            D2 = cat(3, D2, LD2); % Catenate L to stack of I
        else
            D1 = LD1;
            D2 = LD2;
        end
    end
    fprintf('Successfully grabbed all image planes...\n\n');
    
    %% debug to view images
    if viewImage
        % get the max z-projection and view the image
        close all;
        addHigherPixelRange = 1500;
        figure;
        imshow(max(I, [], 3), 'DisplayRange', [min(min(max(I,[],3))) mean(mean(max(I,[],3))) + addHigherPixelRange], 'InitialMagnification', 'fit');
        pause;
        close all;
    end
    
    %warning ('on','all');
    %Return the images
    
end