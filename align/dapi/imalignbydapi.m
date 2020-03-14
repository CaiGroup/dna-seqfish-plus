function [moving, fixed, movingRaw] = imalignbydapi(imagePath, movingImageString, fixedImageString, varargin)
% imalignbydapi aligns the images with dapi and returns the original image
% and aligned image.
%
% Assumptions:
% 1. dapi is the last channle of each image
% 2. images are in the same path
%
% optional Variable channelArray is used to grab specific channels of
% images: ex, user would like channel 1 from the moving (first) image and
% channel 3 from the fixed (second) image 
% >> channelArray = [1 3];
% >> [moving, fixed] = imalignbydapi(imagePath, movingImageString, fixedImageString, channelArray)
% otherwise function will return all the channels of moving and fixed
% images
%
% Author: Nico Pierson
% Date 4/11/2019


    %% Set up optional Parameters
    numvarargs = length(varargin);
    if numvarargs > 1
        error('myfuns:imalignbydapi:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {[]}; % default of using 7 x 7 pixel grid for gaussian function
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Place optional args in memorable variable names
    [channelArray] = optargs{:};
    
    
    
    %% Get the variables of the image
    [channelNumber, zsliceNumber] = getimageinfo(imagePath, fixedImageString);
    numberOfHybChannels = channelNumber - 1;
    
    
    %% Get the images
    movingDapi = grabim(imagePath, movingImageString, channelNumber);
    fixedDapi = grabim(imagePath, fixedImageString, channelNumber);
    
    
    if ~isempty(channelArray)
        numberOfHybChannels = 1;
        movingImage{1} = grabim(imagePath, movingImageString, channelArray(1));
        fixedImage{1} = grabim(imagePath, fixedImageString, channelArray(2));
    else
        movingImage = cell(1, numberOfHybChannels);
        fixedImage = cell(1, numberOfHybChannels);
        for i = 1:numberOfHybChannels
            movingImage{i} = grabim(imagePath, movingImageString, i);
            fixedImage{i} = grabim(imagePath, fixedImageString, i);
        end
    end
    
    
    %% Tform
    tformDapi = grabtform(movingDapi, fixedDapi);
    fprintf('Dapi Tform is: \n');
    [tformDapi.T] % print tform
    
    % reference image
    if zsliceNumber >= 16
        outputRef = imref3d(size(movingImage{1}));
        outputRefDapi = imref3d(size(fixedDapi));
    else
        outputRef = imref2d(size(movingImage{1}(:,:,1)));
        outputRefDapi = imref2d(size(fixedDapi));
    end
    
    %% Only get images
    if length(movingImage) == 1 % if only one hyb image
        moving = imwarp(movingImage{1}, tformDapi, 'OutputView', outputRef);
        fixed = fixedImage{1};
        movingRaw = movingImage{1};
    else % else store in cell array
        moving = movingImage;
        for j = 1:numberOfHybChannels
            moving{j} = imwarp(moving{j}, tformDapi, 'OutputView', outputRef);
        end
        movingRaw = movingImage;
        fixed = fixedImage;
    end
    
    %% show the alignment image
    % use function for saving in imagej
    dapi = cell(1, 2);
    movingDapiTform = imwarp(movingDapi, tformDapi, 'OutputView', outputRefDapi);
    dapi{1} = fixedDapi;
    dapi{2} = movingDapiTform;
    pass = saveimagej(dapi, fullfile(pwd, 'dapiAlignmentCheck.tif'));
    

end