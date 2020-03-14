function [movingAligned, fixed] = imerrorpoints(imagePath, movingImageString, fixedImageString, channelArray)
% function for Yodai to compare two images, align the dapi images, compare
% the error in localization, and use the colocalized points to transform
% the images and view the error betwen the starting points and the
% transformed points.
%
% It is assumed the two images are in the same path with unique strings
% Variables:
% channelArray is used to grab specific channels of images: ex, user would 
% like channel 1 from the moving (first) image and channel 3 from the fixed 
% (second) image:
% >> channelArray = [1 3];
% >> pass = imerrorpoints(imagePath, movingImageString, fixedImageString,
% channelArray);
%
% Dependencies: packages CompareDotsError and FindThreshold
%
% Author: Nico Pierson
% Date: 4/12/2019
    

    %% Align the images using dapi
    [moving, fixed, movingRaw] = imalignbydapi(imagePath, movingImageString, fixedImageString, channelArray);



    %% Colocalize the points 
    % set up options to return points and threshold values
    options.return = true; % should be true
    options.typedetectdots = 'exons'; % toggle different type of dot detection: 'exons' or 'introns'
    options.colocalizationradius = 1.76; % 1.76 is the default value
    options.intarearadius = 3; % if you want to change the calculated area
    options.printoutput = true; % toggle the command window output
    options.processimages = true; % images will be background subtracted before searching for points
    %options.threshold1 = 30000; % use if you already have threshold values
    %options.threshold2 = 32000;
    [colocalizationData, movingPoints, fixedPoints, movingThreshold, fixedThreshold] = imerror(moving, fixed, options);
    


    %% Compare the error between dapi aligned images and point aligned images
    [movingAligned, ~, tformPoints] = imalignhybpoints(moving, fixed, movingRaw, movingPoints, fixedPoints, movingThreshold);

end