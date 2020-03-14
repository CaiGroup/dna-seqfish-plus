function [threshold,points,intensity] = getthreshold(image, options)
% returns the treshold from an image. 
%
%
% Dependencies: 
% - Fiji.app/scripts directory
% - mij.jar: function asks for location for ImageJ
% - FindThreshold package
%
% Example of setting up struct for options:
%    >> options.typedetectdots = 'exons';
%    >> options.savefigs = true;
%    >> options.threshold1 = 800; % set numeric threshold value image1
%    >> options.threshold2 = 'choose'; % default already: choose threshold
%    value for image 2
%    >> options.processimages = true;
%
% Author: Nico Pierson
% nicogpt@caltech.edu
% Date: 1/28/2019
% Modified: 2/12/2019
% Options, try and catch code adapted from 'loadtiff.m' by YoonOh Tak 2012

    % use errcode to output different types of errors
    errcode = -1;
    try
    %% Check options structure for extra data
        if isreal(image) == false
            errcode = 1; assert(false);
        end
        if nargin <= 1
            options.savefigs = false;
            options.processimages = true;
            options.debugimages = false;
            options.threshold1 = 'choose';
            options.threshold2 = 'choose';
            options.typedetectdots = 'exons'; % 'introns' or 'exons'
            options.colocalizationradius = 3; % default 3; colocalize dots in this pixel radius
            options.intarearadius = 3; % default is 3; compare intensity areas in this pixel radius
            options.usefiji = false; % option to use fiji with ImageJ
            options.savedots = false; % option to save dots
            options.pathsavedots = []; % path to save dots
            options.sliding = false;
            options.backradius = 3;
            options.channelLabel = 0;
            options.folderLabel = 0;
        end
        if isfield(options, 'savefigs') == 0
            options.savefigs = false;
        end
        if isfield(options, 'processimages') == 0
            options.processimages = true;
        end
        if isfield(options, 'debugimages') == 0
            options.debugimages = false;
        end
        if isfield(options, 'threshold') == 0
            options.threshold = 'choose';
        end
        if isfield(options, 'typedetectdots') == 0
            options.typedetectdots = 'exons';
        end
        if isfield(options, 'usefiji') == 0
            options.usefiji = false;
        end
        if isfield(options, 'savedots') == 0
            options.savedots = false;
        end
        if isfield(options, 'pathsavedots') == 0
            options.pathsavedots = [];
        end
        if isfield(options, 'sliding') == 0
            options.sliding = false;
        end
        if isfield(options, 'backradius') == 0
            options.backradius = 3;
        end
        if isfield(options, 'channelLabel') == 0
            options.channelLabel = 0;
        end
        if isfield(options, 'folderLabel') == 0
            options.folderLabel = 0;
        end
        
        
        %% Start finding the localization and intensity error
        % Get date to save files for unique string
        %dateStart = datetime;
        %formatDate = 'yyyy-mm-dd';
        %dateSaveString = datestr(dateStart, formatDate);
        %fprintf('Started thresholdbypos.m on %s\n', dateSaveString);

        
        %% Option to Process Images
        % Use ImageJ Background Selection
        if options.processimages

            checkfijipath();
            uniqueString = 'imageTemp-09238409238409';
            % Subtract the background
            image1Process = imagejbackgroundsubtraction(image, ...
                uniqueString, options.backradius, options.sliding);
        end

        % Need to keep processed images for dot finding and original image
        % for intensities: no problems for intron finding dots, but it is a
        % problem for finding exons part

        %% Option to Find the threshold for each Set of Images
        thresholdImage = image;
        if strcmp(options.threshold, 'choose')
            %disp('Choose a threshold for first Image');

            channel = options.channelLabel;
            folder = options.folderLabel;
            %disp('Finding Threshold');
            % Process Image First with regional maxima and a lap filter for
            % 'exons' type of dot detection
            if options.processimages
                thresholdImage = image1Process;
            end
            [regMax, logFish] = preprocessdots(thresholdImage, options.typedetectdots);
            % Get first threshold
            [threshold,points,intensity] = manualthreshold(folder, channel, ... % set channel and folder to 0 for null
               image, regMax, logFish, options.typedetectdots, options.savedots, options.pathsavedots);
           %disp('Completed Threshold for Image 2');
        else
            if isnumeric(options.threshold)
                threshold = options.threshold;
            else
                errcode = 3; assert(false);
            end
            [points, intensity, ~, ~] = detectdotsv2(thresholdImage, threshold, options.typedetectdots);
        end
        
        
    catch exception
        %% catch the exceptions in the function
        % Update the error messages for errors
        switch errcode
            case 0
                error 'Invalide path.';
            case 1
                error 'It does not support complex numbers.';
            case 2
                error '''data'' is empty.';
            case 3
                error 'Threshold option is invalid. Usage: enter numeric value or "choose" for options.threshold';
            otherwise
                rethrow(exception);
        end
    end
end