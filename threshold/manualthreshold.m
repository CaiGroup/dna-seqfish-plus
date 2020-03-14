function [threshold, points, intensity] = manualthreshold(folder, channel, image, regMax, logFish, typedetectdots, varargin)
% function threshold finds the threshold from images based on user input.
%
% Update 1/24/2019: added title for images, and paused twice to allow user 
% to adjest contrast for images.
%
% Future Additions: 
% 1. Update user interface: as of now only has current threshold and number
% of dots.
% 2. Add a machine learning algorithm to automate threshold and have user
% accept or reject threshold value.
% 3. Add an option to save the dots
%
% Author: Nico Pierson
% Date: 1/22/2019
% nicogpt@caltech.edu
% Modified: 1/24/2019

    numvarargs = length(varargin);
    if numvarargs > 2
        error('src:manualthreshold:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs == 1 
        if ~varargin{1} == 1 && ~varargin{1} == 0
            error('myfun:manualthreshold:WrongInput', ...
                'savedots requires type boolean');
        end
    end

    % set defaults for optional inputs
    optargs = {false, []}; % options are 'introns' or 'exons'
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [saveDots, saveDotsPath] = optargs{:};


    %% Set global variables
    newThreshold = []; % use threshold to set new
    folderString = num2str(folder);
    channelString = num2str(channel);
    figureName = ['Hyb ' folderString ' Channel ' channelString];
    returnBoolean = 0;
    highThresholdValue = 99999;
    rawImage = [];
    m = [];
    changeThreshold = false; % switch to tell when the threshold needs to be adjusted

    %% Continuous While Loop
    while(1)
        if isempty(newThreshold)
            % Enter threshold value if first time
            prompt = {'Enter threshold:'};
            title = ['Input for Hyb ' folderString ' Channel ' channelString];
            dims = [1 35];
            definput = {'5000'};
            thresholdString = inputdlg(prompt,title,dims,definput);
            if isempty(thresholdString)
                promptMessage = sprintf('Do you want to Continue processing the same image,\nor Exit the function, \n or Skip to the next image?');
                button = questdlg(promptMessage, 'Continue', 'Continue', 'Exit', 'Skip', 'Continue');
                if strcmp(button, 'Exit')
                    error('exited function in manualthreshold.m');
                elseif strcmp(button, 'Skip')
                    threshold = highThresholdValue; % Apply high threshold value for no dot detection
                    return;
                else
                    continue;
                end
            else
                threshold = round(str2double(thresholdString{1}));  
            end
        else
            changeThreshold = true;
            threshold = newThreshold;
        end
        
        %% Find Dots for exons or introns
        if strcmp(typedetectdots, 'exons') || strcmp(typedetectdots, 'exons2d')
            % Find Dots for exons
            m = regMax & logFish > threshold;
            
        elseif strcmp(typedetectdots, 'log') 
            if changeThreshold
                logFish = rawImage; % assign logFish as raw Image
            else
                rawImage = logFish; % keep the raw Image
            end
            % Find Dots for exons
            msk = true(3,3,3);
            msk(2,2,2) = false;
            apply = logFish < threshold;
            logFish(apply) = 0;
            s_dil = imdilate(logFish,msk);
            m = logFish > s_dil; 
            
        elseif strcmp(typedetectdots, 'introns')
            if changeThreshold
                logFish = rawImage; % assign logFish as raw Image
            else
                rawImage = logFish; % keep the raw Image
            end
            % Find Dots for introns: doesn't use regMax
            regMask = true(3,3,3);
            regMask(2,2,2) = false;
            % assign, to every voxel, the maximum of its neighbors
            belowThresh = logFish < threshold;
            logFish(belowThresh) = 0; % logFish needs to stay the same
            dilate = imdilate(logFish,regMask);
            % create logical matrix of dotsLogical where 1 is the voxel's value is
            % greater than its neighbors
            m = logFish > dilate;
        end
        % Calculate the dots found
        % Get x, y, z coordinates of dots
        [y,x,z] = ind2sub(size(logFish), find(m == 1));
        points = [x y z]; 
        numberOfDots = size(points,1);
        intensity = zeros(length(y), 1);
        removeDots = [];
        for i = 1:length(y)
            intensity(i,1) = rawImage(y(i), x(i), z(i));
            if intensity(i,1) > 30000 % remove very bright dots - max is 65000
                removeDots = cat(2, removeDots, i);
            end
        end
        points(removeDots,:) = [];
        intensity(removeDots) = [];
 
        

        %% Check image
        % Initial Settings
        initMag = 'fit'; % 20% of the image; use 'fit' if too large; 20 is too large on macbooks
        close all;
        
        imageMaxDisplayAdd = 3000;
        % Image 1: raw image 
        figure1Name = ['Raw Image: ' figureName];
        figure('Name', figure1Name)  % test logFish
        imshow(max(image,[],3),[min(min(max(image,[],3))) mean(mean(max(image,[],3)))+imageMaxDisplayAdd], 'InitialMagnification', initMag);
        uicontrol('Position',[20 20 200 30],'String','Adjust Contrast',...
              'Callback','imcontrast(1)');
        movegui(1, 'northwest');
        
        % Image 2: background corrected
        figure2Name = ['Image: ' figureName];
        figure('Name', figure2Name)  % test logFish
        imshow(max(logFish,[],3),[min(min(max(logFish,[],3))) mean(mean(max(logFish,[],3)))+imageMaxDisplayAdd], 'InitialMagnification', initMag);
        uicontrol('Position',[20 20 200 30],'String','Adjust Contrast',...
              'Callback','imcontrast(2)');
        movegui(2, 'northeast');
        
        % Image 3: image with dots
        figure3Name = ['Raw Image with Dots: ' figureName];
        figure('Name', figure3Name)  % test logFish
        imshow(max(image,[],3),[min(min(max(image,[],3))) mean(mean(max(image,[],3)))+imageMaxDisplayAdd], 'InitialMagnification', initMag);
        uicontrol('Position',[20 20 200 30],'String','Adjust Contrast',...
              'Callback','imcontrast(3)');
        hold on;
        [v2,v1]=find(max(m,[],3)==1);
        scatter(v1(:),v2(:),75);
        hold off;
        movegui(3, 'southwest');
        
        % Image 4: subtracted image with dots
        figure4Name = ['Subtracted Image with Dots: ' figureName];
        figure('Name', figure4Name)  % test logFish
        imshow(max(logFish,[],3),[min(min(max(logFish,[],3))) mean(mean(max(logFish,[],3)))+imageMaxDisplayAdd], 'InitialMagnification', initMag);
        uicontrol('Position',[20 20 200 30],'String','Adjust Contrast',...
              'Callback','imcontrast(4)');
        hold on;
        [v2,v1]=find(max(m,[],3)==1);
        scatter(v1(:),v2(:),75);
        hold off;
        movegui(4, 'southeast');
        LinkFigures(1:4);
        
        % Options to accept threshold or enter new value
        figure4Name = ['Accept or Change Threshold for ' figureName];
        f = figure('Name', figure4Name);
        numberOfDotsString = ['Number of Dots Detected: ' num2str(numberOfDots)];
        thresholdAcceptString = ['Accept Current Threshold of ' num2str(threshold)];
        uicontrol('Style', 'text', 'Position', [15 125 110 250], ...
            'String', numberOfDotsString, 'FontSize', 12);
        uicontrol('Style', 'pushbutton', 'Position', [135 290 300 100], ...
            'String', thresholdAcceptString, 'FontSize', 11, ...
              'Callback', @threshreturn);
        uicontrol('Style', 'pushbutton', 'Position', [135 140 300 100], ...
            'String', 'Change to New Threshold Value', 'FontSize', 11, ...
              'Callback', @threshedit);
        uiNewThreshold = uicontrol('Style', 'edit', 'Position',[190 40 200 50], ...
            'String', '5000', 'FontSize', 12);
        movegui(4, 'southeast');
        % wait until the handler or window is closed
        waitfor(f);
        close all;
        if returnBoolean % return if accept threshold
            fprintf('*Selected Threshold is %.0f*\n', threshold);
            if saveDots
                save(saveDotsPath, 'points', 'threshold');
            end
            return;
        end
    end
    
    function threshreturn(source, event)
        returnBoolean = true;
        close(gcf) % close handler and window
    end
    
    function threshedit(source, event)
        newThresholdString = uiNewThreshold.String;
        newThreshold = str2double(newThresholdString); % check if value is a number
        if ~isscalar(newThreshold) || ~isnumeric(newThreshold)
            error('Threshold Value is not Scalar');
        end
        close(gcf)
    end
    
end