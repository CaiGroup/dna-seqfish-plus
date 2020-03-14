function threshold = manualthreshold(folder, channel, image, regMax, logFish)
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
%
% Author: Nico Pierson
% Date: 1/22/2019
% nicogpt@caltech.edu
% Modified: 1/24/2019

    % Set global variables
    newThreshold = []; % use threshold to set new
    folderString = num2str(folder);
    channelString = num2str(channel);
    figureName = ['Hyb ' folderString ' Channel ' channelString];
    returnBoolean = 0;
    highThresholdValue = 99999;
    numberOfDots = []; % Set to null at beginning of function
    m = [];

    % Continuous While Loop
    while(1)
        if isempty(newThreshold)
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
            threshold = newThreshold;
        end
        
        %% Find Dots
        m = regMax & logFish > threshold;
        dotsFound = find(m == 1);
        numberOfDots = size(dotsFound, 1);

        %% Check image
        initMag = 20; % 20% of the image; use 'fit' if too large
        close all;
        figure1Name = ['Original Image: ' figureName];
        figure('Name', figure1Name) 
        imshow(max(image,[],3),[min(min(max(image,[],3))) mean(mean(max(image,[],3)))+5000], 'InitialMagnification', initMag);
        uicontrol('Position',[20 20 200 30],'String','Adjust Contrast',...
              'Callback','imcontrast(1)');
        movegui(1, 'northwest');
        figure2Name = ['Original Image with Dots: ' figureName];
        figure('Name', figure2Name)  
        imshow(max(image,[],3),[min(min(max(image,[],3))) mean(mean(max(image,[],3)))+5000], 'InitialMagnification', initMag);
        uicontrol('Position',[20 20 200 30],'String','Adjust Contrast',...
              'Callback','imcontrast(2)');
        hold on;
        [v2,v1]=find(max(m,[],3)==1);
        scatter(v1(:),v2(:),75);
        hold off;
        movegui(2, 'northeast');
        figure3Name = ['Log Image: ' figureName];
        figure('Name', figure3Name)  % test logFish
        imshow(max(logFish,[],3),[min(min(max(logFish,[],3))) mean(mean(max(logFish,[],3)))+5000], 'InitialMagnification', initMag);
        uicontrol('Position',[20 20 200 30],'String','Adjust Contrast',...
              'Callback','imcontrast(3)');
        movegui(3, 'southwest');
        figure4Name = ['Log Image with Dots: ' figureName];
        figure('Name', figure4Name)  % test logFish
        imshow(max(logFish,[],3),[min(min(max(logFish,[],3))) mean(mean(max(logFish,[],3)))+5000], 'InitialMagnification', initMag);
        uicontrol('Position',[20 20 200 30],'String','Adjust Contrast',...
              'Callback','imcontrast(4)');
        hold on;
        [v2,v1]=find(max(m,[],3)==1);
        scatter(v1(:),v2(:),75);
        hold off;
        movegui(4, 'southeast');
        
        % Options to accept threshold or enter new value
        figure5Name = ['Accept or Change Threshold for ' figureName];
        f = figure('Name', figure5Name);
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
        movegui(5, 'center');
        % wait until the handler or window is closed
        waitfor(f);
        close all;
        if returnBoolean % return if accept threshold
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