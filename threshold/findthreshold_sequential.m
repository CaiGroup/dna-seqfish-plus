function threshold = findthreshold_sequential(folderArray, hybnum, position) % calculate the threshold for each position???
% returns threshold values for each set of images using the folder Array of
% the set and the processed images hybnum, and the position number. The
% final set of threshold values are a hyb number by channel matrix of
% double.
%
% Author: Nico Pierson
% nicogpt@caltech.edu
% Date: 1/22/2019
% Modified: 1/24/2019

    % Initialize Date
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    dateSaveString = datestr(dateStart, formatDate);

    %% Declare Variables
    numOfFolders = length(folderArray);
    numOfChannels = size(hybnum(1).color, 2);
    regMax = cell(numOfFolders, numOfChannels);
    logFish = cell(numOfFolders, numOfChannels);
    threshold = zeros(numOfFolders, numOfChannels);

    %% Preprocess Images Before Finding Threshold
    fprintf('Preprocessing Images: position %.0f\n', position);
    tic
    for folder = 1:numOfFolders
        for channel = 1:numOfChannels
            [regMax{folder, channel}, logFish{folder, channel}] = preprocessdots(hybnum(folder).color{channel});% folder + 1 for this set....make variable that has the starting folder????
        end
    end
    toc

    %% Find the threshold for each Set of Images
    disp('Finding Threshold');
    for folder = 1:numOfFolders 
        for channel = 1:numOfChannels 
           threshold(folder,channel) = manualthreshold(folder, channel, ...
               hybnum(folder).color{channel}, regMax{folder, channel}, logFish{folder, channel}); 
        end
    end


    save(['threshold_Pos' num2str(position) '_' dateSaveString], 'threshold');

end