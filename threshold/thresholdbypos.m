function threshold = thresholdbypos(experimentDir, experimentName, position, ...
    numRounds, numPseudoChannels, typedots)
% function threshold allows the user to assign threshold to all
% hybridizations. Choose the position index and channel index.
%
% Assumptions:
% 1. Image names assumed to be in format: 'MMStack_Pos{number}.ome.tif'.
%
% 2. Folder names assumed to be in format: 'HybCycle_{number}'.
%
% folderCompare gets the folders with this string, strCompare gets the
% filename
%
% Dependencies: AlignImages package ?
%
% Author: Nico Pierson
% Date: 2/20/2019
% nicogpt@caltech.edu
% Modified: 7/29/2019

    % Set up paths in Matlab: CompareDotsError and Fiji packages
    %checkcomparedotserrorpath();
    % Need to check path for bfmatlab path
    checkbfmatlabpath();
    checkfijipath();


    % Set up variables
    options.typedetectdots = typedots; % 'exons' or 'introns'
    options.processimages = false; % do it before in this function
    HIGH_THRESHOLD_VALUE = 99999;
    images = cell(numRounds, numPseudoChannels);
    threshold = ones(numRounds, numPseudoChannels) * HIGH_THRESHOLD_VALUE;
    barcodeRange = 1:numRounds;
    % get number of channles, and z
    imageInitPath = fullfile(experimentDir, 'HybCycle_0', 'MMStack_Pos0.ome.tif');
    positionInit = 0;
    [~, sizeC, sizeZ, ~, ~] = grabimseries(imageInitPath, positionInit);
    numFolderPerBarcode = numPseudoChannels / (sizeC - 1);

    try
        fprintf('Preprocessing Threshold...\n');
        for barcode = barcodeRange
            fprintf('Barcode Round: %.0f ', barcode);
            folderStart = (barcode - 1) * numFolderPerBarcode + 1; % Starts at 10 and ends at 15 for 16 serialhybs            
            folderEnd = folderStart + numFolderPerBarcode - 1;
            folderRange = folderStart:1:folderEnd;
            pseudoIterator = 1;

            fprintf('HybCycle: ');
            for folder = 1:length(folderRange)
                fprintf(' %.0f', folderRange(folder));
                folderName = ['HybCycle_' num2str(folderRange(folder))]; 
                imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
                imagePath = fullfile(experimentDir, folderName, imageName);

                % get the image in the specified channel
                uniqueString = 'imageTemp-barao3290dsfj';
                [image, ~, ~, ~, ~] = grabimseries(imagePath, position);
                for c = 1:(sizeC-1)
                    images{barcode, pseudoIterator} = imagejbackgroundsubtraction(image{c}, uniqueString, experimentDir);
                    pseudoIterator = pseudoIterator + 1;
                end
            end
            fprintf('\n');


        end
        fprintf('Choose Threshold:\n');
        for barcode = barcodeRange
            fprintf('barcode: %.0f\n', barcode);
            for pseudoch = 1:numPseudoChannels
                fprintf('pseuod channel: %.0f, ', pseudoch); 
                options.channelLabel = pseudoch;
                options.folderLabel = barcode;
                threshold(barcode, pseudoch) = getthreshold(images{barcode, pseudoch}, options); % soemething is wrong here
            end
            fprintf('\n');
        end
        
        % Save in the threshold directory
        thresholdDir = fullfile(experimentDir, 'threshold');
        if ~exist(thresholdDir, 'dir')
           mkdir(thresholdDir); 
        end
        savePath = fullfile(thresholdDir, ['threshold-pos' num2str(position) '-R' num2str(numRounds) 'C-' num2str(numPseudoChannels) '-' experimentName '.mat']);
        save(savePath, 'threshold');
    
    catch err
        rethrow(err);
    end

end