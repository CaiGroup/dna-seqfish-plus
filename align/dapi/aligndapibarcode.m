function pass = aligndapibarcode(expDirectory, positionArray, numberOfBarcodes, numberOfSerialHybs)
% aligndapibarcode aligns and saves the images by aligning the dapi. In
% addition the function saves the dapi registration for each barcode round
% and the first dapi images for each barcode round as allHybRegistration
%
% inputs: positionArray, numberOfBarcodes, numberOfSerialHybs, experimentdirecotry, 
%
% 1. Name of images will be assumed to be 'MMStackPos{number}.tif'.
% 2. expDirectory will assume to store images in folders 'HybCycle_{number}'
% and the number starts from 0.
% 3. saveDirectory will be 'organizehybs' in imageDirectory.
% 4. Position is an array, ex: [0, 2, 3] for aligning positions 0, 2 and 3.
% 4. Dapi is assumed to be the last channel.
% 
% Testing:
% 4/10/2019
% So far function only tested on '2019-04-04-E14-DNA-FISH-full-chrom'
% experiment performed by Yodai; details of his experiment: 5 barcode
% rounds, 16 serial hybridizations, 1 position
%
% Things to do: check the auto threshold if it works, the threshold 10000
% returns too many points

    % make sure bfmatlab and Fiji.app is in the path
    % add path of CompareDotsError, and AlignImages
    % Need to make a function for this
    %addpath('C:\github\streamline-seqFISH\src\AlignImages', '-end');
    %addpath(fullfile('C:', 'github', 'streamline-seqFISH', 'src', ...
        %'CompareDotsError'), '-end');
    
    % Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    dateSaveString = datestr(dateStart, formatDate);
    dateStartString = datestr(dateStart);

    % set up directories and filename
    tic
    pass = 0;
    imageTestPath = fullfile(expDirectory, 'HybCycle_0');
    stringCompare = 'MMStack_Pos0.ome.tif';
    [~, numDapi, zSlice, ~, ~] = grabimseries(fullfile(imageTestPath, stringCompare), 0);
    %[numberOfChannels, zSlice] = getimageinfo(imageTestPath, stringCompare);
    numHybChannels = numDapi - 1;
    %dapiChannel = numberOfChannels;
    positionRange = positionArray;
    barcodeRange = 1:numberOfBarcodes;
    basePath = expDirectory;
    folderPath = fullfile(basePath, 'organizehybs');
    if ~exist(folderPath, 'dir')
        mkdir(folderPath);
    end
    tformDapi = [];
    positionCount = 1;
    allHybRegImage = cell(1, length(barcodeRange));

    for position = positionRange

        saveDir = fullfile(folderPath, ['pos' num2str(position)]);
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        dapiref = [];
        for barcode = barcodeRange
            fprintf('Transformation for barcode %.0f...\n', barcode);
            folderStart = (barcode - 1) * numberOfSerialHybs + 1; % Starts at 17 and ends at 33
            folderEnd = folderStart + numberOfSerialHybs - 1;
            folderRange = folderStart:1:folderEnd;


            dapi = cell(1, length(folderRange));
            for folderIndex = 1:length(folderRange) % variables for storing tformDapi
                tformDapi(barcode).barcode(folderIndex).folder = [];
            end
            for channelIndex = 1:numHybChannels % variables for storing images
                hybIms{channelIndex} = cell(1, length(folderRange));
            end


            for folder = 1:length(folderRange)
                folderName = ['HybCycle_' num2str(folderRange(folder)-1)];
                imageDirectory = fullfile(basePath, folderName);
                imageFileName = ['MMStack_Pos' num2str(position) '.ome.tif'];


                % get the images for each channel
                [allIms, numDapi, zSlice, ~, ~] = grabimseries(fullfile(imageDirectory, imageFileName), position);
                for c = 1:numDapi
                    if c < numDapi
                        hybIms{c}{folder} = allIms{c};
                    else
                        dapi{folder} = allIms{c};
                    end
                end
                

                % Apply background nonuniform corrections - maybe add as
                % option
                %image(folder).channel{1} = uint16(double(image1) ./ backcorrections{position+1, 1});

                %% Get the tform and apply the transformation
                % Use the dapi transformations and the chromatic aberrations (from barcoded
                % experiments)
                fprintf('Apply Dapi tform for position %.0f and folder %.0f...\n', position, folder);
                if folder ~= 1 || barcode ~= 1 % if not the first folder
                    tformTemp = grabtform(dapi{folder}, dapiref);

                    if zSlice < 16
                        tformMat = zeros(3, 3, size(zSlice,3));
                        for z = 1:zSlice
                            tformMat(:,:,z) = tformTemp{z}.T;
                        end
                        clearvars tformTemp
                    end

                    % get median of the tform
                    if zSlice < 16
                        tformDapi(barcode).barcode(folder).folder = affine2d(median(tformMat, 3));
                        for c = 1:numHybChannels
                            hybIms{c}{folder} = imwarp(hybIms{c}{folder}, tformDapi(barcode).barcode(folder).folder, 'OutputView', imref2d(size(hybIms{c}{folder})));
                        end
                        dapi{folder} = imwarp(dapi{folder}, tformDapi(barcode).barcode(folder).folder, 'OutputView', imref2d(size(dapi{folder})));
                    else
                        tformDapi(barcode).barcode(folder).folder = tformTemp;
                        for c = 1:numHybChannels
                            hybIms{c}{folder} = imwarp(hybIms{c}{folder}, tformDapi(barcode).barcode(folder).folder, 'OutputView', imref3d(size(hybIms{c}{folder})));
                        end
                        dapi{folder} = imwarp(dapi{folder}, tformDapi(barcode).barcode(folder).folder, 'OutputView', imref3d(size(dapi{folder})));
                    end
                else
                    if zSlice < 16
                        tformDapi(barcode).barcode(folder).folder = affine2d(eye(3));
                    else
                        tformDapi(barcode).barcode(folder).folder = affine3d(eye(4));
                    end
                    dapiref = dapi{1};
                end

                if folder == 1 % get the images for allhybregistration for first image in each barcode round
                    allHybRegImage{barcode} = dapi{folder};
                end

            end



            %% save the dapi files for registration check for each barcode round
            Miji(false);
            for f = 1:numberOfSerialHybs
                    namesh{f} = ['C' num2str(f) '-'  num2str(barcode) '.tif'];
                    MIJ.createImage(namesh{f}, dapi{f}, true);
            end

            str = [];
            for f = 1:numberOfSerialHybs
                    str = [str ' image' num2str(f) '=C' num2str(f) '-' num2str(barcode) '.tif'];
            end

            try
                MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
                MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(numberOfSerialHybs) ' slices=' num2str(zSlice) ' frames=1 display=Grayscale']);
                savePath = fullfile(saveDir, ['dapi_pos' num2str(position) '_bar' num2str(barcode)]);
                MIJ.run('Save', ['save=[' savePath '.tif' ']']);
                MIJ.run('Close All')
            catch
                MIJ.exit;
                error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
            end


            %% Need to save the images per barcode round for each channel
            for c = 1:numHybChannels


                for f = 1:numberOfSerialHybs
                        namesh{f} = ['C' num2str(f) '-'  num2str(barcode) '.tif'];
                        MIJ.createImage(namesh{f}, hybIms{c}{f}, true);
                end

                str = [];
                for f = 1:numberOfSerialHybs
                        str = [str ' image' num2str(f) '=C' num2str(f) '-' num2str(barcode) '.tif'];
                end


                try
                    MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
                    MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(numberOfSerialHybs) ' slices=' num2str(zSlice) ' frames=1 display=Grayscale']);
                    savePath = fullfile(saveDir, ['hyb-pos' num2str(position) '-round' num2str(barcode) '-ch' num2str(c)]);
                    MIJ.run('Save', ['save=[' savePath '.tif' ']']);
                    MIJ.run('Close All')

                catch
                    MIJ.exit;
                    error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
                end

            end
            MIJ.exit;

            save(['imagedata_pos' num2str(position) '_round' num2str(barcode)], 'hybIms');
            clearvars hybIms
        end

        %% Need all hyb registration
        Miji(false);

        for b = 1:length(barcodeRange)
            namesh{b} = ['C' num2str(b) '-'  num2str(1) '.tif'];
            MIJ.createImage(namesh{b}, allHybRegImage{b}, true);
        end

        str = [];
        for b = 1:length(barcodeRange)
            str = [str ' image' num2str(b) '=C' num2str(b) '-' num2str(1) '.tif'];
        end

        try
            MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(length(barcodeRange)) ' slices=' num2str(zSlice) ' frames=1 display=Grayscale']);
            savePath = fullfile(saveDir, 'AllHybRegistrationCheck');
            MIJ.run('Save', ['save=[' savePath '.tif' ']']);
            MIJ.run('Close All')
            MIJ.exit;
        catch
            MIJ.exit;
            error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
        end

        fprintf('Done getting images in position %.0f...\n', position);
        save([saveDir filesep 'yodaiFullChromosome-pos' num2str(position) '.mat'], 'tformDapi');
        positionCount = positionCount + 1; % Iterate position
    end

    pass = 1;
    


    % Save time completed
    dateEnd = datetime;
    dateStartToEnd = [dateStart dateEnd]; 
    dateTotal = diff(dateStartToEnd);
    [h,m,s] = hms(dateTotal);
    dateTotalString = sprintf('Total Elapsed Time: %.0f hours %.0f minutes and %.0f seconds...\n', h, m, s);
    % Save the Data and Time the script started
    fprintf('------------------------------------------------------------------------\n');
    fprintf('------------------------------------------------------------------------\n'); 
    fprintf('Completed Script Successfully\n');
    fprintf(['Date: ' datestr(dateEnd) '\n']);
    fprintf(dateTotalString);
end