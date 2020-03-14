function [hybIms, dapiIms, tform] = aligndapianderson(basePath, folderRange, position, saveDir, splitFactorHybImages, initialRadius)
% Function for Anderson sequential images, using the dapi image from 
% each set of folders, finding the transformation, transformaing the image, then saving
% the dapi images
% each hybcycle and position will have 3 channels with a dapi
% need to align by dapi and save all the images by position; in the
% agression experiment there are 22 HybCycles and 8 positions, so in one
% position 22 dapi images should be registered

% Need to add a script to get make allHybRegistration to check if the first image in each 
% barcode round is aligned....All images are already aligned to the first
% hyb

% make sure bfmatlab and Fiji.app is in the path
% add path of CompareDotsError, and AlignImages
%addpath('C:\github\streamline-seqFISH\src\AlignImages', '-end');
%addpath(fullfile('C:', 'github', 'streamline-seqFISH', 'src', ...
%    'CompareDotsError'), '-end');

% initialze hybnum
numZSlice = [];
tform = [];
hybIms = [];
    dapiref = [];

    dapiIms = cell(1, length(folderRange));


    for folder = 1:length(folderRange)
        fprintf('Transformation for hybcycle %.0f...\n', folder-1);



        % get the images for each channel
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(basePath, ['HybCycle_' num2str(folder-1)], imageName);
        [allIms, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
        numDapiCh = sizeC;
        numHybCh = sizeC - 1;
        numZSlice = sizeZ;
        tform(folder).folder = [];
        hybChannelArray = 1:numHybCh;
        %hybIms(folder).channel = cell(1, 3);
        hybIms(folder).channel = allIms(1:numHybCh);
        dapiIms{folder} = allIms{numDapiCh};
        
        
        %% Get the tform and apply the transformation
        % Use the dapi transformations and the chromatic aberrations (from barcoded
        % experiments)
        fprintf('Apply Dapi tform for position %.0f and folder %.0f...\n', position, folder);
        if folder ~= 1% if not the first folder

                       % find the dapi transformation
            initialRadius = 0.0625; %0.0625 for 3d is default
            numIterations = 100; % 100 is default
            tformTemp = grabtform(dapiIms{folder}, dapiref, initialRadius, numIterations);


            % get median of the tform
            if numZSlice < 16
                tform(folder).folder = tformTemp;
                for f = 1:numHybCh
                    hybIms(folder).channel{f} = imwarp(hybIms(folder).channel{f}, tform(folder).folder, 'OutputView', imref2d(size(hybIms(folder).channel{f})));
                end
                dapiIms{folder} = imwarp(dapiIms{folder}, tform(folder).folder, 'OutputView', imref2d(size(dapiIms{folder})));
            else
                % remove 3rd dim
                tformTemp.T(4,3) = 0;
                tform(folder).folder = tformTemp;
                for f = 1:numHybCh
                    hybIms(folder).channel{f} = imwarp(hybIms(folder).channel{f}, tform(folder).folder, 'OutputView', imref3d(size(hybIms(folder).channel{f})));
                end
                dapiIms{folder} = imwarp(dapiIms{folder}, tform(folder).folder, 'OutputView', imref3d(size(dapiIms{folder})));
            end
        else
            if numZSlice < 16
                tform(folder).folder = affine2d(eye(3));
            else
                tform(folder).folder = affine3d(eye(4));
            end
            dapiref = dapiIms{1};%imagejbackgroundsubtraction(dapiIms{1}, uniqueString, pixelRadius);
        end      

    end
    
    %% Need dapi images for all hyb registration
    Miji(false);

    for f = 1:length(folderRange)
        namesh{f} = ['C' num2str(f) '-'  num2str(1) '.tif'];
        MIJ.createImage(namesh{f}, dapiIms{f}, true);
    end

    str = [];
    for f = 1:length(folderRange)
        str = [str ' image' num2str(f) '=C' num2str(f) '-' num2str(1) '.tif'];
    end

    try
        MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
        MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(length(folderRange)) ' slices=' num2str(numZSlice) ' frames=1 display=Grayscale']);
        savePath = fullfile(saveDir, ['AllHybRegistrationCheck-Pos' num2str(position)]);
        MIJ.run('Save', ['save=[' savePath '.tif' ']']);
        MIJ.run('Close All')
        MIJ.exit;
    catch
        MIJ.exit;
        error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
    end
    
    %% Save the Hyb images as tif images
    %splitFactorHybImages 4; % split image into 4 pieces for aggression, mating and control exp, split into2 pieces for new experiments
    endingString = ''; % no ending string
    savefolchimage(length(folderRange), numHybCh, numZSlice, position, hybIms, saveDir, splitFactorHybImages, endingString)

    %% save .mat files
    fprintf('Saving images for position %.0f\n', position);
    saveDirImages = fullfile(saveDir, ['imagesHybDapi-pos' num2str(position)]);
    save(saveDirImages, 'hybIms', 'dapiIms', '-v7.3');
end