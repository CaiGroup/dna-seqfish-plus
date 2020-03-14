function [I, hybIms, dapiTform] = preprocessimagesRot_v2(experimentName, experimentDir, ...
    position, folderArray, varargin)
% preprocessimages takes raw images from an experiment and preprocesses the
% images: align by dapi, subtract the background, and background subtract
% using the imageJ rolling ball algorithm.
%
% Inputs: experiment name (used for naming files), experiment directory
% (main directory with HybCycle_[i] images, etc.), position or fov, folder array
% (array for number of folders ex. 0:4 for [0,1,2,3,4] folders.
%
% Output: processed images; saves the images in the project folder as a
% tiff file and a mat data file
%
% Requirements: images need at least 4 Z-slices for the alignment
%
% Options: decide to use background images or not - usually used to remove
% tissue background or autofluorescence.
%
% addpath('C:\github\streamline-seqFISH\src\preprocessing\bfmatlab', '-end');
%
% Dependencies: 
% 1. Miji.jar is in the Fiji.app\scripts directory. If not, download and
% install at: ......
% 2. Functions imagejbackgroundsubtraction, grabimseries, shadingcorrection,
% imdivide, grabtform, getdirectory
%
% Time: 3 field of views (without saving hyb processed images) took 5 hours
% - increase speed using parfor loops or increase read/write speed for
% saving.
%
% To Do:
% Add an option to save images; default is on
%
% Options: 
% 1. useBackgroundImages: boolean to use background images for subtraction
% 2. backgroundFolderName: change default folder name in experimentDir
% from 'initial_background' to new folder
% 3. dapiRefPath: path to images that will be used as the dapi reference to
% align images.
% 4. imageJBackSubtract: boolean to use imageJ rolling ball background
% subtraction.
% 5. subtractBackground: boolean to subtract background - useful for high
% bakcground images or autofluorescence.
% 6. saveProcessedHybIms: boolean to save processed Hyb images
% default Values for optional arguments:
% 7. divideIms: divides dapi images into 4 to use gradient descent
% for 3d alignment - good for stacks with less than 16 zslices
% [true, 'initial_background', [], false, false, false]; 
% 
% dapiRefPath should be HybCycle_0 for RNA, not DNA.
% dapiRefPath2 should be HybCycle_0 for DNA.
%
% correct rotation after usual tform alignment to HybCycle_0 RNA 
% and then apply global tfrom between HybCycle_0 RNA and DNA and
% rotation tfrom, calculated from aligned HybCycle_0 RNA and DNA.
%
% Date: 01/07/2020
% Author: Yodai Takei
% ytakei@caltech.edu

    tic
    %% Check if Fiji and bfmatlab is in the path
    fijiDirectory = checkfijipath();
    checkbfmatlabpath();
    
    %% Set up optional Parameters
    numvarargs = length(varargin);
    argsLimit = 9;
    if numvarargs > argsLimit
        error('myfuns:preprocessimages:TooManyInputs', ...
            'requires at most 7 optional inputs');
    end
    % set defaults for optional inputs
    optargs = {true, 'initial_background', '', true, false, false, false, '3d', ''}; 
    optargs(1:numvarargs) = varargin;
    % Place optional args in memorable variable names
    [useBackgroundImages, backgroundFolderName, dapiRefPath, imageJBackSubtract, subtractBackground, ...
        saveProcessedHybIms, divideIms, dim, dapiRefPathDNA] = optargs{:};

    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = datestr(dateStart, formatDate);
    
    
    
    %% Initialize Variables
    if useBackgroundImages
        if ~exist(fullfile(experimentDir, backgroundFolderName), 'dir')
            error('Background Image Folder not found...\n1. Change optional variables to I = preprocessimages(experimentName, experimentDir, fovArray, folderArray, false)\n2.Rename folder I = preprocessimages(experimentName, experimentDir, fovArray, folderArray, true, folderName)');
        end
    end
    folderSaveName = 'processedimages';


    % Make save directory: default 'organizehybs' folder
    saveDir = fullfile(experimentDir, folderSaveName, ['pos' num2str(position)]);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end


    %% Align the Dapi for the hybs
    fprintf('Aligning Dapi for Hyb Images...\n');
    % this function saves the images for the alignment of the images
    [hybIms, dapiIms, dapiTform, numCh] = aligndapiimages(experimentName, experimentDir, folderArray, position, saveDir, dapiRefPath, true, false, divideIms, dim);
    % save dapiTform as a csv
    saveDapiTformPath = fullfile(saveDir, ['dapi-tform-' experimentName '-' endingDateString '.csv']);
    printtform(dapiTform, saveDapiTformPath);
    
    %% Apply rotation correction using hyb1 DAPI vs ref DAPI tform, calculated from already aligned DAPI images above.
    [hybIms, dapiIms, tformDapiRNADNA, tformDapiRot] = aligndapiimagesRotglobal_v2(experimentName, ...
    folderArray, position, saveDir, hybIms, dapiIms, dapiRefPathDNA);
    saveDapiTformPath2 = fullfile(saveDir, ['dapi-rotation-tform-' experimentName '-' endingDateString '.mat']);
    save(saveDapiTformPath2, 'tformDapiRNADNA', 'tformDapiRot');

    %% Use Background Images or just Background Subtract
    I = cell(length(folderArray), numCh);
    if useBackgroundImages % NEEDS DEBUGGIN - ERROR HERE at line ,159
        %% Align Dapi for Background Images for All Positions
        fprintf('Aligning Dapi for Background Images...\n');
        backImBasePath = fullfile(experimentDir, backgroundFolderName);
        listing = dir([experimentDir '\' backgroundFolderName '\*MMStack_Pos' num2str(position) '.ome.tif']);
        backImPath = fullfile(backImBasePath, listing(1).name);
        %backImPath = fullfile(backImBasePath, ['MMStack_Pos' num2str(position) '.ome.tif']);
        if exist(backImPath, 'file') == 2
            [backIms, numDapi, numZSlice, ~, ~] = grabimseries(backImPath, position);
            % dapi is the 4 in the cell
            numCh = numDapi - 1;

            if numZSlice < 16
                % divide image into 4 pieces from 1 image if zslices < 16
                [backImDivide, ~] = imdivideby4(backIms{numDapi});
                [backImRefDivide, numZSliceDivide] = imdivideby4(dapiIms{1});
            else
                backImDivide = backIms{numDapi};
                backImRefDivide = dapiIms{1};
                numZSliceDivide = 15; % arbitrary
            end

            initialRadius = 0.0625; %0.0625 for 3d is default
            numIterations = 100; % 100 is default
            backTform = grabtformRot(backImDivide, backImRefDivide, initialRadius, numIterations);
            % remove the z transformation
            if numZSliceDivide >= 16
                backTform.T(4,3) = 0;
            end
            % apply the tform
            for ch = 1:numCh
                if length(backTform.T) == 3
                    outputView = imref2d(size(backIms{ch}));
                elseif length(backTform.T) == 4
                    outputView = imref3d(size(backIms{ch}));
                end
                backIms{ch} = imwarp(backIms{ch}, backTform, 'OutputView', outputView);
            end



            %% Subtract the background and Multiply by Shading Corrections
            % Get Shading Corrections
            clearvars dapiIms
            shadingcorr = shadingcorrection(backIms(1:numCh));
            for f = 1:length(folderArray)
                for ch = 1:numCh
                    if subtractBackground % option to subract background
                        imageTemp = backsubtract(hybIms{f,ch}, backIms{ch});
                    else
                        imageTemp = hybIms{f,ch};
                    end
                    prct5 = prctile(backIms{ch}(:),5);
                    % Remove Inf double values and set to percentile 5
                    infInd = find(imageTemp == Inf);
                    zeroInd = find(imageTemp == 0);
                    if ~isempty(infInd)
                        imageTemp(ind2sub(size(imageTemp),infInd)) = prct5;
                    end
                    if ~isempty(zeroInd)
                        imageTemp(ind2sub(size(imageTemp),zeroInd)) = prct5;
                    end

                    % Apply the shading correctionsmean
                    imageTemp = double(imageTemp) ./ double(shadingcorr{ch});
                    I{f,ch} = uint16(imageTemp);

                    if imageJBackSubtract
                        % ImageJ Rolling Ball Back Subtract to remove noise using rad 3
                        % replace with deconvolution - need to test first
                        uniqueString = 'imageTempProcess-90jf03j';
                        I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString,...
                            experimentDir);
                    end

                end
            end

            % Save Mat Files
            savePath = fullfile(saveDir, ['preProcessedData-pos' num2str(position) '-' experimentName '-' endingDateString '.mat']);
            save(savePath, 'I', 'backIms', 'backTform', 'dapiTform', 'shadingcorr', 'position','-v7.3');
        else
            error('background image: %s\nNot Found', backImPath);
        end
    else
        %% Don't use Background Images; only background subtract
        for f = 1:length(folderArray)
            for ch = 1:numCh
                I{f,ch} = hybIms{f, ch};

                if imageJBackSubtract
                    % ImageJ Rolling Ball Back Subtract to remove noise using rad 3
                    % replace with deconvolution - need to test first
                    % To do: update imagejbackgroundsubtraction.m
                    uniqueString = 'imageTempProcess-90jf03j';
                    I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString, experimentDir);
                end

            end
        end

        if imageJBackSubtract % no need to save if no background subtraction is used; only need hybIms from aligndapiimages
            % Save Mat Files
            savePath = fullfile(saveDir, ['preProcessedData-pos' num2str(position) '-' experimentName '-' endingDateString '.mat']);
            save(savePath, 'I', 'dapiTform', 'position', '-v7.3');
        end
    end


    % Save Hyb Images as Tiff files
    if saveProcessedHybIms
        startString = ['preprocessedImages-' experimentName '-'];
        savefolchimage(position, I, saveDir, startString, endingDateString)
    end


    toc
end