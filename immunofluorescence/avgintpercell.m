function [] = avgintpercell(experimentDir, experimentName, experimentLabel, ...
    I, L, position, folderArray, chArray, segoption, binSize)
% function to bin pixels based on the segmentation
%
% try getting the images of 42-43 - NEED to add the illumination correction
% to make sure the intensities of all the images are correct and non-biased
% 1. Align by dapi
% 2. get the cell boundaries in 3d
% 3. Get the intensities
% 4. Put into the Matrix
%
% Dependencies:
% 1. Processed Images
% 2. Segmentation in 2d or 3d
%
% Path Dependencies:
% addpath('C:\github\streamline-seqFISH\src\preprocessing', '-end');
% addpath('C:\github\streamline-seqFISH\src\preprocessing\bfmatlab\', '-end');
% addpath('C:\github\streamline-seqFISH\src\AlignImages\testscripts', '-end');
%
% Bug Report:
% 8/26/2019 - Valid File ID not found when using fopen....only decoded up
% to hyb19 from a folder array of 0:19, which should go until hyb20
% Error Message -----------------------------------------------------------
% Error using fprintf
% Invalid file identifier. Use fopen to generate a valid file
% identifier.
% Error in binpixelspercell (line 18)
%     fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s\n', ...
% Error in avgintpercell (line 81)
%             binpixelspercell(binSize, statsFilt, shape2d, image,
%             position, folder, channel, savePath, segoption,
%             endingString);
%-------------------------------------------------------------------------- 
%
% Example:
% ran preprocessimages.m to get I
%
% experimentDir = 'I:\2019-07-25-E14-DNA-seqFISH+rep2-2-DNAFISH-plate2';
% experimentName = 'ImmunoFluorescenc-2019-07-25-E14-DNA-seqFISH+rep2-2-DNAFISH-plate2';
% position = 0
% chArray = 1:2;
% folderArray = 4:19; % choose hybs to look at; ex: 4:19
% segoption = '2d';
%
% avgintpercell(experimentDir, experimentName, I, position, folderArray, chArray, segoption);
%
% Date: 8/26/2019
% Author: Nico Pierson
% Email: nicogpt@caltec.edu
    

    %% need to get voxel size from the image
    voxelSize = 0.1108; % in um


    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    dateString = [datestr(dateStart, formatDate)];


    
    %% Use the Segmentation to get cells
    saveFolder = 'analysis';
    savePath = fullfile(experimentDir, saveFolder, experimentLabel);
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end

    segPath = fullfile(experimentDir, 'segmentation');
    switch segoption
        case '3d'
            %multiName = 'testHyb38dapi_Multicut Segmentation.h5'; % need to change this %%%%%
            %simpleName = 'testHyb38dapi_Simple Segmentation.h5';
            %stainName = 'dapi';
            %[statsFilt, boundaries, shape2d] = getsegmentation(segPath, multiName, simpleName, position, savePath, experimentName, stainName);
            [statsFilt, boundaries, shape2d] = getsegmentationavgintdapi(L, I{1,4}, position, savePath, experimentName);
        case '2d'
            roiSegPath = fullfile(segPath, ['Pos' num2str(position)], 'RoiSet.zip');
            shape2d = selfsegzip(roiSegPath);
            statsFilt = [];
            boundaries = []; % later add boundaries using polygon
            %% save data per position
            numCells = numel(fieldnames(shape2d));
            saveFileName = ['boundaries2d-Data-pos' num2str(position)];
            saveDataPath = fullfile(savePath, saveFileName);
            save(saveDataPath, 'numCells', 'statsFilt', 'boundaries', 'shape2d');
    end
    
    
    
    % set up csv file
    endingString = ['-' dateString '.csv'];%['-preData-' dateString '.csv'];
    sumcellxhybs = cell(length(chArray), 1);
    chIter = 1;
    for channel = chArray % want vector for each channel for each cell
        
        % Get specific image based on folder and channel
        image = I(:,channel);


        % calculate the average pixels
        %binSize = 2; % bin by 2
    
        % save for each cell
        % add 1 to folder because it starts at 0
        %sumcellxhybs{chIter} = binpixelspercell(binSize, statsFilt, shape2d, image, savePath, endingString, position, folderArray + 1, channel, segoption);
        sumcellxhybs{chIter} = binpixelspercellparallel(binSize, statsFilt, shape2d, image, savePath, endingString, position, folderArray + 1, channel, segoption);
        %endingString = ['-preliminaryData-' dateString];
        %binpixelspercell(binSize, statsFilt, shape2d, image, position, folder, channel, savePath, segoption, endingString);
        chIter = chIter + 1;
    end
    
    %% print the sum of the pixels using the sumcellxhybs structure
    % set up csv file
    listSavePath = fullfile(savePath, ['sumData-ImmunoFluorescence-Pos' num2str(position) endingString]);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID,'%s,%s,%s,%s,%s\n', 'chID', 'cellID', 'hybID', 'sumint', 'voxelsize');
    
    
    chIter = 1;
    for ch = chArray
        numCells = size(sumcellxhybs{chIter}, 1);
        for c = 1:numCells
            hybIter = 1;
            for h = folderArray + 1 % add 1 because index starts at 0 for input
                % for non parallel function bin pixels per cell
                %fprintf(fileID,'%d,%d,%d,%.2f,%.4f\n', ch, c, h, sumcellxhybs{chIter}{c,hybIter}, voxelSize);
                fprintf(fileID,'%d,%d,%d,%.2f,%.4f\n', ch, c, h, sumcellxhybs{chIter}{c}{hybIter}, voxelSize);
                hybIter = hybIter + 1;
            end
        end
        chIter = chIter + 1;
    end
                
    fclose(fileID);
    
    
end
