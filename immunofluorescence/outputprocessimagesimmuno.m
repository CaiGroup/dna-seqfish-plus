function [I, points, intensity, shadingcorr, adjustedThreshold] = outputprocessimagesimmuno(experimentDir, ...
    position, folderArray, numCh, typedots, superres, firstThreshold, numRefPoints, ...
        experimentName, experimentLabel, roimask, filtersigma, refSaveName, ...
        githubDir, longestEdges, minEdgeMatch, minDotMatch, chaTform, physicaltform, ...
        threshold, dnafishPath)
% output the processeed images depending on the position to output the ch3
% reference points using the number of beads from the beads

    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end

    saveFigDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', ['points_debug_pos' num2str(position)]);
    if exist(saveFigDir, 'dir') ~= 7
        mkdir(saveFigDir);
    end
    
    saveImageDir = fullfile(experimentDir, 'analysis', experimentLabel);
    
    
    %% Declare Variables
    selectedCh = 3;
    points = cell(length(folderArray),1);
    intensity = cell(length(folderArray), 1);
    intcheck = cell(length(folderArray), 1);
    sigma = cell(length(folderArray), 1);
    xyPixSize = 1;
    zPixSize = 1;
    maxXY = 2048;
    maxZ = [];
    min = 1;
    adjustedPoints = cell(length(folderArray),1);
    adjustedIntensity = adjustedPoints;
    medianError = adjustedPoints;
    adjustedThreshold = ones(length(folderArray),1) * 999999;
    numPointError = ones(length(folderArray),1) * 999999;
            


    %% load pos1-4 images, grab points using same treshold vales
    % get images and grab points
    sizeC = [];
    I = cell(length(folderArray), numCh);
    dapiI = cell(length(folderArray), 1);
    allIms = cell(length(folderArray),1);
    folderArrayIdx = folderArray + 1;
    tic
    parfor fidx = folderArrayIdx
        fprintf('Retrieving Position %.0f Folder %.0f images\n', position, fidx-1);
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(experimentDir, ['HybCycle_' num2str(fidx-1)], imageName);
        [allIms{fidx}, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
    end
    imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
    imagePath = fullfile(experimentDir, ['HybCycle_' num2str(0)], imageName);
    [~, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
    dapiRef = allIms{1}{sizeC};


    %% Align Dapi for Background Images for All Positions
    fprintf('Aligning Dapi for Background Images...\n');
    backgroundFolderName = 'initial_background';
    backImBasePath = fullfile(experimentDir, backgroundFolderName);
    backImPath = fullfile(backImBasePath, ['MMStack_Pos' num2str(position) '.ome.tif']);
    if exist(backImPath, 'file') == 2
        [backIms, numDapi, maxZ, backTform] = alignbackimages(dapiRef, backImPath, position);


        %% Subtract the background and Multiply by Shading Corrections
        % Get Shading Corrections
        shadingcorr = shadingcorrection(backIms(1:numDapi-1));
    else
        backTform = [];
        shadingcorr = [];
    end
    for f = 1:length(folderArray)
        for ch = 1:numCh
            fprintf('folder: %.0f; channel %.0f\n', f, ch);
            % Apply the shading correctionsmean
            if ~isempty(shadingcorr)
                I{f,ch} = uint16(double(allIms{f}{ch}) ./ double(shadingcorr{ch}));
            else
                I{f,ch} = uint16(allIms{f}{ch});
            end
            if ch == selectedCh
                % if there is a threshold
                if isempty(threshold)
                    %% Get the points using the number of reference points
                    [~, ~, adjustedThreshold(f), numPointError(f)] = ...
                    autothresholdbynumberpoints(I{f,ch}, numRefPoints, firstThreshold, ...
                    typedots, roimask);
                else
                    adjustedThreshold = threshold;
                end
                
                
                %% grab the points
                saveFigPath = fullfile(saveFigDir, ['figure-Pos' num2str(position) '-Folder' num2str(f) '-Ch' num2str(ch) '.fig']);
                [pointsTemp2, intcheck{f}, ~, ~] = detectdotsv2(I{f,ch}, adjustedThreshold(f), typedots, true, saveFigPath, 1);
                switch superres
                    case 'radial3d'
                    [points{f}, intensity{f}, sigma{f}] = SuperResPoints(pointsTemp2,I{f,ch},xyPixSize,zPixSize);
                    case 'gaussian'
                    [points{f}, intensity{f}] = getgaussian(pointsTemp2, I{f,ch});
                    [~, ~, sigma{f}] = SuperResPoints(pointsTemp2,I{f,ch},xyPixSize,zPixSize, false);
                end
                % filter using sigma
                if filtersigma
                    removeidx = [];
                    for i = 1:length(points{f})
                        sigcheck = sigma{f}(i);
                        int = intensity{f}(i);
                        if sigcheck < 0.6 || int < 500
                            removeidx = cat(1, removeidx, i);
                        end
                    end
                    points{f}(removeidx, :) = [];
                    intensity{f}(removeidx, :) = [];
                    sigma{f}(removeidx, :) = [];
                end
            end
        end
    end
    

    numAlignCh = 1;
    % save the points in a csv file
    [rawpoints, intensity, hybSaveName] = outputrawpointsch3(points, intensity, adjustedThreshold, ...
        folderArray, chaTform, experimentDir, experimentName, position, ...
        numAlignCh, shadingcorr, experimentLabel);
    
    
    % get the offsets from the points - debug this
    saveName = 'point-ch3-offsets_pos%d_';
    saveNameOut = sprintf(saveName, position);
    dnafishSaveName = sprintf('point-offsets_pos%d_', position);
    fiduciaryalign(experimentDir, experimentLabel, refSaveName, hybSaveName, saveNameOut, ...
        githubDir, position, longestEdges, minEdgeMatch, minDotMatch, numAlignCh, length(folderArray), dnafishPath);
    
    % get the offsets 
    offsetsPath = getfile(fullfile(saveDir, 'positions'), saveNameOut, 'match');
    offsetsT = readtable(offsetsPath);
    dnafishOffsetsPath = getfile(fullfile(dnafishPath, 'points', 'positions'), dnafishSaveName, 'match');
    dnafishOffsetsT = readtable(dnafishOffsetsPath);
    
    % add offsets and save the images as processed images
    fiduciaryOffsets = cell(length(folderArray),numCh);
    finalOffset = cell(length(folderArray),numCh);
    offsets = offsetsT(offsetsT.ch == numAlignCh,:);
    
    % reference offsets from hyb 1 channel 1
    dnafishOffsets = dnafishOffsetsT(dnafishOffsetsT.ch == 1,:); 
    offsets.row(:) = offsets.row(:) - dnafishOffsets.row(1);
    offsets.col(:) = offsets.col(:) - dnafishOffsets.col(1);
    offsets.z(:) = offsets.z(:) - dnafishOffsets.z(1);

    
    % get the tform for channel 1, 2 and 3
    if length(chaTform) < numCh
        numCh = length(chaTform);
    end

    
    % apply the final tform to the images
    idx = 1;
    for f = folderArray
        for ch = 1:numCh
            
            fiduciaryOffsets{idx,ch} = maketform2(-offsets.col(idx), -offsets.row(idx), -offsets.z(idx));
            finalOffset{idx,ch} = addtform(chaTform{ch}, fiduciaryOffsets{idx,ch});
            % apply to each image
            I{idx,ch} = imwarp(I{idx,ch}, finalOffset{idx,ch}, 'OutputView', imref3d(size(I{idx,ch})));
            
        end
        idx = idx + 1;
    end
    
    % save the processed and aligned images for ch1 and ch2
    saveImName = ['aligned-processed-ch1-I-pos' num2str(position) '-' experimentLabel '-' experimentName '.mat'];
    saveImPath = fullfile(saveImageDir, saveImName);
    I1 = I(:,1);
    save(saveImPath, 'I1', '-v7.3');
    clearvars I1
    saveImName = ['aligned-processed-ch2-I-pos' num2str(position) '-' experimentLabel '-' experimentName '.mat'];
    saveImPath = fullfile(saveImageDir, saveImName);
    I2 = I(:,2);
    save(saveImPath, 'I2', '-v7.3');
    clearvars I2
    
    % save the variables
    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'points');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end
    savePath = fullfile(saveDirPath, ['points-int-thresh-ch3-pos' num2str(position) '-' experimentLabel '.mat']);
    save(savePath, 'points', 'intensity', 'adjustedThreshold', 'shadingcorr', 'sigma', 'backTform', 'fiduciaryOffsets', ...
        'chaTform', 'finalOffset');
    toc

    % check images only 10-15 zslices
    idx = 1;
    Icut = cell(length(folderArray), numCh);
    for f = folderArray
        for ch = 1:numCh
            
            Icut{idx,ch} = I{idx,ch}(:,:,10:15);
            
        end
        idx = idx + 1;
    end 
    savefolchimage(position, Icut(:,3), saveImageDir, 'alignmentcheck-', 'ch3')
end