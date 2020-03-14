function [I, points, intensity, shadingcorr, adjustedThreshold] = outputprocessimagesimmunoUseDapiTform2(experimentDir, ...
    position, folderArray, numCh, typedots, superres, firstThreshold, numRefPoints, ...
        experimentName, experimentLabel, roimask, filtersigma, refSaveName, ...
        githubDir, longestEdges, minEdgeMatch, minDotMatch, chaTform, physicaltform, ...
        threshold, dnafishPath, imageJBackSubtract)
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
    sizeC = numCh + 1;
    sizeCRef = [];
    I = cell(length(folderArray), numCh);
    dapiI = cell(length(folderArray), 1);
    allIms = cell(length(folderArray),1);
    folderArrayIdx = folderArray + 1;
    tic
    
    imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
    imagePath = fullfile(experimentDir, ['HybCycle_' num2str(0)], imageName);
    [~, sizeCRef, sizeZ, ~, ~] = grabimseries(imagePath, position, sizeC);
    
    
    parfor fidx = folderArrayIdx
        fprintf('Retrieving Position %.0f Folder %.0f images\n', position, fidx-1);
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(experimentDir, ['HybCycle_' num2str(fidx-1)], imageName);
        [allIms{fidx}, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position, sizeCRef);
    end
    dapiRef = allIms{1}{sizeCRef};


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
            
            
            if imageJBackSubtract
                % ImageJ Rolling Ball Back Subtract to remove noise using rad 3
                % replace with deconvolution - need to test first
                uniqueString = 'imageTempProcess-3535fsg';
                I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString,...
                    experimentDir);
            end
        end
    end
    


    % get the tform for channel 1, 2 and 3
    if length(chaTform) < numCh
        numCh = length(chaTform);
    end

    newOffsetsPath = getfile(fullfile(experimentDir, 'processedimages',['pos' num2str(position)]), 'dapi-tform-', 'match');
    offsets = readtable(newOffsetsPath);
    offsets.Properties.VariableNames([1 2 3]) = {'hyb' 'col' 'row'};
    fiduciaryOffsets = cell(length(folderArray),numCh);
    
    % apply the final tform to the images
    idx = 1;
    for f = folderArray
        for ch = 1:numCh
            
            fiduciaryOffsets{idx,ch} = maketform2(offsets.col(idx), offsets.row(idx), offsets.z(idx));
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
    for ch = 1:numCh
        savefolchimage(position, Icut(:,ch), saveImageDir, 'alignmentcheck-', ['ch' num2str(ch)])
    end
end