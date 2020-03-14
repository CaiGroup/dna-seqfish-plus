function [] = findpointsfromthreshold(experimentDir, experimentName, position, folderArray, threshold)
% function to get points from the images after being preprocessed

% take out multiplier for other experimetns; only using multiplier because
% the spot detection algorithms were different.

addpath('C:\github\streamline-seqFISH\src\preprocessing', '-end');
addpath('C:\github\streamline-seqFISH\src\preprocessing\bfmatlab', '-end');
saveDir = fullfile(experimentDir, 'points');
if exist(saveDir, 'dir') ~= 7
    mkdir(saveDir);
end
useBackgroundImages = true;
imageJBackSubtract = true;
numCh = 3;
backgroundFolderName = 'initial_background';
subtractBackground = false;
I = cell(length(folderArray), numCh);
dapiI = cell(length(folderArray),1);
% get the images for each channel
 for folder = 1:length(folderArray)

     %% Use Background Images or just Background Subtract
        
        fprintf('Retrieving Position %.0f Folder %.0f images\n', position, folderArray(folder));
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(experimentDir, ['HybCycle_' num2str(folderArray(folder))], imageName);
        [allIms, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
        I(folder,:) = allIms(1:numCh);
        dapiI{folder} = allIms(sizeC);
        if folder == 1
            dapiIms = allIms{sizeC};
        end

 end

%% Align Dapi for Background Images for All Positions
fprintf('Aligning Dapi for Background Images...\n');
backImBasePath = fullfile(experimentDir, backgroundFolderName);
backImPath = fullfile(backImBasePath, ['MMStack_Pos' num2str(position) '.ome.tif']);
if exist(backImPath, 'file') == 2
    [backIms, numDapi, numZSlice, ~, ~] = grabimseries(backImPath, position);
    % dapi is the 4 in the cell
    numCh = numDapi - 1;

    if numZSlice < 16
        % divide image into 4 pieces from 1 image if zslices < 16
        [backImDivide, ~] = imdivideby4(backIms{numDapi});
        [backImRefDivide, numZSliceDivide] = imdivideby4(dapiIms);
    else
        backImDivide = backIms{numDapi};
        backImRefDivide = dapiIms;
        numZSliceDivide = 15; % arbitrary
    end

    initialRadius = 0.0625; %0.0625 for 3d is default
    numIterations = 100; % 100 is default
    backTform = grabtform(backImDivide, backImRefDivide, initialRadius, numIterations);
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
    shadingcorr = shadingcorrection(backIms(1:numCh));
    for f = 1:length(folderArray)
        for ch = 1:numCh
            if subtractBackground % option to subract background
                imageTemp = backsubtract(I{f,ch}, backIms{ch});
            else
                imageTemp = I{f,ch};
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
else
    error 'background directory or images were not found';
end


%% get the points
points = cell(length(folderArray),numCh);
intensity = cell(length(folderArray), numCh);
sqrtintensity = cell(length(folderArray), numCh);
xyPixSize = 1;
zPixSize = 1;
maxXY = 2048;
maxZ = size(I{f,ch}, 3);
min = 1;
for f = 1:length(folderArray)
    for ch = 1:numCh
        [pointsTemp, ~, ~, ~] = detectdotsv2(I{f, ch}, threshold(f,ch), 'log', false, '', 0.72);
        [points{f,ch}, sqrtintensity{f,ch}] = SuperResPoints(pointsTemp,I{f,ch},xyPixSize,zPixSize);
        %points{f,ch} = getgaussian(pointsTemp, I{f,ch});
        %intensity = zeros(length(points{f,ch}), 1);
        %removeDots = [];
        
        % round points and get intensity
        tempPoints = points{f,ch};
        for i = 1:length(points{f,ch})
            y = round(tempPoints(i,2));
            x = round(tempPoints(i,1));
            z = round(tempPoints(i,3));
            if z > maxZ
                z = maxZ;
            elseif z < min
                z = min;
            end
            if x > maxXY
                x = maxXY;
            elseif x < min
                x = min;
            end
            if y > maxXY
                y = maxXY;
            elseif y < min
                y = min;
            end
            intensity{f,ch}(i) = I{f,ch}(y, x, z);
        end
    end
end
        
savePath = fullfile(saveDir, ['hyb-points-0_80-ForBeadAlignment-pos' num2str(position) '-radialcenter3d-' experimentName '.mat']);
save(savePath, 'points', 'intensity', 'backTform', 'shadingcorr');

listSavePath = fullfile(saveDir, ['hyb-points-0_80-ForBeadAlignment-pos' num2str(position) '-radialcenter3d-' experimentName '.csv']);
fileID = fopen(listSavePath,'w');
fprintf(fileID, '%s,%s,%s,%s,%s,%s\n', 'ch', 'hyb', 'x', 'y', 'z', 'int');

for hyb = 1:length(folderArray)
    for ch = 1:numCh
        pointsSize = length(points{hyb,ch});
        for i = 1:pointsSize
            x = points{hyb,ch}(i,1);
            y = points{hyb,ch}(i,2);
            z = points{hyb,ch}(i,3);
            int = intensity{hyb,ch}(i);
            fprintf(fileID, '%.0f,%.0f,%.3f,%.3f,%.3f,%.0f\n', ch, hyb, x, y, z, int);
        end
    end
end

%{
% save the dapi images
saveDapiPath = fullfile(saveDir, ['dapiIms-pos' num2str(position) '-' experimentName '.mat']);
save(saveDapiPath, 'dapiI');
% save the images for each channel
for ch = 1:numCh
    saveImPath = fullfile(saveDir, ['hybIms-pos' num2str(position) '-ch' num2str(ch) '-' experimentName '.mat']);
    chI = I(:,ch);
    save(saveImPath, 'chI');
end
%}


end
