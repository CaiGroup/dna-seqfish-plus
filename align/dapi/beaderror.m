function rawThreshold = beaderror(saveFilePath, saveDir, hybIndex, rawHybIms, ...
    beadAlignedImage, beadPoints, fixedPoints, movingMatchPoints, tformBeads, ...
    tformDapi, position, typedots, backradius, sliding, rawThreshold, targetNumPoints)
% gets the error between the dapi aligned image and the beads aligned image
% and outputs the final tform for dapi and beads.
%
% inputs: refPoints, rawImage, dapiAlignedImage
%
% Assumptions:
% 1. fixed match points and moving match points are match already
% 2. 
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 8/12/2019

    
    % open csv file to be rewritable and add to the csv file
    fileID = fopen(saveFilePath, 'a+'); % open and append
    fprintf(fileID, '\n\n--------------------\n');
    fprintf(fileID, 'HybCycle %.0f:\n', hybIndex);
    
    %% Get Points from the Dapi Images
    processHybImsDapi = true;
    %threshold = []; % use the autothreshold
    %numRefPoints = size(fixedPoints, 1);
    %targetNumPoints = numRefPoints * 5;
    matchDapiRadius = 3;
    [hybDapiPoints, rawThreshold] = matchpointsims(rawHybIms, rawThreshold, ...
        typedots, processHybImsDapi, backradius, sliding, targetNumPoints, matchDapiRadius);
    
    
    %% Before Dapi Alignment Error 
    matchRadius = 5;
    [~, ~, rawPoints, ~] = matchpoints(fixedPoints, hybDapiPoints(1).channels, matchRadius);
    rawError1 = mean(rawPoints.distance);
    rawErrorX1 = mean(abs(rawPoints.ref(:,1) - rawPoints.points(:,1)));
    rawErrorY1 = mean(abs(rawPoints.ref(:,2) - rawPoints.points(:,2)));
    rawErrorZ1 = mean(abs(rawPoints.ref(:,3) - rawPoints.points(:,3)));
    fprintf(fileID, 'Before Dapi Alignment Error:\n');
    fprintf(fileID, 'Beads Colocalized:,%.0f\n', size(rawPoints.ref,1));
    fprintf(fileID, 'Distance:,%.3f \n', rawError1);
    fprintf(fileID, 'X:,%.3f \n', rawErrorX1);
    fprintf(fileID, 'Y:,%.3f \n', rawErrorY1);
    fprintf(fileID, 'Z:,%.3f \n', rawErrorZ1);
    
    
    
    %% After Dapi Alignment Error 
    rawError2 = mean(movingMatchPoints.distance);
    rawErrorX2 = mean(abs(movingMatchPoints.ref(:,1) - movingMatchPoints.points(:,1)));
    rawErrorY2 = mean(abs(movingMatchPoints.ref(:,2) - movingMatchPoints.points(:,2)));
    rawErrorZ2 = mean(abs(movingMatchPoints.ref(:,3) - movingMatchPoints.points(:,3)));
    fprintf(fileID, 'After Dapi Alignment Error:\n');
    fprintf(fileID, 'Dapi Tform is:\n');
    fprintf(fileID, '%.2f,%.2f,%.2f,%.2f\n', tformDapi.T'); % invert because fprintf prints columnwise
    fprintf(fileID, 'Beads Colocalized:,%.0f\n', size(movingMatchPoints.ref,1));
    fprintf(fileID, 'Distance:,%.3f\n', rawError2);
    fprintf(fileID, 'X:,%.3f \n', rawErrorX2);
    fprintf(fileID, 'Y:,%.3f \n', rawErrorY2);
    fprintf(fileID, 'Z:,%.3f \n', rawErrorZ2);
    
    
    
    %% After Dapi and Bead Alignment Error 
    rawError3 = mean(sqrt(sum((beadPoints - movingMatchPoints.points) .^ 2, 2))); % movingMatchPoints.points is the fixed points
    rawErrorX3 = mean(abs(beadPoints(:,1) - movingMatchPoints.points(:,1)));
    rawErrorY3 = mean(abs(beadPoints(:,2) - movingMatchPoints.points(:,2)));
    rawErrorZ3 = mean(abs(beadPoints(:,3) - movingMatchPoints.points(:,3)));
    fprintf(fileID, 'After Dapi and Bead Alignment Error:\n');
    fprintf(fileID, 'Bead Tform is:\n');
    fprintf(fileID, '%.2f,%.2f,%.2f,%.2f\n', tformBeads.T');
    fprintf(fileID, 'Beads Colocalized:,%.0f\n', size(beadPoints,1));
    fprintf(fileID, 'Distance:,%.3f\n', rawError3);
    fprintf(fileID, 'X:,%.3f \n', rawErrorX3);
    fprintf(fileID, 'Y:,%.3f \n', rawErrorY3);
    fprintf(fileID, 'Z:,%.3f \n', rawErrorZ3);
    % close file
    fclose(fileID);
    
    
    
    %% Figure for all three points
    close all
    H1 = figure;
    imshow(max(beadAlignedImage{1}, [], 3), 'DisplayRange', [0 2000], 'InitialMagnification', 'fit');
    hold on
    scatter(hybDapiPoints(1).channels(:,1), hybDapiPoints(1).channels(:,2), 150, 'gs', 'LineWidth', 1);
    scatter(fixedPoints(:,1), fixedPoints(:,2), 125, 'cx', 'LineWidth', 1);
    scatter(beadPoints(:,1), beadPoints(:,2), 150, 'm.', 'LineWidth', 1);
    legend('none','dapi','bead','Location','southeast');
    hold off
    save3Name = fullfile(saveDir, ['Hyb-' num2str(hybIndex) '-Pos' num2str(position) '-raw-dapi-bead-error.fig']);
    savefig(H1, save3Name);
    close all;

end