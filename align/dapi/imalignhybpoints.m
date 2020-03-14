function [movingFinalImage, fixed, tformPoints] = imalignhybpoints(moving, fixed, movingRawImage, movingPoints, fixedPoints, movingThreshold)
% function imalignhybpoints takes the images and aligned points and outputs
% the error of the current point locations, calculates the error using
% gaussian fitting, then uses closest points to make a transformatin and
% calculates the error
%
% Images are assumed to be in the raw format
% points1 and 2 are expected to have the same number of points

    
    %% Grab points for the raw moving image
    viewPoints = false;
    radius = 3;
    typedots = 'introns';
    [movingRawPoints, ~, ~] = detectdots(movingRawImage, movingThreshold, 'introns', viewPoints);
    
    
    
    %% Gaussian Adjust the dots using the points and the raw images
    movingGauss = getgaussian(movingPoints, moving);
    fixedGauss = getgaussian(fixedPoints, fixed);
    movingRawGauss = getgaussian(movingRawPoints.location, movingRawImage);
    
    
    
    %% Calculate and print the error befor the dapi
    [~, ~, pointsRawMatch, ~] = matchpoints(fixedGauss, movingRawGauss, radius);
    % Print the Error
    rawError = mean(pointsRawMatch.distance);
    rawErrorX = mean(abs(pointsRawMatch.ref(:,1) - pointsRawMatch.points(:,1)));
    rawErrorY = mean(abs(pointsRawMatch.ref(:,2) - pointsRawMatch.points(:,2)));
    rawErrorZ = mean(abs(pointsRawMatch.ref(:,3) - pointsRawMatch.points(:,3)));
    fprintf('\nGaussian Points Before Dapi Alignment:\n');
    fprintf('Gaussian point distance error: %.3f \n', rawError);
    fprintf('X gaussian point error: %.3f \n', rawErrorX);
    fprintf('Y gaussian point error: %.3f \n', rawErrorY);
    fprintf('Z gaussian point error: %.3f \n', rawErrorZ);
    
    
    
    %% Calculate and print the error
    [~, ~, pointsDapiMatch, ~] = matchpoints(fixedGauss, movingGauss, radius);
    gaussError = mean(pointsDapiMatch.distance);
    gaussErrorX = mean(abs(pointsDapiMatch.ref(:,1) - pointsDapiMatch.points(:,1)));
    gaussErrorY = mean(abs(pointsDapiMatch.ref(:,2) - pointsDapiMatch.points(:,2)));
    gaussErrorZ = mean(abs(pointsDapiMatch.ref(:,3) - pointsDapiMatch.points(:,3)));
    fprintf('\nGaussian Points After Dapi Alignment:\n');
    fprintf('Gaussian point distance error: %.3f \n', gaussError);
    fprintf('X gaussian point error: %.3f \n', gaussErrorX);
    fprintf('Y gaussian point error: %.3f \n', gaussErrorY);
    fprintf('Z gaussian point error: %.3f \n', gaussErrorZ);
    
    
    
    %% Transform image2 or the moving image based on points
    tformPoints = getglobaltform(fixedGauss, movingGauss);
    fprintf('Points Tform is: \n');
    [tformPoints.T] % print tform
    if size(moving, 3) >= 16
        outputRef = imref3d(size(moving)); 
    else
        outputRef = imref2d(size(moving));
    end
    movingFinalImage = imwarp(moving, tformPoints, 'OutputView', outputRef);
    
    
    
    %% Get points in the moving image 
    viewPoints = true;
    [movingTform, ~, ~] = detectdots(movingFinalImage, movingThreshold, 'exons', viewPoints);
    movingTformGauss = getgaussian(movingTform.location, movingFinalImage);
    
    
    
    %% Calculate the new error between the points
    [~, ~, tformPointsMatch, ~] = matchpoints(fixedGauss, movingTformGauss, radius);
    tformPointsError = mean(tformPointsMatch.distance);
    tformPointsErrorX = mean(abs(tformPointsMatch.ref(:,1) - tformPointsMatch.points(:,1)));
    tformPointsErrorY = mean(abs(tformPointsMatch.ref(:,2) - tformPointsMatch.points(:,2)));
    tformPointsErrorZ = mean(abs(tformPointsMatch.ref(:,3) - tformPointsMatch.points(:,3)));
    fprintf('\nGaussian Points After Bead Alignment:\n');
    fprintf('Gaussian Point distance error: %.3f \n', tformPointsError);
    fprintf('X gaussian point error: %.3f \n', tformPointsErrorX);
    fprintf('Y gaussian point error: %.3f \n', tformPointsErrorY);
    fprintf('Z gaussian point errorn: %.3f \n', tformPointsErrorZ);
    
    %% Figure for all three points
    close all
    H1 = figure;
    imshow(max(fixed, [], 3), [0 3000]);
    hold on
    scatter(pointsRawMatch.ref(:,1), pointsRawMatch.ref(:,2), 150, 'gs', 'LineWidth', 1);
    scatter(pointsDapiMatch.ref(:,1), pointsDapiMatch.ref(:,2), 125, 'cx', 'LineWidth', 1);
    scatter(tformPointsMatch.ref(:,1), tformPointsMatch.ref(:,2), 150, 'm.', 'LineWidth', 1);
    legend('none','dapi','bead','Location','southeast');
    hold off
    save3Name = [pwd filesep 'fig-compare-transformations.fig'];
    savefig(H1, save3Name);

end