function [fixed, gauss, radial] = imerrorgaussianradial(image, refPoints, fixedPoints, hyb, ch, textFilePath, searchradius)
% function compare super resolved points gaussian and radial center
%
% Date: 10/2/2019

    
    
    %% Gaussian Adjust the dots using the points and the raw images
    gauss = getgaussian(fixedPoints, image);
    radial = getradialcenter(fixedPoints, image);
    
    
    %% set up file to write
    fileID = fopen(textFilePath,'a+');
    if hyb == 1
        fprintf(fileID, '\nChannel %.0f: \n', ch);
    end
    fprintf(fileID, '\nHyb %.0f: \n', hyb);
    
    
    %% Colocalization rate
    numRefPoints = size(refPoints, 1);
    
    
    %% Calculate and print the error for fixed Points
    [~, ~, pointsFixed, ~] = matchpoints(fixedPoints, refPoints, searchradius);
    % Print the Error
    rawError = mean(pointsFixed.distance);
    rawErrorX = mean(abs(pointsFixed.ref(:,1) - pointsFixed.points(:,1)));
    rawErrorY = mean(abs(pointsFixed.ref(:,2) - pointsFixed.points(:,2)));
    rawErrorZ = mean(abs(pointsFixed.ref(:,3) - pointsFixed.points(:,3)));
    numRawPoints = size(pointsFixed.ref,1);
    rawCRate = numRawPoints / numRefPoints;
    fprintf(fileID, '\nPoints :\n');
    fprintf(fileID, 'Colocalization Rate: %.3f \n', rawCRate);
    fprintf(fileID, 'Point distance error: %.3f \n', rawError);
    fprintf(fileID, 'X point error: %.3f \n', rawErrorX);
    fprintf(fileID, 'Y point error: %.3f \n', rawErrorY);
    fprintf(fileID, 'Z point error: %.3f \n', rawErrorZ);
    
    
    
    %% Calculate and print the error for gauss fitted points
    [~, ~, pointsGauss, ~] = matchpoints(gauss, refPoints, searchradius);
    % Print the Error
    gaussError = mean(pointsGauss.distance);
    gaussErrorX = mean(abs(pointsGauss.ref(:,1) - pointsGauss.points(:,1)));
    gaussErrorY = mean(abs(pointsGauss.ref(:,2) - pointsGauss.points(:,2)));
    gaussErrorZ = mean(abs(pointsGauss.ref(:,3) - pointsGauss.points(:,3)));
    numGaussPoints = size(pointsGauss.ref,1);
    gaussCRate = numGaussPoints / numRefPoints;
    fprintf(fileID, '\nGaussian Points:\n');
    fprintf(fileID, 'Colocalization Rate: %.3f \n', gaussCRate);
    fprintf(fileID, 'Gaussian point distance error: %.3f \n', gaussError);
    fprintf(fileID, 'X gaussian point error: %.3f \n', gaussErrorX);
    fprintf(fileID, 'Y gaussian point error: %.3f \n', gaussErrorY);
    fprintf(fileID, 'Z gaussian point error: %.3f \n', gaussErrorZ);
    
    
    
    %% Calculate and print the error for radial center points
    [~, ~, pointsRadial, ~] = matchpoints(radial, refPoints, searchradius);
    radError = mean(pointsRadial.distance);
    radErrorX = mean(abs(pointsRadial.ref(:,1) - pointsRadial.points(:,1)));
    radErrorY = mean(abs(pointsRadial.ref(:,2) - pointsRadial.points(:,2)));
    radErrorZ = mean(abs(pointsRadial.ref(:,3) - pointsRadial.points(:,3)));
    numRadialPoints = size(pointsRadial.ref,1);
    radialCRate = numRadialPoints / numRefPoints;
    fprintf(fileID, '\nRadial Center Points:\n');
    fprintf(fileID, 'Colocalization Rate: %.3f \n', radialCRate);
    fprintf(fileID, 'Radial Center point distance error: %.3f \n', radError);
    fprintf(fileID, 'X radial Center point error: %.3f \n', radErrorX);
    fprintf(fileID, 'Y radial Center point error: %.3f \n', radErrorY);
    fprintf(fileID, 'Z radial Center point error: %.3f \n', radErrorZ);
    
    fclose(fileID);
    
    %% assign variables
    %fixed = pointsFixed.points;
    %gauss = pointsGauss.points;
    %radial = pointsRadial.points;
    fixed = fixedPoints;
    
end