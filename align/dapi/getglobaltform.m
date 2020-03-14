function tform = getglobaltform(refPoints,matchPoints)
% function calculates the tforms between two images by comparing points
% found.
%
% Inputs:
% filename1 is the moving image and filename2 is reference image, and i is
% the barcode round
%
% Outputs: tforms and a control image for debugging
%
% Author: Sheel Shah
% Date: 2018
% Modified: Nico Pierson
% 4/2/2019


    indexPairs = knnsearch(matchPoints,refPoints,'NSMethod','Exhaustive');
    matchedPoints1 = matchPoints(indexPairs,:);
    matchedPoints2 = refPoints;
    dist = sum((matchedPoints1 - matchedPoints2).^2,2);
    dist = dist<70;
    matchedPoints1 = matchedPoints1(dist,:);
    matchedPoints2 = matchedPoints2(dist,:);
    tform = fitgeotrans(matchedPoints1(:,1:2),matchedPoints2(:,1:2),'affine'); % change affine to nonreflectivesimilarity
    tformv = fitgeotrans(matchedPoints1(:,2:3),matchedPoints2(:,2:3),'nonreflectivesimilarity');
    tformvv = fitgeotrans(matchedPoints1(:,[1 3]),matchedPoints2(:,[1 3]),'nonreflectivesimilarity');
   

    z = mean([tformv.T(3,2) tformvv.T(3,2)]);

    tformfinal = affine3d([tform.T(1,1) tform.T(1,2) 0 0;tform.T(2,1) tform.T(2,2) 0 0; 0 0 1 0; tform.T(3,1) tform.T(3,2) z 1]);
    tform = tformfinal;

end
