function [tform, controlimage] = findtformV3(filename1,filename2,i)

%first use imageJ to do a 10 pixel background subtraction
%save individual channels
%Filename2 is reference image
close all;
img = loadtiff(filename1);
img2 = loadtiff(filename2);

%logFish = img;
logFish = max(img,[],3);
%logFish = double(max(fish,[],3));
thresh = multithresh(logFish,2);
%# s = 3D array
msk = true(3,3,3);
msk(2,2,2) = false;
%# assign, to every voxel, the maximum of its neighbors
logFish = img;
apply = logFish < thresh(2);
logFish(apply) = 0;
s_dil = imdilate(logFish,msk);
m = logFish > s_dil; %# M is 1 wherever a voxel's value is greater than its neighbors

[y,x,z] = ind2sub(size(m),find(m == 1));
        figure 
        imshow(max(img,[],3));
        hold on;
        [v2,v1]=find(max(m,[],3)==1);
        scatter(v1(:),v2(:),75);
        hold off;
centroids = [x,y,z];
logFish = max(img2,[],3);
%logFish = double(max(fish,[],3));
thresh = multithresh(logFish,2);
%# s = 3D array
msk = true(3,3,3);
msk(2,2,2) = false;
%# assign, to every voxel, the maximum of its neighbors
logFish = img2;
apply = logFish < thresh(2);
logFish(apply) = 0;
s_dil = imdilate(logFish,msk);
m = logFish > s_dil; %# M is 1 wherever a voxel's value is greater than its neighbors

[y,x,z] = ind2sub(size(m),find(m == 1));

centroids2 = [x,y,z];

        figure 
        imshow(max(img2,[],3));
        hold on;
        [v2,v1]=find(max(m,[],3)==1);
        scatter(v1(:),v2(:),75);
        hold off;

    indexPairs = knnsearch(centroids,centroids2,'NSMethod','Exhaustive');
    matchedPoints1 = centroids(indexPairs,:);
    matchedPoints2 = centroids2;
    dist = sum((matchedPoints1 - matchedPoints2).^2,2);
    dist = dist<70;
    matchedPoints1 = matchedPoints1(dist,:);
    matchedPoints2 = matchedPoints2(dist,:);
    %moving_pts_adj= cpcorr(matchedPoints1(:,2), matchedPoints2(:,2), max(img,[],3), max(img2,[],3));
    %close all;
    figure; showMatchedFeatures(max(img,[],3), max(img2,[],3), matchedPoints1(:,1:2), matchedPoints2(:,1:2));
    tform = fitgeotrans(matchedPoints1(:,1:2),matchedPoints2(:,1:2),'affine');
    tformv = fitgeotrans(matchedPoints1(:,2:3),matchedPoints2(:,2:3),'nonreflectivesimilarity');
    tformvv = fitgeotrans(matchedPoints1(:,[1 3]),matchedPoints2(:,[1 3]),'nonreflectivesimilarity');
    
%    img3 = imwarp(img,tform,'OutputView',imref2d(size(img2)));
    
% yz = permute(max(img3,[],2),[1 3 2]);
% yzREF = permute(max(img2,[],2),[1 3 2]);
% c = normxcorr2(yz,yzREF); 
% [~, xpeak] = find(c==max(c(:)));
% z1 = xpeak-size(yz,2);
% 
% xz = permute(max(img3,[],1),[2 3 1]);
% xzREF = permute(max(img2,[],1),[2 3 1]);
% c = normxcorr2(xz,xzREF); 
% [~, xpeak] = find(c==max(c(:)));
% z2 = xpeak-size(xz,2);

z = mean([tformv.T(3,2) tformvv.T(3,2)]);

tformfinal = affine3d([tform.T(1,1) tform.T(1,2) 0 0;tform.T(2,1) tform.T(2,2) 0 0; 0 0 1 0; tform.T(3,1) tform.T(3,2) z 1]);
tform = tformfinal;
controlimage = imwarp(img,tform,'OutputView',imref3d(size(img2)));

falsecolorOverlay1 = imfuse(max(img2,[],3),max(controlimage,[],3));
falsecolorOverlay2 = imfuse(max(img2,[],3),max(img,[],3));
figure;
imshow(falsecolorOverlay1,'InitialMagnification','fit');
figure;
imshow(falsecolorOverlay2,'InitialMagnification','fit');

LinkFigures(4:5)
xz = permute(max(img,[],1),[3 2 1]);
xzREF = permute(max(img2,[],1),[3 2 1]);
xzshow = repmat(xz,100,1);
xzREFshow = repmat(xzREF,100,1);
falsecolorOverlay1 = imfuse(xzREFshow,repmat(permute(max(controlimage,[],1),[3 2 1]),100,1));
falsecolorOverlay2 = imfuse(xzREFshow,xzshow);
figure;
imshow(falsecolorOverlay1,'InitialMagnification','fit');
figure;
imshow(falsecolorOverlay2,'InitialMagnification','fit');

LinkFigures(6:7)
    C = strsplit(filename1, filesep);
    FileName = C{end};
    k = strfind(filename1, filesep);
    PathName = filename1(1:k(end)-1);
saveastiff(controlimage, [PathName filesep 'channel' i 'Registered.tif'])
