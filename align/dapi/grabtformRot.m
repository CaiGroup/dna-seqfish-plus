function tform = grabtformRot(moving, fixed, varargin)

im2 = max(moving,[],3);
im1 = max(fixed,[],3);

[optimizer,metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 500;
tform2d = imregtform(im2,im1,'rigid',optimizer,metric);
movingr= imwarp(moving,tform2d,'OutputView',imref2d(size(im1)));

maxim2 = permute(max(movingr,[],1),[2,3,1]);
maxim1 = permute(max(fixed,[],1),[2,3,1]);
maxim2(maxim2 == 2^16 - 1) = 150;
maxim1(maxim1 == 2^16 - 1) = 150;
c = normxcorr2(maxim2,maxim1);
[~, zpeak1] = find(c==max(c(:)));
zoffset1 = zpeak1-size(maxim2,2);

maxim2 = permute(max(movingr,[],2),[1,3,2]);
maxim1 = permute(max(fixed,[],2),[1,3,2]);
maxim2(maxim2 == 2^16 - 1) = 150;
maxim1(maxim1 == 2^16 - 1) = 150;
c = normxcorr2(maxim2,maxim1);
[~, zpeak2] = find(c==max(c(:)));
zoffset2 = zpeak2-size(maxim2,2);

zoffsetfinal = round((zoffset1+zoffset2)/2);

newT = [tform2d.T(1,1) tform2d.T(1,2)       0        0;...
        tform2d.T(2,1) tform2d.T(2,2)       0        0;...
             0              0               1        0;...
        tform2d.T(3,1) tform2d.T(3,2)  zoffsetfinal  1];
    
tform = affine3d(newT);
    

