function points = detectintrons(image)
% Get Intron dots and fiduciary dots
%
% Author: Sheel Shah
% Date: March 2018
% Modified By: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 3/29/2019

%{
    logFish = max(image,[],3);
    thresh = multithresh(logFish,2);
    msk = true(3,3,3);
    msk(2,2,2) = false;
    % assign, to every voxel, the maximum of its neighbors
    logFish = image;
    apply = logFish < thresh(2); %thresh(2) - 1000 % minus 1000 to pick up slightly more points
    logFish(apply) = 0;
    s_dil = imdilate(logFish,msk);
    m = logFish > s_dil; % M is 1 wherever a voxel's value is greater than its neighbors
    %}
    
    % See if this new method of thresholding can grab introns
    thresh = adaptthresh(image,'neigh',[3 3 3],'Fore','bright','Statistic','median');
    highthresh = thresh >= 0.005; % 0.02 for without background subtract and 0.005 with
    imageThresh = thresh .* highthresh;
    m = imregionalmax(imageThresh);


    %% Remove border dots
    bordSize = 5;
    bord = ones(size(m));
    bord(1:bordSize,:,:) = 0;
    bord(end-bordSize:end,:,:) = 0;
    bord(:,end-bordSize:end,:) = 0;
    bord(:,1:bordSize,:) = 0;
    m = m.*logical(bord);
    
    [y,x,z] = ind2sub(size(m),find(m == 1));
            figure 
            imshow(max(image,[],3), 'InitialMagnification', 'fit','DisplayRange', [0 1000]);
            hold on;
            [v2,v1]=find(max(m,[],3)==1);
            scatter(v1(:),v2(:),75);
            hold off;

     points = [x, y, z];
end