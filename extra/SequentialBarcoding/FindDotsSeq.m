function [color, m, dots, Data] = FindDotsSeq(roiPath, hybs, channels, multiplier, segmentation, debug, hybnum)
% FindDotsSeq finds the dots for each set of hyb rounds and save the copy
% number and the positions of the points as Data.
%
% Author: Sheel and Yodai
% Email: ytakei@caltech.edu
% Date: 2018
% Modified By: Nico Pierson
% Email: nicopgt@caltech.edu
% Date: 1/29/2019


for num = 1:hybs
    
    % Have to put in piece of code to adjust for extra hybnum in the
    % beginning - take out after adjusting function
    color = hybnum(num).color;

    for dee = 1:length(channels)
        % Set Threshold
        threshold = multiplier(num, dee);
        
        fish = color{channels(dee)}(:,:,1:end); % get image
        logFish=[];
        fish=double(fish);
        for i=1:size(fish,3)
            logFish(:,:,i)=logMask(fish(:,:,i)); % use laplacian filter
        end
        cands=imregionalmax(logFish);
        m{dee}= cands & logFish > threshold; % select dots with regional maxima and those in the laplacian filter that are above the threshold


        if debug == 1
            figure 
            imshow(max(color{channels(dee)},[],3),[min(min(max(color{channels(dee)},[],3))) mean(mean(max(color{channels(dee)},[],3)))+5000]);
            hold on;
            [v2,v1]=find(max(m{dee},[],3)==1);
            scatter(v1(:),v2(:),75);
            hold off;
            %txt_h = labelpoints(v1+.05, v2+.05, ind2sub(size(v1),v1), 'NE',.01 ,'Color','y');
        end

        [y,x,z] = ind2sub(size(m{dee}),find(m{dee} == 1));

        dots(dee).channels = [x y z];
        im = max(color{channels(dee)},[],3);
        for i = 1:length(y)
            dots(dee).intensity(i,1) = im(y(i),x(i));
        end

    end

    if strcmp(segmentation, 'roi')
        vertex = selfseg(roiPath);
        for i = 1:length(vertex)
                for k = 1:length(channels)
                    include = inpolygon(dots(k).channels(:,1),dots(k).channels(:,2),vertex(i).x,vertex(i).y);
                    copy(k,i) = sum(include);
                    celldata.Positions{k,i} = dots(k).channels(include,:);
                    if isfield(dots,'intensity')
                        celldata.Intensity{k,i} = dots(k).intensity(include);
                    end
                end
        end
    else
        copy = [];
        celldata = [];
    end
    Data(num).copy = copy;
    Data(num).celldata = celldata;
end