function [mTest, dotsTest, DataTest] = checkthreshold(PathNameRoi, image, threshold, segmentation, position, channels, debug)
    %load([PathName '\Pos' num2str(posnum) '\Pos' num2str(posnum) 'SequentialHybsNewRegis.mat']);


     %% Find Dots
     %{
     logFish=[];
     fish=double(image);
            %fish = double(max(fish,[],3));
     for i=1:size(fish,3)
         logFish(:,:,i)=logMask(fish(:,:,i));
     end
     cands=imregionalmax(logFish);
     thresh=threshold;
     mTest = cands & logFish > thresh;   
     %}
     
    %% Try other part of code
    %{
    fish=double(image);
    thresh = multithresh(max(fish,[],3),1); 
    % * threshold
    if thresh == 0
        thresh = mean(mean(double(max(fish,[],3))));
    end
    msk = true(3,3,3);
    msk(2,2,2) = false;
    % assign, to every voxel, the maximum of its neighbors
    apply = fish < threshold;
    fish(apply) = 0;
    s_dil = imdilate(fish,msk);
    mTest = fish > s_dil;
    %}
    
    %% Try third option to find dots
         logFish=[];
     fish=double(image);
            %fish = double(max(fish,[],3));
     for i=1:size(fish,3)
         logFish(:,:,i)=logMask(fish(:,:,i));
     end
            cands=imregionalmax(logFish);
        sortedValues = unique(logFish(cands));
        p10 = round(length(sortedValues)*.1);
        maxValues = sortedValues(end-p10:end);  
        maxIndex = ismember(logFish,maxValues);  
        logFish2 = logFish;
        logFish2(maxIndex) = mean(mean(mean(logFish)));
        %thresh=multithresh(mat2gray(logFish2(cands)),2)*max(logFish2(cands))*multiplier(dee);
        mTest= cands & logFish > threshold;

    %% Show Images
    if debug == 1
        figure 
        imshow(max(fish,[],3),[min(min(max(fish,[],3))) mean(mean(max(fish,[],3)))+5000]);
        hold on;
        [v2,v1]=find(max(mTest,[],3)==1);
        scatter(v1(:),v2(:),75);
        hold off;
        %txt_h = labelpoints(v1+.05, v2+.05, ind2sub(size(v1),v1), 'NE',.01 ,'Color','y');
    end

    [y,x,z] = ind2sub(size(mTest),find(mTest == 1));

    dotsTest.channels = [x y z];
    im = max(fish,[],3);
    for i = 1:length(y) % Get the intensity
        dotsTest.intensity(i,1) = im(y(i),x(i));
    end

    if strcmp(segmentation, 'roi')
        vertex = selfseg([PathNameRoi '\Pos' num2str(position) '\RoiSet']);
        for i = 1:length(vertex)
            %for k = 1:length(channels)
            if i == 1
                include1 = inpolygon(dotsTest.channels(:,1),dotsTest.channels(:,2),vertex(i).x,vertex(i).y);
            end
                include = inpolygon(dotsTest.channels(:,1),dotsTest.channels(:,2),vertex(i).x,vertex(i).y);
                fprintf('Number of Dots in %.0f ROI: %.0f\n', i, sum(include));
            %end
        end
        %{
        for i = 1:length(vertex)
                for k = 1:length(channels)
                    include = inpolygon(dots(k).channels(:,1),dots(k).channels(:,2),vertex(i).x+regvec(1),vertex(i).y+regvec(2));
                    copy(k,i) = sum(include);
                    celldata.Positions{k,i} = dots(k).channels(include,:);
                    if isfield(dots,'intensity')
                        celldata.Intensity{k,i} = dots(k).intensity(include);
                    end
                end
        end
        %}
    end
    
    if debug == 1
        figure 
        imshow(max(fish,[],3),[min(min(max(fish,[],3))) mean(mean(max(fish,[],3)))+5000]);
        hold on;
        [v2,v1]=find(max(mTest,[],3)==1);
        scatter(x(include1),y(include1),75);
        hold off;
        %txt_h = labelpoints(v1+.05, v2+.05, ind2sub(size(v1),v1), 'NE',.01 ,'Color','y');
    end
    
    
    close all;
    
    DataTest(num).copy = copy;
    DataTest(num).celldata = celldata;
end