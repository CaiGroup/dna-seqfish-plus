function [m,dots] = findDotsBarcodeV3(images, multiplier, HCRorFISH,debug)
%#codegen

m = cell(1,length(images));
s = struct('channels',[0 0 0],'intensity',0);
dots = repmat(s,length(images),1);
for dee = 1:length(images)
    
    fish = images{dee};
    
    if HCRorFISH(dee) == 1
        h = fspecial3('log',7,1);
        BW = imfilter(double(fish), -20*h, 'replicate');
        msk = true(3,3,3);
        msk(2,2,2) = false;
        apply = BW < multiplier(dee);
        BW(apply) = 0;
        s_dil = imdilate(BW,msk);
        m{dee} = BW > s_dil; 
    else
        msk = true(3,3,3);
        msk(2,2,2) = false;
        apply = fish < multiplier(dee);
        fish(apply) = 0;
        s_dil = imdilate(fish,msk);
        m{dee} = fish > s_dil; 
    end   

    
    if debug == 1
        figure 
        imshow(max(fish,[],3),[multiplier(dee)-150 multiplier(dee)]);
        hold on;
        [v2,v1]=find(max(m{dee},[],3)==1);
        scatter(v1(:),v2(:),75);
        hold off;
        %txt_h = labelpoints(v1+.05, v2+.05, ind2sub(size(v1),v1), 'NE',.01 ,'Color','y');
    end
    
    [y,x,z] = ind2sub(size(m{dee}),find(m{dee} == 1));
    
    dots(dee).channels = [x y z];
    im = max(images{dee},[],3);
    for i = 1:length(y)
        dots(dee).intensity(i,1) = im(y(i),x(i));
    end
end