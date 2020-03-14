function [dotlocations] = PointLocations(hyb, channum, points, consensuscell,copynumfinal,radius)
%#codegen

alllen = 0;
for j = 1:hyb
    for k = 1:channum
        alllen = alllen+length(consensuscell{j,k});
    end
end
allcodes = zeros(alllen,1);
allpoints = zeros(alllen,3);
ptchannel = allcodes;
hybnum = allcodes;
ptintabs = allcodes;
start = 1;
for j = 1:hyb
    for k = 1:channum
        if ~isempty(consensuscell{j,k})
            stop = length(consensuscell{j,k})+start-1;
            allcodes(start:stop) = consensuscell{j,k};
            allpoints(start:stop,1:3) = points{j}(k).channels;
            ptchannel(start:stop) =  k*ones(size(consensuscell{j,k},1),1);
            ptintabs(start:stop) = points{j}(k).intensity;
            hybnum(start:stop) =  j*ones(size(consensuscell{j,k},1),1);
            start = stop+1;
        end
    end
end
fun = unique(allcodes);
dotlocations = repmat({0},length(fun)-1,5);
for l = 2:length(fun) 
    dotlocations{fun(l),1} = allpoints(allcodes==fun(l),:); 
    dotlocations{fun(l),2} = ptchannel(allcodes==fun(l),:);
    dotlocations{fun(l),3} = hybnum(allcodes==fun(l),:);
    dotlocations{fun(l),5} = ptintabs(allcodes==fun(l),:);
    clusnum = copynumfinal(fun(l));
    if size(dotlocations{fun(l),1},1) < clusnum
        clusnum = size(dotlocations{fun(l),1},1);
    end
    bobobo = 1;
    keepgoing = 1;
    sprev = Inf;
    C=1;
    while keepgoing == 1
        if clusnum <2
            C = mean(dotlocations{fun(l),1},1);
            idx = ones(size(dotlocations{fun(l),1},1),1);
        else
            [idx,C] = kmeans(dotlocations{fun(l),1},clusnum);
        end
        d = zeros(size(C,1),1);
        for m = 1:size(C,1) 
            d(m) = sum(pdist2(C(m,:),dotlocations{fun(l),1}(idx == m,:))>sqrt(radius)+.0001)>0;
        end
        s = sum(d);
        if s > 0 && sprev > s
            clusnum = clusnum + s;
            sprev = s;
        elseif pdist(C)<sqrt(radius)+.00001
            clusnum = clusnum - 1;
        else
            keepgoing = 0;
        end
        bobobo = bobobo + 1;
        if bobobo == 100
            keepgoing = 0;
        end
    end
    dotlocations{fun(l),4} = C;
end