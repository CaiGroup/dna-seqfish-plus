function [consensuscell, copynumfinal ] = BarcodeNoMiji_v8( channels, points, hyb, barcodekey,radius,alloweddiff)
%#codegen

%channels = length(channels);

%find all dots in all images of all hybs
alloweddiff = alloweddiff + 1;

% points = cell(1,hyb);
% for i = 1:hyb
%     [~,points{i}] = findDotsBarcodeV2(hybnum(i).color, multiplier(i,:), HCRorFISH);
% end

%Initializing
channelcell = cell(hyb,channels);
idxcell = cell(hyb,channels);
for i = 1:hyb
    for j = 1:channels
        c = size(points{i}(j).channels,1);
        filler = repmat({0},c,hyb);
        for k = 1:c; filler{k,i} = j;end 
        idxfill = repmat({zeros(1,channels)},c,hyb);
        for k = 1:c; makevec = zeros(1,channels); makevec(j) = k; idxfill{k,i} = makevec; end
        channelcell{i,j} = filler;
        idxcell{i,j} = idxfill;
    end
end

%colocalize channels and output found codes and dots indices
for i = 1:hyb
    for k = 1:hyb
        if k ~= i
            for j = 1:channels
                calledall = cell(1,channels);
                calledidx = cell(1,channels);
                for l = 1:channels
                    [~, ~, called1, ~, pair ] = colocalizeBarcodeV3(points, i,k, [j l],radius);
                    calledall{l} = l*double(called1);
                    brat = zeros(length(called1),1);
                    if sum(pair) ~= 0
                        brat(pair(:,1)) = pair(:,2);
                    end
                    calledidx{l} = brat;
                end
                calledall = [calledall{:}];
                calledidx = [calledidx{:}];
                for iter = 1:size(calledall,1)
                    xiter = calledall(iter,:)~=0;
                    if sum(xiter) > 0
                        channelcell{i,j}{iter,k} = calledall(iter,xiter);
                    end
                end
                for iter = 1:size(calledidx,1)
                    idxcell{i,j}{iter,k} = calledidx(iter,:);
                end
            end
        end
    end
end

%call Barcodes
calledcell = cell(hyb,channels);
for k = 1:hyb
    for j = 1:channels
        lentemp = size(channelcell{k,j},1);
        br = zeros(lentemp,1);
        len = zeros(lentemp,1);
        for i = 1:lentemp
            coder.varsize('bobo',[1 Inf], [0 1]);
            bobo = zeros(1,0);
            for bobiter = 1:hyb
                bobo = [bobo channelcell{k,j}{i,bobiter}];
            end
            br(i,1) = sum(bobo==0,2);
            drop = br(i,1) > alloweddiff -1;
            if drop ~= 1
                br2 = bobo;
                br2(br2 == 0) = [];
                len(i,1) = length(br2);
            else 
                for inter = 1:hyb
                    channelcell{k,j}{i,inter} =0;
                    idxcell{k,j}{i,inter} =zeros(1,channels);
                end
                len(i,1) = 0;
            end
        end
        
        %reduce multiple matches and fill in drops
        multi = len+br>hyb | len > hyb-alloweddiff+1;
        [rows] = find(multi ==1);
        for l = 1:length(rows)
            dummycell = cell(1,hyb);
            for cis = 1:hyb
                dummycell{cis} = channelcell{k,j}{rows(l),cis}(1,:);
            end
            possbarcodes = allcomb(dummycell);
            Dposs = pdist2(possbarcodes,barcodekey,'hamming');
            A = Dposs == min(min(Dposs,[],1),[],2);
            if size(A,1)==1
                rvec = A(1,:);
                [code, gene] = find(rvec == 1);
                g1 = sum(Dposs(rvec)*hyb < alloweddiff);
            else
                [code, gene] = find(A == 1);
                g1 = sum(Dposs(A)*hyb < alloweddiff);
            end
            if length(code) == 1 && g1 > 0
                for iter = 1:hyb; channelcell{k,j}{rows(l),iter} = barcodekey(gene(1),iter); end
            elseif g1 > 0 && length(code) > 1
                possreal = possbarcodes(code,:);
                ind = cell(1,hyb);
                for cis = 1:hyb
                    ind{cis} = idxcell{k,j}{rows(l),cis};
                end
                vari = zeros(size(possreal,1),1);
                for p = 1:size(possreal,1)
                    set = nan(hyb,3);
                    for h = 1:hyb
                        if possreal(p,h) > 0
                            set(h,:) = points{h}(possreal(p,h)).channels(ind{h}(possreal(p,h)),:);
                        end
                    end
                    vari(p) = sum(var(set,'omitnan'));
                end
                vernoi = vari == min(vari);
                if sum(vernoi) == 1
                    possreal = possreal(vernoi',:);
                    for iter = 1:hyb; channelcell{k,j}{rows(l),iter} = possreal(1,iter); end
                else
                    for inter = 1:hyb
                        channelcell{k,j}{rows(l),inter} =0;
                        idxcell{k,j}{rows(l),inter} =zeros(1,channels);
                    end
                end
            else
                for inter = 1:hyb
                    channelcell{k,j}{rows(l),inter} =0;
                    idxcell{k,j}{rows(l),inter} =zeros(1,channels);
                end
            end
        end
        
        % call barcodes
        if ~isempty(channelcell{k,j})
            %call dots
            diga = zeros(size(channelcell{k,j}));
            for dumrow = 1:size(channelcell{k,j},1)
                for dumdum = 1:hyb 
                    diga(dumrow,dumdum) = channelcell{k,j}{dumrow,dumdum}(1); 
                end
            end
            D = pdist2(diga,barcodekey,'hamming');
            called = zeros(size(D,1),1);
            hordor = find(min(D,[],2)<=((alloweddiff-1)/hyb));
            for iter = 1:length(hordor)
                [~,called(hordor(iter))] = min(D(hordor(iter),:));
                for ider = 1:hyb
                    mu = zeros(1,channels); 
                    mu(barcodekey(called(hordor(iter)),ider))=1;
                    idxcell{k,j}{hordor(iter),ider}=idxcell{k,j}{hordor(iter),ider}.*mu;
                end
            end
            calledcell{k,j} = called;
        else
            idxcell{k,j} = {0};
            calledcell{k,j} = zeros(size(channelcell{k,j},1),1);
        end        
    end
end

%consensus point calling
compiledcell = cell(hyb,channels);
for i = 1:hyb
     for j = 1:channels
         compiledcell{i,j} = zeros(size(calledcell{i,j},1),hyb);
    end
end

for i = 1:hyb
    for j = 1:channels
        calledrows = find(calledcell{i,j}>0);
        if ~isempty(calledrows)
            for k = 1:length(calledrows)
                for l = 1:hyb
                    [~,c,v] = find(idxcell{i,j}{calledrows(k),l});
                    if ~isempty(c)
                        compiledcell{l,c(1)}(v(1),i) = calledcell{i,j}(calledrows(k));
                    end
                end
            end
        end
    end
end

%Dropping Ambiguous Matches
consensuscell = cell(hyb, channels);
for i = 1:hyb
    for j = 1:channels
        % Drop Ambigous Matches
        M = zeros(size(compiledcell{i,j},1),1);
        eq = M;
        sums = find(sum(compiledcell{i,j},2)>0);
        for k = 1:size(sums,1)
            hordor = compiledcell{i,j}(sums(k),:);
            hordor(compiledcell{i,j}(sums(k),:)==0)=[];
            if isscalar(unique(hordor)) == 1
                M(sums(k)) = unique(hordor);
                eq(sums(k)) = 1;
            else 
                [M(sums(k))] = mode(hordor);
                [a]=hist(hordor,unique(hordor));
                eq(sums(k)) = sum(ismember(a,max(a)));
            end
        end
        M(eq>1) = 0;
        consensuscell{i,j} = M;
    end
end

copynum = zeros(size(barcodekey,1),hyb);

for j = 1:hyb    
    alllen = 0;
    for k = 1:channels
    	alllen = alllen+length(consensuscell{j,k});
    end
    allcalled = zeros(1,alllen);
    start = 1;
    for k = 1:channels
        stop = length(consensuscell{j,k})+start-1;
        allcalled(start:stop) = consensuscell{j,k};
        start = stop+1;
    end
    copy = histc(allcalled,0:size(barcodekey,1));
    copynum(:,j) = copy(2:end);
end

copynumfinal = zeros(size(barcodekey,1),1);
for i = 1:size(barcodekey,1)
    copynumfinal(i) = max(copynum(i,:));
end