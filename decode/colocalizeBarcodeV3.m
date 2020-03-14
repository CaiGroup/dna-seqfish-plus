function [ pts1, pts2, called1, called2, pair ] = colocalizeBarcodeV3(points, hybnum,hybnum2, channels,radius)
%#codegen

    conc1 = points{hybnum}(channels(1)).channels;
    conc2 = points{hybnum2}(channels(2)).channels;

[idx,D]= rangesearch(conc1,conc2,sqrt(radius)+.00001);
[idx2,D2] = rangesearch(conc2,conc1,sqrt(radius)+.00001);

%coder.varsize('pair', [Inf 2], [1 0]);
pair = zeros(max([size(conc1,1) size(conc2,1)]),2);
%pair = [0 0];

if ~isempty(idx) && ~isempty(idx2)
    %Match Unique Points symmetrically
    idx2temp = idx2;
    D2temp = D2;
    idxtemp = idx;
    Dtemp = D;
    %rows = find(idx1==1);
    for i = 1:length(idx)
        if length(idx{i}) == 1
            if length(idx2{idx{i}(1)}) == 1
                idxtemp{i,1} = idx{i};
                Dtemp{i,1} = D{i};
            else
                dummy = idx2{idx{i}(1)};
                dummyD = D2{idx{i}(1)};
                %matches = cellfun(@length,idx(dummy))';
                matches = zeros(1,length(dummy)); 
                for iter = 1:length(dummy) 
                    matches(iter) = length(idx{dummy(iter)}); 
                end
                keepers = matches == 1;
                idxtemp{i,1} = idx{i};
                Dtemp{i,1} = D{i};
                idx2temp{idx{i}} = dummy(keepers);
                D2temp{idx{i}} = dummyD(keepers);
            end
        end
    end
    
    idxunique = idxtemp;
    Dunique = Dtemp;
    idx2unique = idx2temp;
    D2unique = D2temp;
    for i = 1:length(idx2temp)
        if length(idx2temp{i}) == 1
            if length(idxtemp{idx2temp{i}(1)}) == 1
                idx2unique{i,1} = idx2temp{i};
                D2unique{i,1} = D2temp{i};
            else
                dummy = idxtemp{idx2temp{i}(1)};
                dummyD = Dtemp{idx2temp{i}(1)};
                matches = zeros(1,length(dummy)); 
                for iter = 1:length(dummy) 
                    matches(iter) = length(idx2temp{dummy(iter)}); 
                end
                keepers = matches == 1;
                idx2unique{i,1} = idx2temp{i};
                D2unique{i,1} = D2temp{i};
                idxunique{idx2temp{i}} = dummy(keepers);
                Dunique{idx2temp{i}} = dummyD(keepers);
            end
        end
    end

    % Place unique matches in curated points list
    % Remove points from not uniquely matched list
    idxcurated = repmat({0},length(idxunique),1);
    idx2curated = repmat({0},length(idx2unique),1);
    counter = 1;
    if length(idx2unique)>length(idxunique)  
        for i = 1:length(idx2unique)
            if length(idx2unique{i})==1
                if logical(sum(idxunique{idx2unique{i}(1)} == i)) && length(idxunique{idx2unique{i}(1)}) == 1
                    idx2curated{i} = idx2unique{i};
                    idxcurated{idx2unique{i}(1)} = idxunique{idx2unique{i}(1)};
                    D2unique{i} = zeros(1,0);
                    Dunique{idx2unique{i}} = zeros(1,0);
                    idxunique{idx2unique{i}} = zeros(1,0);
                    idx2unique{i} = zeros(1,0);
                    pair(counter,2) = idx2curated{i}(1);
                    pair(counter,1) = i;
                    counter = counter +1;
                end
            end
        end
    else
        for i = 1:length(idxunique)
            if length(idxunique{i})==1 
                if logical(sum(idx2unique{idxunique{i}(1)} == i)) && length(idx2unique{idxunique{i}(1)}) == 1
                    idxcurated{i} = idxunique{i};
                    idx2curated{idxunique{i}(1)} = idx2unique{idxunique{i}(1)};
                    Dunique{i} = zeros(1,0);
                    D2unique{idxunique{i}} = zeros(1,0);
                    idx2unique{idxunique{i}} = zeros(1,0);
                    idxunique{i} = zeros(1,0);
                    pair(counter,1) = idxcurated{i}(1); %conc1
                    pair(counter,2) = i; %conc2
                    counter = counter + 1;
                end
            end
        end
    end

    %Distance Minimization

    idxmin = idxunique;
    Dmin = Dunique;
    for i = 1:length(idxunique) 
        if length(Dunique{i}) > 1
            smallest = Dunique{i} == min(Dunique{i});
            if sum(smallest) == 1
                idxmin{i,1} = idxunique{i}(smallest);
                Dmin{i,1} = Dunique{i}(smallest);
            end
        end
    end

    idx2min = idx2unique;
    D2min = D2unique;
    for i = 1:length(idx2unique) 
        if length(D2unique{i}) > 1
            smallest = D2unique{i} == min(D2unique{i});
            if sum(smallest) == 1
                idx2min{i,1} = idx2unique{i}(smallest);
                D2min{i,1} = D2unique{i}(smallest);
            end
        end
    end

    if length(idx2min)>length(idxmin)
        for i = 1:length(idx2min)
            if length(idx2min{i})==1 && length(idxmin{idx2min{i}(1)}) == 1
                if idxmin{idx2min{i}(1)} == i
                    idx2curated{i} = idx2min{i};
                    idxcurated{idx2min{i}(1)} = idxmin{idx2min{i}(1)};
                    D2min{i} = zeros(1,0);
                    Dmin{idx2min{i}} = zeros(1,0);
                    idxmin{idx2min{i}} = zeros(1,0);
                    idx2min{i} = zeros(1,0);
                    pair(counter,2) = idx2curated{i}(1);
                    pair(counter,1) = i;
                    counter = counter +1;
                else
                    idx2min{i} = zeros(1,0);
                    D2min{i} = zeros(1,0);
                end
            elseif length(idx2min{i})==1 && isempty(idxmin{idx2min{i}(1)})
                idx2min{i} = zeros(1,0);
                D2min{i} = zeros(1,0);
            end
        end
    else
        for i = 1:length(idxmin)
             if length(idxmin{i})==1 && length(idx2min{idxmin{i}(1)}) == 1
                if idx2min{idxmin{i}(1)} == i
                    idxcurated{i} = idxmin{i};
                    idx2curated{idxmin{i}(1)} = idx2min{idxmin{i}(1)};
                    Dmin{i} = zeros(1,0);
                    D2min{idxmin{i}} = zeros(1,0);
                    idx2min{idxmin{i}} = zeros(1,0);
                    idxmin{i} = zeros(1,0);
                    pair(counter,2) = i;
                    pair(counter,1) = idxcurated{i}(1);
                    counter = counter +1;
                else
                    idxmin{i} = zeros(1,0);
                    Dmin{i} = zeros(1,0);
                end
             elseif length(idxmin{i})==1 && isempty(idx2min{idxmin{i}(1)})
                idxmin{i} = zeros(1,0);
                Dmin{i} = zeros(1,0);
            end
        end
    end


    idxz = idxmin;
    for i = 1:length(idxz) 
        if length(idxz{i}) > 1 
            zref = conc2(i,3);
            zdots = conc1(idxmin{i}, 3);
            keepers = zdots == zref;
            dummy = idxz{i};
            idxz{i} = dummy(keepers);
        end 
    end

    idx2z = idx2min;
    for i = 1:length(idx2z) 
        if length(idx2z{i}) > 1 
            zref = conc1(i,3);
            zdots = conc2(idx2min{i}, 3);
            keepers = zdots == zref;
            dummy = idx2min{i};
            idx2z{i} = dummy(keepers);
        end 
    end
 
    if length(idx2z)>length(idxz)
        for i = 1:length(idx2z)
            if length(idx2z{i})==1 &&  length(idxz{idx2z{i}(1)}) == 1
                if idxz{idx2z{i}(1)} == i
                    idx2curated{i} = idx2z{i};
                    idxcurated{idx2z{i}(1)} = idxz{idx2z{i}(1)};
                    idxz{idx2z{i}} = zeros(1,0);
                    idx2z{i} = zeros(1,0);
                    pair(counter,2) = idx2curated{i}(1);
                    pair(counter,1) = i;
                    counter = counter +1;
                else
                    idx2z{i} = zeros(1,0);
                end
            elseif length(idx2z{i})==1 && isempty(idxz{idx2z{i}(1)})
                idx2z{i} = zeros(1,0);
            elseif length(idx2z{i})>1
                %len = cellfun(@length,idxz(idx2z{i}));
                len = zeros(length(idx2z{i}),1);
                for iter = 1:length(idx2z{i})
                    len(iter) = length(idxz{idx2z{i}(iter)});
                end
                if sum(len) == 1
                    keepers = len' == 1;
                    vals = idx2z{i}(keepers);
                    idx2curated{i} = vals;
                    idxcurated{vals} = idxz{vals(1)};
                    idx2z{i} = zeros(1,0);
                    idxz{vals} = zeros(1,0);
                    pair(counter,2) = idx2curated{i}(1);
                    pair(counter,1) = i;
                    counter = counter +1;
                end     
            end
        end
    else
        for i = 1:length(idxz)
             if length(idxz{i})==1 && length(idx2z{idxz{i}(1)}) == 1
                if idx2z{idxz{i}(1)} == i
                    idxcurated{i} = idxz{i};
                    idx2curated{idxz{i}(1)} = idx2z{idxz{i}(1)};
                    idx2z{idxz{i}} = zeros(1,0);
                    idxz{i} = zeros(1,0);
                    pair(counter,2) = i;
                    pair(counter,1) = idxcurated{i}(1);
                    counter = counter +1;
                else
                    idxz{i} = zeros(1,0);
                end
            elseif length(idxz{i})==1 && isempty(idx2z{idxz{i}(1)})
                idxz{i} = zeros(1,0);
            elseif length(idxz{i})>1
                %len = cellfun(@length,idx2z(idxz{i})); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                len = zeros(length(idxz{i}),1);
                for iter = 1:length(idxz{i})
                    len(iter) = length(idx2z{idxz{i}(iter)});
                end
                if sum(len) == 1
                    keepers = len' == 1; 
                    vals = idxz{i}(keepers);
                    idxcurated{i} = vals;
                    idx2curated{vals} = idx2z{vals(1)};
                    idxz{i} = zeros(1,0);
                    idx2z{vals} = zeros(1,0);
                    pair(counter,2) = i;
                    pair(counter,1) = idxcurated{i}(1);
                    counter = counter +1;
                end
            end
        end
    end
    
    %%%% added to remove multiple matching
    for i = 1:length(idxz)
        if length(idxz{i})>1
            a = repmat({zeros(1,0)},length(idxz{i}),1);
            for cis = 1:length(idxz{i})
                a{cis} = idx2z{idxz{i}(cis)};
            end
            for j = length(a):-1:1
                if length(a{j}) == 1 && sum(a{j} ~= i)>0
                    idxz{i}(j) = [];
                elseif isempty(a{j})
                    idxz{i}(j) = [];
                end
            end
        end
    end
    
    for i = 1:length(idx2z)
        if length(idx2z{i})>1
            a = repmat({zeros(1,0)},length(idx2z{i}),1);
            for cis = 1:length(idx2z{i})
                a{cis} = idxz{idx2z{i}(cis)};
            end
            for j = length(a):-1:1
                if length(a{j}) == 1 && sum(a{j} ~= i)>0
                    idx2z{i}(j) = [];
                elseif isempty(a{j})
                    idx2z{i}(j) = [];
                end
            end
        end
    end
    %%%%%
    
    if length(idxz)>length(idx2z)
        for i = 1:length(idxz)
            if ~isempty(idxz{i}) && length(idx2z{idxz{i}(1)}) == 1 && length(idxz{i}) ==1 && sum(idx2z{idxz{i}(1)} == i)>0
                pair(counter,2) = i;
                pair(counter,1) = idxz{i}(1);
                idx2z{idxz{i}} = zeros(1,0);
                idxz{i} = zeros(1,0);
                counter = counter + 1;
            end
        end
    else
        for i = 1:length(idx2z)
            if ~isempty(idx2z{i}) && length(idxz{idx2z{i}(1)}) == 1 && length(idx2z{i}) ==1 && sum(idxz{idx2z{i}(1)} == i)>0
                pair(counter,2) = idx2z{i}(1);
                pair(counter,1) = i;
                idxz{idx2z{i}} = zeros(1,0);
                idx2z{i} = zeros(1,0);
                counter = counter + 1;
            end
        end
    end
    
    for i = 1:length(idxz)
       if length(idxz{i}) > 1
           bart = points{hybnum}(channels(1)).intensity(idxz{i}) == max(points{hybnum}(channels(1)).intensity(idxz{i}));
           if sum(bart) > 1
               bart = false(length(bart),1);
               bart(1) = 1;
           end
           pair(counter,2) = i;
           A = idxz{i}(bart');
           pair(counter,1) = A(1);
           idxz{i} = zeros(1,0);
           counter = counter + 1;
       end
    end

    for i = 1:length(idx2z)
       if length(idx2z{i})==1 && isempty(idxz{idx2z{i}(1)})
           idx2z{i} = zeros(1,0);
       elseif length(idx2z{i}) > 1
           bart = points{hybnum2}(channels(2)).intensity(idx2z{i}) == max(points{hybnum2}(channels(2)).intensity(idx2z{i}));
           if sum(bart) > 1
               bart = false(length(bart),1);
               bart(1) = 1;
           end
           A = idx2z{i}(bart');
           pair(counter,2) = A(1);
           pair(counter,1) = i;
           idx2z{i} = zeros(1,0);
           counter = counter + 1;
       end
    end
    
    if sum(pair) ~= 0
        pair(counter:end,:) = [];
        called1 = zeros(size(conc1,1),1);
        called1(pair(:,1)) = 1;
        called1 = logical(called1);
        called2 = zeros(size(conc2,1),1);
        called2(pair(:,2)) = 1;
        called2 = logical(called2);
        pts1 = conc1(pair(:,1),:);
        pts2 = conc2(pair(:,2),:);
    else
        called1 = false(size(conc1,1),1);
        called2 = false(size(conc2,1),1);
        pts1 = zeros(1,0);
        pts2 = zeros(1,0);
    end
else
    called1 = false(size(conc1,1),1);
    called2 = false(size(conc2,1),1);
    pts1 = zeros(1,0);
    pts2 = zeros(1,0);
end

