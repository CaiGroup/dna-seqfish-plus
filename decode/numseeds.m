function [seeds] = numseeds(dotlocations)
%#codegen

seeds = repmat({0},size(dotlocations,1),1);
% for each cell i
    for j = 1:size(dotlocations,1)
        % find the conversion from dotlocation id to global geneID (ind)
        %ind = find(strcmpi(dotlocations(i).cell{j,1},PosList));
        % check if there's only one such gene called (i+1 because PosList
        % has gene names in the first column
        if size(dotlocations{j,4},1) == 1
            % if this is the only time this gene is called, then all
            % seedings are associated with this one call and the number of
            % seeds for this RNA call = the number of seedings called for
            % this gene
            seeds{j} = size(dotlocations{j,3},1);
        elseif size(dotlocations{j,4},1) == 0
            seeds{j} = [];
        else
            % if a gene is called multiple times (has multiple centroids),
            % then for each FISH dot [dotlocations(i).cell{j,2}] it finds
            % the closest RNA centroid [PosList{ind,i+1}]
            idx = knnsearch(dotlocations{j,4},dotlocations{j,1});
            % uses hist to bin the centroid ids associated with each FISH
            % dot, so the result is the number of FISH dots associated with
            % each centroid
            seeds{j} = histc(idx(:,1),unique(idx(:,1)));
        end
    end
