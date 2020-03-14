function threshold = ch2wholethreshold(thresh)
% change the organization of the threshold from round by hyb matrix to hyb
% by 1 matrix

    threshold = [];
    [row, col] = size(thresh);
    for r = 1:row
        threshrow = thresh(1,:)';
        threshold = cat(1, threshold, threshrow);
    end
    
end