function threshold = threshcsv2mat(csvPath, chArray)
% function to extract threshold from a csv file and convert to a matrix
%
%

    T = readtable(csvPath);
    numHybs = max(T.hyb);
    threshold = ones(numHybs, length(chArray)) * 999999;
    
    for c = chArray
        for h = 1:numHybs
            rows = and(T.channel == c, T.hyb == h);
            threshold(h,c) = T.threshold(rows);
        end
    end
end