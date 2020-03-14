function [points, intensity] = csv2hybxcell(pointsPath, chArray)
% converts the csv to a hyb by cell points structure

    T = readtable(pointsPath);
    numHybs = max(T.hyb);
    points = cell(numHybs, length(chArray));
    intensity = cell(numHybs, length(chArray));
    
    for c = chArray
        for h = 1:numHybs
            rows = and(T.ch == c, T.hyb == h);
            points{h,c} = [T.x(rows), T.y(rows), T.z(rows)];
            intensity{h,c} = T.int(rows);
        end
    end

end