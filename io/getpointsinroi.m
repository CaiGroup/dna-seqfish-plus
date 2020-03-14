function keepPoints = getpointsinroi(points, roiPath)
% getpointsinroi gets points and uses the roi to return the points
%
% Author: Nico Pierson
% Date: 4/4/2019
% Email: nicogpt@caltech.edu
% Modified:


    rows = size(points, 1); % size of rows
    keepLog = false(size(rows, 1), 1); % logical for extracting whith rows
    vertex = selfsegzip(roiPath);
    
    for i = 1:length(vertex)
        keepLog = or(keepLog, inpolygon(points(:,1),points(:,2),vertex(i).x,vertex(i).y));
    end

    keepPoints = points(keepLog,:);

end
