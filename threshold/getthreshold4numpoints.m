function [d, logImage, regMax] = getthreshold4numpoints(threshold, target, I, typedots, logImage, regMax)
% get the minimum number of points needed based on the threshold and number
% of points
%
% Inputs: threshold, target, image, typedots ['exons', 'introns']
%
% Outputs: min distance between numberOfPoints and desiredNumPoints
% (target)
%
% Author: Nico Pierson
% Date: August 2018

    [numPoints, logImage, regMax] = numdetectdotsv2(I, threshold, logImage, regMax, typedots);

    % Get minimum of this distance
    d =  abs(target - numPoints);

end