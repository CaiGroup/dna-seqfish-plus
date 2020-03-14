function [radialspots, intensity] = getradialcenter(rawSpots, rawImage)
% super resolves spots with the radial center algorithm
%
% Dependencies: radialcenter.m
%
% Author: Mike Lawson
% Date: 2018
% Adapted by: Nico Pierson
% Date:7/16/2019
% Email: nicogpt@caltech.edu

    % Variables
    Pixels = 3;
    xc = [];
    yc = [];
    zc = [];
    [ys,xs,zs] = size(rawImage);
    intensity = zeros(length(rawSpots),1);
    for dotID = 1:length(rawSpots)
        zc(dotID) = rawSpots(dotID, 3);
        I = double((rawImage([-Pixels:Pixels]+rawSpots(dotID,2),[-Pixels:Pixels]+rawSpots(dotID,1), zc(dotID))));
        [xc(dotID),yc(dotID)] = radialcenter(I);
    end
    radialspots = rawSpots(:,1:2) + [xc; yc]' - [Pixels, Pixels]-1;
    radialspots(:,3) = zc;
    for idx = 1:size(radialspots,1)
        x = round(radialspots(idx,1));
        if x > xs
            x = xs;
        elseif x < 1
            x = 1;
        end
        y = round(radialspots(idx,2));
        if y > ys
            y = ys;
        elseif y < 1
            y = 1;
        end
        z = radialspots(idx,3);
        intensity(idx) = rawImage(y,x,z);
    end
    

end