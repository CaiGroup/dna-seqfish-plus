function [globaltform, threshold, refPoints] = beadglobaltform(beadImPath, position)
% grabs global tform from the bead images
%
% Assumes the images are named 'MMStack_Pos0.ome.tif' where position number
% varies
%
% Dependencies; 
% Fiji Path
% addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end'); 
% bfmatlab path to open images
% addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\');


    %% initial variables
    processBeads = 1;
    typedots = 'log';
    sliding = 0;
    backradius = 3;
    searchradius = 2.45;
    
    %% get the images
    beadImPath = fullfile(beadImPath, ['MMStack_Pos' num2str(position) '.ome.tif']);
    [I, sizeC, sizeZ, ~] = grabimseries(beadImPath, position);
    
    %% Calculate the autothreshold
    threshold = ones(1, sizeC) * 99999;
    points = cell(1, sizeC);
    for ch = 1:sizeC
        logFish = max(I{ch},[],3);
        thresh = multithresh(logFish,2);
        threshold(ch) = thresh(2);
        %[points{ch}, ~, ~, ~] = detectdotsv2(I{ch}, threshold(ch), typedots);
    end

    %% All points are in all channels
    [refPoints, fixedThreshold] = matchpointsims(I(1:sizeC), ...
    threshold, typedots, processBeads, backradius, sliding, ...
    [], searchradius, 0, pwd, typedots);

    %% calculate the global tform
    globaltform = loadchabberationtform(1:sizeC, refPoints);

end