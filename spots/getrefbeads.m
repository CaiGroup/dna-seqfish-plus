function [chaTform, refPointsAligned, intensity] = getrefbeads(experimentDir, position, beadfoldername) 
% gets bead ref points from bead images and transform them to the hyb0 ref
% dapi image
%
% Path Dependencies
%addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
%addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
%addpath('C:\github\streamline-seqFISH\src\AlignImages\', '-end');
%addpath('C:\github\streamline-seqFISH\src\FindThreshold\', '-end');

    saveDir = fullfile(experimentDir, 'points');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end

    backradius = 3;
    sliding = false;
    searchradius = 3;
    processBeads = true;
    typedotsFixed = 'log';
    


    % initial beads
    %beadfoldername = 'initial_fiducial_markers';
    imBeadsPath = fullfile(experimentDir, beadfoldername, ['MMStack_Pos' num2str(position) '.ome.tif']);
    [beadImsTemp, numBeadDapi, numBeadZSlice, ~, ~] = grabimseries(imBeadsPath, position);
    numCh = numBeadDapi -1;
    fixedThreshold = ones(numCh, 1) * 99999;
    for ch = 1:numCh
        logFish = max(beadImsTemp{ch},[],3);
        thresh = multithresh(logFish,2);
        if thresh(2) < 550 || thresh(2) > 15000
            fixedThreshold(ch) = fixedThreshold(ch-1);
        else
            fixedThreshold(ch) = thresh(2);
        end
    end
    
    [refPoints, fixedThreshold] = matchpointsims(beadImsTemp(1:numCh), ...
    fixedThreshold, typedotsFixed, processBeads, backradius, sliding, ...
    [], searchradius, 0, pwd, 'log');
    maxZ = size(beadImsTemp{1},3);
    maxXY = 2048;
    min = 1;

    % Apply the tform of dapi from final beads to the dapi of Hyb0
    % get Dapi from Hyb0
    dapiPath = fullfile(experimentDir, 'HybCycle_0', ['MMStack_Pos' num2str(position) '.ome.tif']);
    [allhyb0, numhyb0Dapi, numhyb0ZSlice, ~, ~] = grabimseries(dapiPath, position);
    dapiRef = allhyb0{numhyb0Dapi};
    % find tform from the bead dapi to the dapiRef
    beadDapi = beadImsTemp{numBeadDapi};
    bead2dapiTform = grabtform(beadDapi, dapiRef);



    intensity = cell(1,numCh);
    for ch = 1:numCh
            % apply transformation to the points
            tempPoints = refPoints(ch).channels;
            for i = 1:length(tempPoints)
                y = round(tempPoints(i,2));
                x = round(tempPoints(i,1));
                z = round(tempPoints(i,3));
                if z > maxZ
                    z = maxZ;
                elseif z < min
                    z = min;
                elseif isnan(z)
                    z = round(tempPoints(i-1,3));
                end
                if x > maxXY
                    x = maxXY;
                elseif x < min
                    x = min;
                end
                if y > maxXY
                    y = maxXY;
                elseif y < min
                    y = min;
                end
                intensity{ch} = cat(1, intensity{ch}, beadImsTemp{ch}(y, x, z));
            end
    end

    refPointsAligned = refPoints;
    
    for ch = 1:numCh
        % apply transformation to the points
        refPointsAligned(ch).channels = transformPointsForward(bead2dapiTform, refPoints(ch).channels);
    end

    
    %% get chromatic abberations
    chaTform = loadchabberationtform(1:numCh, refPointsAligned);
end