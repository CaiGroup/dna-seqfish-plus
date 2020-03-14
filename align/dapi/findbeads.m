function [] = findbeads(experimentDir, experimentName, position)
% function prints the initial beads and final beads

% addpaths
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');

saveDir = fullfile(experimentDir, 'points');
if exist(saveDir, 'dir') ~= 7
    mkdir(saveDir);
end

% variables
backradius = 3;
sliding = false;
searchradius = 3;
processBeads = true;
typedotsFixed = 'log';
addpath('C:\github\streamline-seqFISH\src\AlignImages\', '-end');
addpath('C:\github\streamline-seqFISH\src\FindThreshold\', '-end');
addpath('C:\github\streamline-seqFISH\src\AlignImages\bfmatlab\', '-end');


% initial beads
beadfoldername = 'initial_fiducial_markers';
imBeadsPath = fullfile(experimentDir, beadfoldername, ['MMStack_Pos' num2str(position) '.ome.tif']);
[beadImsTemp, numBeadDapi, numBeadZSlice, ~, ~] = grabimseries(imBeadsPath, position);
numHyb = numBeadDapi -1;
fixedThreshold = []; %fixedThreshold = [10000 20000 15000]; %fixedThreshold = [5000 12000 18000]; % ch1 %fixedThreshold = [10000 20000 15000]; % Can add threshold if already set
[refPoints, fixedThreshold] = matchpointsims(beadImsTemp(1:numHyb), ...
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



intensity = cell(1,3);
for ch = 1:3
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
listSavePath = fullfile(saveDir, ['ref-points-ForBeadAlignment-pos' num2str(position) '-radial3d-initialBeads-raw-intensity-' experimentName '.csv']);
fileID = fopen(listSavePath,'w');
fprintf(fileID, '%s,%s,%s,%s,%s\n', 'ch', 'x', 'y', 'z', 'int');
for ch = 1:3
    pointsSize = length(refPoints(ch).channels);
    % apply transformation to the points
    refPointsAligned(ch).channels = transformPointsForward(bead2dapiTform, refPoints(ch).channels);
    pointsTemp = refPointsAligned(ch).channels;
    for i = 1:pointsSize
        x = pointsTemp(i,1);
        y = pointsTemp(i,2);
        z = pointsTemp(i,3);
        int = intensity{ch}(i);
        fprintf(fileID, '%.0f,%.3f,%.3f,%.3f,%.0f\n', ch, x, y, z, int);
    end
end
fclose(fileID);



% final beads
beadfoldername = 'final_fiducial_markers';
imBeadsPath = fullfile(experimentDir, beadfoldername, ['MMStack_Pos' num2str(position) '.ome.tif']);
[beadImsTemp, numBeadDapi, numBeadZSlice, ~, ~] = grabimseries(imBeadsPath, position);
numHyb = numBeadDapi -1;
fixedThreshold = []; %fixedThreshold = [7000 12000 10000]; %fixedThreshold = []; %
[refPoints, fixedThreshold] = matchpointsims(beadImsTemp(1:numHyb), ...
fixedThreshold, typedotsFixed, processBeads, backradius, sliding, ...
[], searchradius, 0, pwd, 'log');

% Apply the tform of dapi from final beads to the dapi of Hyb0
% get Dapi from Hyb0
dapiPath = fullfile(experimentDir, 'HybCycle_0', ['MMStack_Pos' num2str(position) '.ome.tif']);
[allhyb0, numhyb0Dapi, numhyb0ZSlice, ~, ~] = grabimseries(dapiPath, position);
dapiRef = allhyb0{numhyb0Dapi};
% find tform from the bead dapi to the dapiRef
beadDapi = beadImsTemp{numBeadDapi};
bead2dapiTform = grabtform(beadDapi, dapiRef); 
maxZ = size(beadImsTemp{1},3);
maxXY = 2048;
min = 1;


intensity = cell(1,3);
for ch = 1:3
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
listSavePath = fullfile(saveDir, ['ref-points-ForBeadAlignment-pos' num2str(position) '-radial3d-finalBeads-raw-intensity-' experimentName '.csv']);
fileID = fopen(listSavePath,'w');
fprintf(fileID, '%s,%s,%s,%s,%s\n', 'ch', 'x', 'y', 'z', 'int');
for ch = 1:3
    pointsSize = length(refPoints(ch).channels);
    % apply transformation to the points
    refPointsAligned(ch).channels = transformPointsForward(bead2dapiTform, refPoints(ch).channels);
    %refPointsAligned(ch).channels = refPoints(ch).channels;
    
    pointsTemp = refPointsAligned(ch).channels;
    for i = 1:pointsSize
        x = pointsTemp(i,1);
        y = pointsTemp(i,2);
        z = pointsTemp(i,3);
        int = intensity{ch}(i);
        fprintf(fileID, '%.0f,%.3f,%.3f,%.3f,%.0f\n', ch, x, y, z, int);
    end
end
fclose(fileID);



% get the chatforms from final beads
chArray = 1:3;
chaTform = loadchabberationtform(chArray, refPoints);
save(fullfile(saveDir, ['beadInitialRefPoints_chaTforms-pos' num2str(position) '-' experimentName '.mat']), 'refPointsAligned', 'chaTform', 'refPoints', 'fixedThreshold', 'bead2dapiTform');

end