% function to adjust the threshold for all positions by using a multiplier
% to compare the median intensity of the captured points - use E14 rep3-1
% as the reference

% get the images and points using the reference threshold and calculate the
% median reference intensity across all hybridizations....

% create a table and compare to the reference and adjust the multiplier to
% find a closer median intensity across all positions

addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
experimentDir = 'I:\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH';
experimentName = '2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH';


%% reference median intensity for capture points in each hybridization
pointsRefPath = 'I:\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\points\hyb-points-0_80-ForBeadAlignment-pos0-radialcenter3d-2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH.csv';
chArray = 1:3;
[refPoints, refIntensity] = csv2hybxcell(pointsRefPath, chArray);
% need function to get the points and the intensity for the firs position
% that is manually thresholded


% get the median ref intensity
medianIntensity = refIntensity;
for i = 1:size(medianIntensity,1)
    for ch = 1:3
        medianIntensity{i,ch} = median(medianIntensity{i,ch});
    end
end

%% load pos1-4 images, grab points using same treshold vales
% load threshold
load('I:\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\threshold\thresholdAllCh-pos0-numHybCycles85-numCh-3-E14-DNA-seqFISH+rep3-1-DNAFISH.mat');
% get images and grab points
position = 1;
folderArray = 0:79;
numCh = 3;
I = cell(length(folderArray), numCh);
dapiI = cell(length(folderArray), 1);
allIms = cell(length(folderArray),1);
folderArrayIdx = folderArray + 1;
tic
parfor f = folderArrayIdx
    folder = f + 1;
    fprintf('Retrieving Position %.0f Folder %.0f images\n', position, f-1);
    imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
    imagePath = fullfile(experimentDir, ['HybCycle_' num2str(f-1)], imageName);
    [allIms{f}, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
end
toc
I = allIms(:,1:sizeC-1);
dapiI = allIms(:,sizeC);


%% Align Dapi for Background Images for All Positions
fprintf('Aligning Dapi for Background Images...\n');
backgroundFolderName = 'initial_background';
backImBasePath = fullfile(experimentDir, backgroundFolderName);
backImPath = fullfile(backImBasePath, ['MMStack_Pos' num2str(position) '.ome.tif']);
if exist(backImPath, 'file') == 2
    [backIms, numDapi, numZSlice, ~, ~] = grabimseries(backImPath, position);
    numCh = numDapi - 1;

    %% Subtract the background and Multiply by Shading Corrections
    % Get Shading Corrections
    subtractBackground = false;
    imageJBackSubtract = true;
    shadingcorr = shadingcorrection(backIms(1:numCh));
    parfor f = 1:length(folderArray)
        for ch = 1:numCh

            % Apply the shading correctionsmean
            I{f,ch} = uint16(double(allIms{f}{ch}) ./ double(shadingcorr{ch}));

            if imageJBackSubtract
                % ImageJ Rolling Ball Back Subtract to remove noise using rad 3
                % replace with deconvolution - need to test first
                uniqueString = 'imageTempProcess-90jf03j';
                I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString,...
                    experimentDir);
            end

        end
    end
else
    error 'background directory or images were not found';
end



%% Get the Ref beads and chromatic aberrations from the initial and final beads
beadFolderName = 'initial_fiducial_markers';
[chaTform, refPointsAligned] = getrefbeads(experimentDir, position, beadFolderName);
%% QC check for chromatic aberration corrections
% if the initial fiducial marker chaTforms will line up the final fiducial markers



%% get the points
points = cell(length(folderArray),numCh);
intensity = points;
medianError = points;
xyPixSize = 1;
zPixSize = 1;
adjustedThreshold = ones(length(folderArray),numCh) * 999999;
for ch = 1:numCh
    parfor f = 1:length(folderArray)
        [points{f,ch}, intensity{f,ch}, adjustedThreshold(f,ch), medianError{f,ch}] = adjustthreshold(I{f,ch}, medianIntensity{f,ch}, threshold(f,ch), typedots);
    end
end

%% save points to Csv files 
saveDir = fullfile(experimentDir, 'points');
if exist(saveDir, 'dir') ~= 7
mkdir(saveDir);
end


listSavePath = fullfile(saveDir, ['hyb-points-0_80-ForBeadAlignment-pos' num2str(position) '-radialcenter3d-' experimentName '.csv']);
fileID = fopen(listSavePath,'w');
fprintf(fileID, '%s,%s,%s,%s,%s,%s\n', 'ch', 'hyb', 'x', 'y', 'z', 'int');

for hyb = 1:length(folderArray)
    for ch = 1:numCh
        pointsSize = length(points{hyb,ch});
        %% apply chromatic aberration corrections to the points
        points{hyb,ch} = transformPointsForward(chaTform{ch}, points{hyb,ch});
        for i = 1:pointsSize
            x = points{hyb,ch}(i,1);
            y = points{hyb,ch}(i,2);
            z = points{hyb,ch}(i,3);
            int = intensity{hyb,ch}(i);
            fprintf(fileID, '%.0f,%.0f,%.3f,%.3f,%.3f,%.0f\n', ch, hyb, x, y, z, int);
        end
    end
end

savePath = fullfile(saveDir, ['hyb-points-0_80-ForBeadAlignment-pos' num2str(position) '-radialcenter3d-' experimentName '.mat']);
save(savePath, 'points', 'intensity', 'shadingcorr');
        
        
%end


%% use bead alignment code - separate the code here in the bash script or batch file



%% Decode from points
% Set up Variables
numRounds = 5;
numChannels = 16;
sqrtradius = 6;
segment = 'whole';
experimentName = '2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH'; 
minseeds = 3;
alloweddiff = 2;
chArray = 1:2;
typedots = 'log';
superres = 'radial3d';

offsetsBaseDir = [];  
pointsDir = 'I:\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\points';
[pointsch, offsets] = alignpointswrapper(chArray, position, pointsDir, offsetsBaseDir, chaTform); % need to adjust this to fix the pipeline

% save data
save([experimentDir '\points\pointsAlignedBeadsJonathan-pos' num2str(position) '-chs1_2.mat'], 'pointsch', 'offsets', 'experimentDir', 'position', 'numRounds', ...
    'segment', 'sqrtradius', 'minseeds', 'alloweddiff', 'experimentName', 'chArray', 'chaTform');




%% Decode the points


% To load the processed images
%dataDir = fullfile(experimentDir, 'processedimages\pos0\preProcessedData-pos0-2019-09-09-brain-rep2-2-DNAFISH-2019-09-28.mat');
%load(dataDir, 'I');


for channel = chArray
    [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
    ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
    position, numRounds, numChannels, pointsch{channel}, segment, sqrtradius, typedots, superres, alloweddiff, ...
    channel, minseeds);
end

