% Example Script to run the seqFISH pipleline
% Tested on E14-rep2-2
%
% To Do:
% 1. threshold function is specific to grab column 'Pos0'
% 2. add option to use 3d segmentation
% 3. update parameters for bead alignment
% 4. debug brain dataset
%
% Date: 12/19/2019


%% Variables
experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentLabel = '12-19-2019-pipeline-test';
superres = 'radial3d'; % 'radial3d', 'gaussian'
thresholdadjust = true; %false; % option to adjust threshold
posArray = 0;%[0 1 2 3]; %0:4;
posArrayChaTform = 0:4;
folderArray = 0:79;
numCh = 2;
sizeZ = 25;
numRounds = 5;
numChannels = 16;
sqrtradius = 6;
minseeds = 3;
alloweddiff = 2;
usechabboffsets = true;
usephysicaloffsets = true;
filtersigma = false;
longestEdges = 20; 
minEdgeMatch = 3; 
minDotMatch = 4;   
numRefPoints = [];
medRefIntensity = [];



%% save global variables
saveGlobalDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'variables');
if exist(saveGlobalDir, 'dir') ~= 7
    mkdir(saveGlobalDir);
end
dateStart = datetime;
formatDate = 'yyyy-mm-dd';
dateString = datestr(dateStart, formatDate);
saveGlobalName = ['global-variables-seqfish-pipeline_' experimentLabel '_' dateString];
saveGlobalPath = fullfile(saveGlobalDir, saveGlobalName);
save(saveGlobalPath, 'experimentDir', 'experimentName', 'experimentLabel', 'superres', ...
    'thresholdadjust', 'posArray', 'folderArray', 'numCh', 'sizeZ', 'numRounds', ...
    'numChannels', 'sqrtradius', 'alloweddiff', 'usechabboffsets', 'usephysicaloffsets', ...
    'filtersigma', 'numRefPoints', 'medRefIntensity', 'longestEdges', 'minEdgeMatch', ...
    'minDotMatch');
       


%% Chromatic aberration and physical offsets
[chaTform, ~, physicalTform, refPointsAlignedInitial, refInitialSaveName] ...
    = debugbeads(experimentDir, posArrayChaTform, numCh, sizeZ, superres, experimentLabel);


for position = posArray

    %% process - ref position assumed to be 0
    [rawpoints, intensity, pointsCsvName, numRefPoints, medRefIntensity] ...
        = seqfishprocess(experimentDir, experimentName, experimentLabel, ...
        position, folderArray, thresholdadjust, superres, numChannels, ...
        numRounds, sizeZ, chaTform, numRefPoints, medRefIntensity, ...
        filtersigma);



    %% decode
    seqfishdecode(experimentDir, experimentName, experimentLabel, position, ...
        numRounds, numChannels, sqrtradius, minseeds, alloweddiff, numCh, ...
        folderArray, chaTform, physicalTform{position+1}, usechabboffsets, usephysicaloffsets, ...
        refInitialSaveName{position+1}, pointsCsvName, longestEdges, minEdgeMatch, minDotMatch); 
end