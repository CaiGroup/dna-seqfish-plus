
% edit by Yodai Takei
% Last updated on 01/16/20.

%% Experiment dependent Variables
addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts');
experimentDir = 'E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH';
experimentLabel = '2020-02-11-sequential-control-hyb81-84-hyb80rehyb';
posArray = 0:4; % for decoding. needs Pos0 anyway. Don't start from Pos1 or later.
posArrayAll = 0:4;% for fiducial markers.
sizeZ = 25; % number of z slices. in the future, this should be automatically calculted.
folderNames = ["HybCycle_0","HybCycle_initial_hyb80readouts","HybCycle_79","HybCycle_80","HybCycle_81","HybCycle_82",...
               "HybCycle_83","initial_fiducial_markers","final_fiducial_markers"];
thresh_folder = '2020-01-16-sequential-control-hyb81-84-test\threshold'; % threshold folder name in the experimentDir

%% Other fixed variables for DNA seqFISH+. Some may still need to be changed to compare to.
numChAll = 3;
chArrayAll = 1:numChAll;
numCh = 3;
chArray = 1:numCh;
folderArray = 0:length(folderNames)-1;
numRounds = 1; % barcoding rounds
numChannels = length(folderArray); % barcoding pseudo-channels.
typedots = 'log';
superres = 'radial3d'; % 'radial3d', 'gaussian'
usechabboffsets = true; % should be always true.
usephysicaloffsets = true; % default is true.
filtersigma = true;
longestEdges = 20; 
minEdgeMatch = 3; 
minDotMatch = 4;   
segment = 'roi';
saveprocessedImages = true;

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
save(saveGlobalPath, 'experimentDir', 'experimentName', 'experimentLabel', 'typedots', 'superres', ...
    'posArray', 'folderNames', 'folderArray', 'numChAll','sizeZ', 'numRounds', ...
    'numChannels', 'usechabboffsets', 'usephysicaloffsets', ...
    'filtersigma', 'longestEdges', 'minEdgeMatch', ...
    'minDotMatch');

load('E:\Yodai\DNAFISH+\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\analysis\2020-01-24-pipeline-full-decoding\points\variables\refpoints-chromaticaberrations-initial-final-2020-01-24.mat',...
    'chaTformInitial','physicalTform');
chaTform = chaTformInitial;

%% Step1-6: for DNA seqFISH+ analysis.
%% Step1: Chromatic aberration and physical offsets
%[chaTform, ~, physicalTform, refPointsAlignedInitial, refInitialSaveName] ...
%    = debugbeads(experimentDir, posArrayAll, numChAll, sizeZ, superres, experimentLabel);
% comment this only when you run second time, when it stopped inside the pos for loop.

for position = posArray
    
    %% Step2: Process - grab the raw points from preprocessed images.
    [points, intensity, pointsCsvName, I] ...
        = seqfishprocess_initialpoints_specificfolder(experimentDir, experimentName, experimentLabel, ...
        position, folderNames, typedots, superres, chaTform, filtersigma, numChAll, thresh_folder);
    
    %if position == 0
    %    saveDirPos0= fullfile(experimentDir, 'preprocessedimages_unaligned',  experimentLabel);
    %    if exist(saveDirPos0, 'dir') ~= 7
    %        mkdir(saveDirPos0);
    %    end
    %    savePathPos0 = fullfile(saveDirPos0, ['preprocessedimages_unaligned-pos' num2str(position) '-' experimentLabel '.mat']);
    %    save(savePathPos0, 'I', '-v7.3');
    %end
    
    %% Step3: Compute tforms for all hybs with fiducial alignment and output aligned/corrected points and offsets
    [pointsch, offsets, corrpoints] = seqfishtformalign(experimentDir, experimentLabel, ...
        position, numRounds, numChAll, folderArray, chaTform, physicalTform, ...
        usechabboffsets, usephysicaloffsets, refInitialSaveName{position+1}, pointsCsvName, longestEdges, minEdgeMatch, minDotMatch);
    % pointsch: already barcoded cell format
    % corrpoints: sequential cell format, cell(numFolders, numCh)
    % offsets: reference image/points are hyb1 ch1.
    
    % save chromatic aberration corrected preprocessed/aligned images.
    I = seqfish_outputpreprocessedimagesAll(experimentDir, experimentName, experimentLabel, I, offsets, chaTform, physicalTform{position+1}, position, usechabboffsets, usephysicaloffsets, saveprocessedImages);
    clearvars I
    
    pointsfinal = cell(numChannels, 1);
    for f = 1:numChannels
        for ch = 1:numCh
            pointsfinal{f}(ch).channels = pointsch{ch,1}{1,1}(f).channels;
            pointsfinal{f}(ch).intensity = pointsch{ch,1}{1,1}(f).intensity;
        end
    end
    
    [pointschpercell, numCells] = segment2dcells_outputpointspercell(experimentDir, experimentLabel, chArrayAll, position, pointsfinal);
    saveDir= fullfile(experimentDir, 'analysis',  experimentLabel, 'variables');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    savePath = fullfile(saveDir, ['outputVariables-pos' num2str(position) '-' experimentLabel '.mat']);
    save(savePath, 'points', 'pointsch', 'pointsfinal','intensity', 'pointschpercell', 'numCells', 'chArray', 'folderNames', 'position');
    
end