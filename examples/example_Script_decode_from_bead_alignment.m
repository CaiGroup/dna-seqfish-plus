%script to decode the points after the transformations

addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');


experimentDir = 'I:\2019-09-09-brain-rep2-2-DNAFISH';
% Set up Variables
position = 3;
numRounds = 5;
numChannels = 16;
sqrtradius = 6;
segment = 'whole';
experimentName = '2019-09-09-brain-rep2-2-DNAFISH'; 
minseeds = 3;
alloweddiff = 2;
chArray = 1:2;
typedots = 'log';
superres = 'radial3d';


% load chromatic aberrations
chaDir = 'I:\2019-09-09-brain-rep2-2-DNAFISH\points';
chaPath = getfile(chaDir, ['beadFinalRefPoints_chaTforms-pos' num2str(position) '.mat'], 'match');
load(chaPath, 'chaTform');


offsetsBaseDir = ['C:\Users\Long Cai - 2\Box\Jonathan\2019-09-09-brain-rep2-2-DNAFISH\phase2_many_hybs\positions\pos' num2str(position)];  
pointsDir = 'C:\Users\Long Cai - 2\Box\Jonathan\2019-09-09-brain-rep2-2-DNAFISH\phase2_many_hybs\pre_formated_aligned2hyb0';
[pointsch, offsets] = alignpointswrapper(chArray, position, pointsDir, offsetsBaseDir, chaTform);

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