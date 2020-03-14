% 2nd phase for decoding points

addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentLabel = 'nochabb-nophysical';

%% variables
posArray = [2 3 4]; %1;%[0 1 2 3]; % total is 0:4
numRounds = 5;
numChannels = 16;
sqrtradius = 6;
segment = 'roi';%'whole';
minseeds = 3;
alloweddiff = 2;
chArray = 1:2;
typedots = 'log';%'exons';
superres = 'radial3d';
folderArray = 0:79;
physicalTforms = [];
usechabboffsets = false;
usephysicaloffsets = false;



% load chaTorm
load('I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points\pre_formated\refpoints-chromaticaberrations.mat', 'chaTform')


for position = posArray
    offsetsName = ['20191125_pos' num2str(position) '_offsets_initial_roithresholdfiltered_UpdatedThresholdPointRatio-logdots-Pos0vs3_UpdatedThresholdPos0-12-11-201920191211_xyse3_zse3_xyte2_zte2_xyme1_zme1.csv'];
    pointsName = ['points-int-thresh-pos' num2str(position) '-AdjustedThreshold-12-12-2019-MedianIntRatio-logdot.mat'];
    %['hyb-points-0_80-ForBeadAlignment-pos' num2str(position) '-radial3d-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH-UpdatedThresholdPointRatio-logdots-Pos0vs3-NoAdjustment.csv'];
    pointsDir = fullfile(experimentDir, 'points');%, 'pre_formated');
    offsetsDir = fullfile(experimentDir, 'points', 'positions');
    offsetsPath = fullfile(offsetsDir, offsetsName);
    pointsPath = getfile(pointsDir, pointsName, 'strict');

    % get points and apply offsets
    [pointsch, offsets] = alignpointswrapper(chArray, pointsPath, offsetsPath, ...
        chaTform, numRounds, folderArray, physicalTforms, usechabboffsets, usephysicaloffsets);


    %% Decode the points
    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, pointsch{channel}, segment, sqrtradius, typedots, superres, alloweddiff, ...
        channel, minseeds, experimentLabel);
    end
end