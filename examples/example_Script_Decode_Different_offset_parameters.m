% 2nd phase for decoding points

addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';







%% variables
posArray = 4;%[1 2 3]; % total is 0:4
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


% get physical offsets
beadDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points\pre_formated';
finalFile = 'ref-points-ForBeadAlignment-pos%.0f-final_fiducial_markers-FilteredByROI.csv';
initialFile = 'ref-points-ForBeadAlignment-pos%.0f-initial_fiducial_markers-FilteredByROI.csv';
numCh = 2;


% load chaTorm
load('I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points\pre_formated\refpoints-chromaticaberrations.mat', 'chaTform')


for position = posArray
    %% Align the Beads
    home_dir = '"I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points"';
    ref_fname = '"ref-points-ForBeadAlignment-pos%d-radial3d-initial_fiducial_markers-raw-intensity.csv"';
    %ref_fname = '"ref-points-ForBeadAlignment-pos%d-initial_fiducial_markers-FilteredByROI-manualThreshold-log3d.csv"';
    ro_fname = '"hyb-points-0_80-ForBeadAlignment-pos%d-radialcenter3d-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH-AdjustedThreshold-12-12-2019-MedianIntRatio-logdot.csv"';
    savefname = '"20191125_pos%d_offsets_initial_AdjustedPointRatio-logdots-"';
    savefnameMat = '20191125_pos%d_offsets_initial_AdjustedPointRatio-logdots-';
    % Run the python comand for bead alignment
    addpath('C:\github\streamline-seqFISH\src\beadalignment\', '-end');
    pythonCommand = ['python C:\github\streamline-seqFISH\src\beadalignment\20191127_main_par_matlab.py ' home_dir ' ' ref_fname ' ' ro_fname ' ' savefname ' ' num2str(position)];
    system(pythonCommand);
    
    
    
    %% no offsets
    experimentLabel = 'nochabb-nophysical';
    usechabboffsets = false;
    usephysicaloffsets = false;
    
    offsetsName = sprintf(savefnameMat, position);%['20191125_pos' num2str(position) '_offsets_initial_AdjustedPointRatio-logdots-20191213_xyse3_zse3_xyte2_zte2_xyme1_zme1.csv'];
    
    pointsName = ['points-int-thresh-pos' num2str(position) '-AdjustedThreshold-12-12-2019-MedianIntRatio-logdot.mat'];
    %['hyb-points-0_80-ForBeadAlignment-pos' num2str(position) '-radial3d-2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH-UpdatedThresholdPointRatio-logdots-Pos0vs3-NoAdjustment.csv'];
    pointsDir = fullfile(experimentDir, 'points');%, 'pre_formated');
    offsetsDir = fullfile(experimentDir, 'points', 'positions');
    
    offsetsPath = getfile(offsetsDir, offsetsName, 'match');
    pointsPath = getfile(pointsDir, pointsName, 'strict');

    % get points and apply offsets
    [pointsch, offsets] = alignpointswrapper(chArray, pointsPath, offsetsPath, ...
        chaTform, numRounds, folderArray, physicalTforms{position+1}, usechabboffsets, usephysicaloffsets);


    %% Decode the points
    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, pointsch{channel}, segment, sqrtradius, alloweddiff, ...
        channel, minseeds, experimentLabel);
    end
    %[sqrtradius, alloweddiff, channel, minseeds, experimentLabel, removePoints, dotlocations, removeInd]
    
    
    %% test using chromatic aberration
    experimentLabel = 'chabb-nophysical';
    usechabboffsets = true;
    usephysicaloffsets = false;
    
       % get points and apply offsets
    [pointsch, offsets] = alignpointswrapper(chArray, pointsPath, offsetsPath, ...
        chaTform, numRounds, folderArray, physicalTforms{position+1}, usechabboffsets, usephysicaloffsets);

    
    %% Decode the points
    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, pointsch{channel}, segment, sqrtradius, alloweddiff, ...
        channel, minseeds, experimentLabel);
    end
    
    
    
    %% test using chromatic aberration
    experimentLabel = 'chabb-physical';
    usechabboffsets = true;
    usephysicaloffsets = true;
    
       % get points and apply offsets
    [pointsch, offsets] = alignpointswrapper(chArray, pointsPath, offsetsPath, ...
        chaTform, numRounds, folderArray, physicalTforms{position+1}, usechabboffsets, usephysicaloffsets);

    
    %% Decode the points
    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, pointsch{channel}, segment, sqrtradius, alloweddiff, ...
        channel, minseeds, experimentLabel);
    end
    
end