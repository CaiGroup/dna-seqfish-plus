function [] = seqfishdecode(experimentDir, experimentName, experimentLabel, ...
    position, numRounds, numChannels, sqrtradius, minseeds, alloweddiff, numCh, ...
    folderArray, chaTform, physicalTform, usechabboffsets, usephysicaloffsets, ...
    refSaveName, hybSaveName, longestEdges, minEdgeMatch, minDotMatch)   
    
    % 2nd phase for decoding points
    addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    addpath('C:\github\streamline-seqFISH\src\beadalignment\', '-end');


    
    %% variables
    segment = 'roi';%'whole';
    chArray = 1:numCh;
    pointsDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points');%, 'pre_formated');
    preformatDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'pre_formated');
    if exist(preformatDir, 'dir') ~= 7
        mkdir(preformatDir);
    end
    offsetsDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'positions');
    if exist(offsetsDir, 'dir') ~= 7
        mkdir(offsetsDir);
    end



    %% Align using Beads
    home_dir = strcat('"', pointsDir, '"');
    ref_fname = strcat('"', refSaveName, '"');
    ro_fname = strcat('"', hybSaveName , '"');
    savefname = '"point-offsets_pos%d_"';
    savefnameMat = 'point-offsets_pos%d_';
    % Run the python comand for bead alignment
    pythonCommand = ['python "C:\Users\Long Cai - 1\Desktop\code\streamline-seqFISH-master-20191219-2\streamline-seqFISH-master\src\beadalignment\20191127_main_par_matlab.py" ' ...
        home_dir ' ' ref_fname ' ' ro_fname ' ' savefname ' ' num2str(position) ' ' ...
        num2str(longestEdges) ' ' num2str(minEdgeMatch) ' ' num2str(minDotMatch)];
    system(pythonCommand);

    

    %% get offsets from mat file     
    offsetsName = sprintf(savefnameMat, position);%['20191125_pos' num2str(position) '_offsets_initial_AdjustedPointRatio-logdots-20191213_xyse3_zse3_xyte2_zte2_xyme1_zme1.csv'];
    pointsName = sprintf('points-int-thresh-pos%d', position);
    offsetsPath = getfile(offsetsDir, offsetsName, 'match');
    pointsPath = getfile(pointsDir, pointsName, 'match');

    
    
    % get points and apply offsets
    [pointsch, offsets] = alignpointswrapper(chArray, pointsPath, offsetsPath, ...
        chaTform, numRounds, folderArray, physicalTform, usechabboffsets, usephysicaloffsets);

    

    %% Decode the points
    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, pointsch{channel}, segment, sqrtradius, alloweddiff, ...
        channel, minseeds, experimentLabel);
    end


end