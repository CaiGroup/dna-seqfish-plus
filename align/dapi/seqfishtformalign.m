function [pointsch, offsets, corrpoints] = seqfishtformalign(experimentDir, experimentLabel, ...
    position, numRounds, numCh, ...
    folderArray, chaTform, physicalTform, usechabboffsets, usephysicaloffsets, ...
    refSaveName, hybSaveName, longestEdges, minEdgeMatch, minDotMatch)   
    
    % 2nd phase for decoding points
    %addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
    %addpath('C:\Users\Long Cai - 1\Desktop\Fiji.app\scripts\', '-end');
    %addpath('C:\github\streamline-seqFISH\src\beadalignment\', '-end');


    
    %% variables
    segment = 'roi';%'whole';
    chArray = 1:numCh;
    numHybs = length(folderArray);
    pointsDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points');%, 'pre_formated');
    preformatDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'pre_formated');
    if exist(preformatDir, 'dir') ~= 7
        mkdir(preformatDir);
    end
    offsetsDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'positions');
    if exist(offsetsDir, 'dir') ~= 7
        mkdir(offsetsDir);
    end
    
    pos = 0; % always set to zero
    hybSaveNameNew = 'hybridization-points-pos%d.csv';
    hybCsvOldPath = fullfile(preformatDir, hybSaveName);
    hybCsvNewPath = fullfile(preformatDir, sprintf(hybSaveNameNew, pos));
    copyfile(hybCsvOldPath, hybCsvNewPath);
    refSaveNameNew = 'reference-points-pos%d.csv';
    refCsvOldPath = fullfile(preformatDir, refSaveName);
    refCsvNewPath = fullfile(preformatDir, sprintf(refSaveNameNew, pos));
    copyfile(refCsvOldPath, refCsvNewPath);



    %% Align using Beads
    home_dir = strcat('"', pointsDir, '"');
    %ref_fname = strcat('"', refSaveName, '"');
    %ro_fname = strcat('"', hybSaveName , '"');
    ref_fname = strcat('"', refSaveNameNew, '"');
    ro_fname = strcat('"', hybSaveNameNew , '"');
    savefnameMat = 'point-offsets_pos%d_';
    savefname = strcat('"', sprintf(savefnameMat, position), '"');
    %savefnameMat = 'point-offsets_pos%d_';
    %savefname = '"point-offsets_pos%d_"';
    %savefnameMat = 'point-offsets_pos%d_';
    % Run the python comand for bead alignment
    
    % previous version
    %pythonCommand = ['python "C:\Users\Long Cai - 1\Desktop\code\streamline-seqFISH-master-20191219-2\streamline-seqFISH-master\src\beadalignment\20191127_main_par_matlab.py" ' ...
    %    home_dir ' ' ref_fname ' ' ro_fname ' ' savefname ' ' num2str(pos) ' ' ...
    %    num2str(longestEdges) ' ' num2str(minEdgeMatch) ' ' num2str(minDotMatch)];
    
    % new version for more robust alignment correction.
    pythonCommand = ['python "C:\Users\Long Cai - 1\Desktop\code\streamline-seqFISH-master-20191219-2\streamline-seqFISH-master\src\beadalignment\20191218_align_fill_ch_matlab.py" ' ...
        home_dir ' ' ref_fname ' ' ro_fname ' ' savefname ' ' num2str(position) ' ' ...
        num2str(longestEdges) ' ' num2str(minEdgeMatch) ' ' num2str(minDotMatch) ' ' ...
        num2str(numCh) ' ' num2str(numHybs)];
    
    system(pythonCommand);

    %% get offsets from mat file     
    offsetsName = sprintf(savefnameMat, position);%['20191125_pos' num2str(position) '_offsets_initial_AdjustedPointRatio-logdots-20191213_xyse3_zse3_xyte2_zte2_xyme1_zme1.csv'];
    pointsName = sprintf('points-int-thresh-pos%d', position); % need to grab initial points
    offsetsPath = getfile(offsetsDir, offsetsName, 'match');
    pointsPath = getfile(pointsDir, pointsName, 'match');

    % run alignment code until offsets is output
    limit = 10;
    limit_iter = 1;
    while(exist(offsetsPath, 'file') ~= 2)
        % get offsets
        system(pythonCommand);
        offsetsPath = getfile(offsetsDir, offsetsName, 'match');
        limit_iter = limit_iter + 1;
        if limit_iter >= limit
            error 'alignmnet failed to output offsets';
        end
    end
    
    % get points and apply offsets
    [pointsch, offsets] = alignpointswrapper(chArray, pointsPath, offsetsPath, ...
        chaTform, numRounds, folderArray, physicalTform{position+1,1}, usechabboffsets, usephysicaloffsets);
    
    numChannels = round(numHybs / numRounds);% barcoding pseudo-channels.
    corrpoints = cell(numHybs,numCh);
    for ch = chArray % actual imaging channels.
        for ro = 1:numRounds
            for bch = 1:numChannels
                corrpoints{(ro-1)*numChannels+bch,ch} = pointsch{ch,1}{ro,1}(bch).channels;
            end
        end
    end
    

end