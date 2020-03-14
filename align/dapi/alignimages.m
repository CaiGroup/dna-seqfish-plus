function alignedI = alignimages(I, folderArray, offsetsPath)
% get the offsets, grab the images from the experiment and align the
% images for each channel

    %% add package
    addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');

    %% get the offset from the csv
    offsets = readtable(offsetsPath); % table of values: row_offset, col_offset, z_offset

    
    %% initialize variables
    numCh = size(I,2);
    alignedI = cell(length(folderArray),numCh);
    for folder = folderArray
        idx = folder + 1;
        for ch = 1:numCh
        
            tform = maketform(offsets.col_offset(idx), offsets.row_offset(idx), offsets.z_offset(idx));
            alignedI{idx,ch} = applydapitform(I{idx,ch}, tform);
        end
        
    end

end