function points = alignpoints(folderArray, initialpoints, intensity, offsetsTable, chromaticAbbTforms, numRounds, physicalTforms, usechabboffsets, usephysicaloffsets)
% function to get the tforms from the csv offsets from each channel
% apply the chromatic abberations for each chanel...and save the points as
% .mat file and a csv

    %% add package
    addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');

    
    %% initialize variables
    numFolders = length(folderArray);
    numCh = round(numFolders / numRounds);
    beadOffsets = cell(numFolders, 1);
    % initialize points
    points = cell(numRounds, 1);
    for i = 1:numRounds
        points{i} = struct('channels', cell(numCh,1), 'intensity', cell(numCh, 1), 'scaledIntensity', cell(numCh,1));
    end

    for folder = folderArray
        idx = folder + 1;
        r = ceil(idx / numCh);
        ch = mod(idx, numCh);
        if ch == 0
            ch = numCh;
        end

        % offsets in python are reversed for x and z
        beadOffsets{idx} = maketform2(-offsetsTable.col(idx), offsetsTable.row(idx), -offsetsTable.z(idx));
        %beadOffsets{idx} = maketform(-offsetsTable.aligned_col(idx), offsetsTable.aligned_row(idx), -offsetsTable.aligned_z(idx));
        pointsTemp = transformPointsForward(beadOffsets{idx},  initialpoints{idx});
        if usechabboffsets
            points{r}(ch).channels = transformPointsForward(chromaticAbbTforms,  pointsTemp);
        else
            points{r}(ch).channels = pointsTemp;
        end
        if usephysicaloffsets
            points{r}(ch).channels = transformPointsForward(physicalTforms,  points{r}(ch).channels);
        end
        points{r}(ch).intensity = intensity{idx};
        points{r}(ch).scaledIntensity = intensity{idx};
        
    end

end