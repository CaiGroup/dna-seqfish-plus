function [pointsch, offsets] = alignpointswrapper2(chArray, pointsPath, offsetsPath, chaTform, numRounds, folderArray)%points, experimentDir, intensity, chaTform, chArray)

    offsets = cell(length(chArray), 1);
    pointsch = cell(length(chArray), 1);

    for c = chArray
        offsetsPathCh = [offsetsPath num2str(c) '.csv'];
        offsets{c} = readtable(offsetsPathCh);
    end
    
    
   [points, intensity] = csv2hybxcell(pointsPath, chArray);
    

    
    
    for ch = chArray    
        pointsch{ch} = alignpoints(folderArray, points(:,ch), intensity(:,ch), offsets{ch}, chaTform{ch}, numRounds); % load the chabberation transformations
    end

end
