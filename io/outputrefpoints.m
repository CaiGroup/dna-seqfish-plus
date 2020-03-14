function [refPoints, refSaveName] = outputrefpoints(refPointsAligned, chaTform, intensity, position, endingString, saveDir, numCh)
% output final refPoints using the global chromatic aberrations for the
% other channels

    refPoints = refPointsAligned;
    refSaveName = ['ref-points-pos' num2str(position) '-' endingString '.csv'];
    listSavePath = fullfile(saveDir, refSaveName);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID, '%s,%s,%s,%s,%s\n', 'ch', 'x', 'y', 'z', 'int');
    for ch = 1:numCh
        pointsSize = length(refPoints(ch).channels);
        % apply transformation to the points
        %refPoints(ch).channels = transformPointsForward(chaTform{ch}, refPointsAligned(ch).channels);
        pointsTemp = refPoints(ch).channels;
        for i = 1:pointsSize
            x = pointsTemp(i,1);
            y = pointsTemp(i,2);
            z = pointsTemp(i,3);
            int = intensity{ch}(i);
            fprintf(fileID, '%.0f,%.3f,%.3f,%.3f,%.0f\n', ch, x, y, z, int);
        end
    end
    fclose(fileID);

end