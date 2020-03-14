function [rawpoints, intensity, pointsCsvName] = outputrawpoints(rawpoints, intensity, adjustedThreshold, ...
    folderArray, chaTform, experimentDir, experimentName, position, ...
    numChannels, shadingcorr, experimentLabel)
% output raw points to csv files


    %% save points to Csv files 
    numCh = numChannels;
    saveDir = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'pre_formated');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end
    
    
    
    %% Print points and intensity in csv file
    pointsCsvName = ['hyb-points-pos' num2str(position) '-' experimentName '-' experimentLabel '.csv'];
    listSavePath = fullfile(saveDir, pointsCsvName);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID, '%s,%s,%s,%s,%s,%s\n', 'ch', 'hyb', 'x', 'y', 'z', 'int');

    for hyb = 1:length(folderArray)
        for ch = 1:numCh
            pointsSize = length(rawpoints{hyb,ch});
            %% apply chromatic aberration corrections to the points
            %points{hyb,ch} = transformPointsForward(chaTform{ch}, rawpoints{hyb,ch});
            for i = 1:pointsSize
                x = rawpoints{hyb,ch}(i,1);
                y = rawpoints{hyb,ch}(i,2);
                z = rawpoints{hyb,ch}(i,3);
                int = intensity{hyb,ch}(i);
                fprintf(fileID, '%.0f,%.0f,%.3f,%.3f,%.3f,%.0f\n', ch, hyb, x, y, z, int);
            end
        end
    end

    pointsName = ['hyb-points-pos' num2str(position) '-' experimentName '-' experimentLabel '.mat'];
    savePath = fullfile(saveDir, pointsName);
    save(savePath, 'rawpoints', 'intensity', 'shadingcorr', 'chaTform', 'adjustedThreshold');
    
end