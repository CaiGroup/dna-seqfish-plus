function [] = savethresholdimages(experimentDir, experimentName, I, numCh, position)

       
        savePathImages = fullfile(experimentDir, 'processedimages');
        if exist(savePathImages, 'dir') ~= 7
            mkdir(savePathImages);
        end
        if numCh > 0
            I1 = I(:,1);
            save(fullfile(savePathImages, ['processedImages-' experimentName ...
                '-pos' num2str(position) '-ch' num2str(1) '.mat']), 'I1', '-v7.3');
            clearvars I1
        end
        if numCh > 1
            I2 = I(:,2);
            save(fullfile(savePathImages, ['processedImages-' experimentName ...
                '-pos' num2str(position) '-ch' num2str(3) '.mat']), 'I2', '-v7.3');
            clearvars I2
        end
        if numCh == 3
            I3 = I(:,3);
            save(fullfile(savePathImages, ['processedImages-' experimentName ...
                '-pos' num2str(position) '-ch' num2str(3) '.mat']), 'I3', '-v7.3');
            clearvars I3
        end

end