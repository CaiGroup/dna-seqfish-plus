function I = seqfish_outputpreprocessedimagesAll(experimentDir, experimentName, experimentLabel, I, offsets, chaTform, physicalTform, position, usechabboffsets, usephysicaloffsets, saveprocessedImages)
    
    numHybs = size(I,1);
    numCh = size(I,2);
    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'processedImages');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end
    
    % get the tform for channel 1, 2 and 3
    if length(chaTform) < numCh
        numCh = length(chaTform);
    end
    chabbAndPhysicalTform = cell(1, numCh);
    if usechabboffsets && usephysicaloffsets
        for ch = 1:numCh
            chabbAndPhysicalTform{ch} = addtform(chaTform{ch}, physicalTform{ch});
        end
    elseif usechabboffsets
        for ch = 1:numCh
            chabbAndPhysicalTform{ch} = chaTform{ch};
        end
    else
        error('check chromatic aberration');
    end
    
    % apply the final tform to the images
    fiduciaryOffsets = cell(numHybs, numCh);
    finalOffset = cell(numHybs, numCh);
    for f = 1:numHybs
        for ch = 1:numCh
            fiduciaryOffsets{f,ch} = maketform2(-offsets{ch}.col(f), offsets{ch}.row(f), -offsets{ch}.z(f));
            finalOffset{f,ch} = addtform(chabbAndPhysicalTform{ch}, fiduciaryOffsets{f,ch});
            I{f,ch} = imwarp(I{f,ch}, finalOffset{f,ch}, 'OutputView', imref3d(size(I{f,ch})));
        end
    end
    
    if saveprocessedImages
        % save the processed and aligned images for I (all channels).
        saveImName = ['alignedcorr-processed-I-pos' num2str(position) '-' experimentName '.mat'];
        saveImPath = fullfile(saveDirPath, saveImName);
        save(saveImPath, 'I', '-v7.3');
     
        savePath = fullfile(saveDirPath, ['alignedcorr-processed-variables-pos' num2str(position) '-' experimentName '.mat']);
        save(savePath, 'offsets', 'chaTform', 'physicalTform', 'position', 'usechabboffsets', 'usephysicaloffsets', 'finalOffset');

    end

end

