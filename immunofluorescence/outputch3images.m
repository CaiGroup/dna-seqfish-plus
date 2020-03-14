function [] = outputch3images(experimentDir, position, folderArray, ...
        experimentName, experimentLabel, saveName, chaTform, ...
        physicaltform, I)


    % get the offsets
    saveDir = fullfile(experimentDir, 'analysis', experimentLabel);
    ch = 3;
    saveNameOut = sprintf(saveName, position);
    offsetsPath = getfile(fullfile(saveDir, 'positions'), saveNameOut, 'match');
    offsetsT = readtable(offsetsPath);
    offsets = offsetsT(offsetsT.ch == ch,:);
    
    % add offsets and save the images as processed images
    fiduciaryOffsets = cell(length(folderArray),1);
    finalOffset = cell(length(folderArray),1);
    chabbAndPhysicalTform = addtform(chaTform{ch}, physicaltform{position+1}{ch});
    idx = 1; 
    for f = folderArray
        fiduciaryOffsets{idx} = maketform2(-offsets.col(idx), offsets.row(idx), -offsets.z(idx));
        finalOffset{idx} = addtform(chabbAndPhysicalTform, fiduciaryOffsets{idx});
        % apply to each image
        I{idx} = imwarp(I{idx}, finalOffset{idx}, 'OutputView', imref3d(size(I{idx})));
        idx = idx + 1;
    end
    
    % save the processed and aligned images for ch3
    saveImName = ['aligned-processed-ch3-I-pos' num2str(position) '-' experimentLabel '-' experimentName '.mat'];
    saveImPath = fullfile(saveDir, saveImName);
    save(saveImPath, 'I', '-v7.3');
    
    % save the variables
    saveDirPath = fullfile(experimentDir, 'analysis', experimentLabel, 'points');
    if exist(saveDirPath, 'dir') ~= 7
        mkdir(saveDirPath);
    end
    %savePath = fullfile(saveDirPath, ['points-int-thresh-ch3-pos' num2str(position) '-' experimentLabel '.mat']);
    %save(savePath, 'points', 'intensity', 'adjustedThreshold', 'shadingcorr', 'sigma');

    savefolchimage(position, I(1:10), saveDir, 'alignmentcheck-', 'ch3');
    
end