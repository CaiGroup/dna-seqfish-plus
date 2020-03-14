function [polya, movingI, fixedI, tform] = alignfinalRNAimages(movingDir, fixedDir, experimentName, saveDir, position)
% function to align final alignment images in brain rep2-2 RNAFISH final
% alignment to the brain rep2-2 DNAFISH dapi images, and output polya to be
% used for ilastik training.
%
% saves polya images to the polya folder in the movingDir
%
% Date: 10/31/2019

    
    % get tform 
    [movingI, fixedI, tform] = getdapitform(movingDir, fixedDir, position);


    % store dapi images 
    numFixedCh = length(fixedI);
    numMovingCh = length(movingI);
    dapiAlign = cell(2,1);
    dapiAlign{1} = fixedI{numFixedCh};
    dapiAlign{2} = movingI{numMovingCh};
    
    
    % make z-slices same 
    fixedZ = size(fixedI{1},3);
    movingZ = size(movingI{1},3);
    if fixedZ < movingZ
        addZ = movingZ - fixedZ;
        for i = fixedZ+1:fixedZ+addZ
            dapiAlign{1}(:,:,i) = zeros(2048,2048);
        end
    elseif movingz < fixedZ
        addZ = fixedZ - movingZ;
        for i = movingZ+1:movingZ+addZ
            dapiAlign{2}(:,:,i) = zeros(2048,2048);
        end
    end
    
    
    % save dapi images to compare
    startSaveString = 'dapi-alignment-check';
    endSaveString = experimentName;
    savefolchimage(position, dapiAlign, saveDir, startSaveString, endSaveString);


    % get polya channel
    polyaCh = 3;
    polya = movingI{polyaCh};
    savePolyAPath = fullfile(saveDir, ['polya-image-pos' num2str(position) '-' experimentName '.tif']);
    saveastiff(polya, savePolyAPath);

end