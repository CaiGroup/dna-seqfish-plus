function [movingI, fixedI, tform] = getdapitform(movingDir, fixedDir, position)
% assumes dapi images are the last channels
%
% returns the fixed images and the aligned moving images

    % fixed images
    fixedImPath = fullfile(fixedDir, ['MMStack_Pos' num2str(position) '.ome.tif']);
    [fixedI, fixedsizeC, ~, ~, ~] = grabimseries(fixedImPath, position);
    dapiFixedI = fixedI{fixedsizeC};

    % moving images
    movImPath = fullfile(movingDir, ['MMStack_Pos' num2str(position) '.ome.tif']);
    [movingI, movsizeC, ~, ~, ~] = grabimseries(movImPath, position);
    dapiMovingI = movingI{movsizeC};

    % tform from dapi
    tform = grabtform(dapiMovingI, dapiFixedI);
    
    % apply tform
    movingI = applydapitform(movingI, tform);

end