function [pointsch, offsets] = alignpointswrapper(chArray, pointsPath, offsetsPath, chaTform, numRounds, folderArray, physicalTforms, usechabboffsets, usephysicaloffsets)
    % add options to apply chatforms or physical shift
    % tforms....automatically uses csv files to form tforms with is already
    % corrected for chabb corrections

    %folderArray = 0:79; 
    %numRounds = 5; 
    offsets = cell(length(chArray), 1);
    pointsch = cell(length(chArray), 1);

    offsetsT = readtable(offsetsPath);
    for c = chArray
        %offsetsDir = fullfile(offsetsBaseDir, ['ch' num2str(c)], ['ch' num2str(c) '_offsets']);
        %offsetsPath = getfile(offsetsDir, ['offsets_ch' num2str(c)], 'match');
        %offsetsPath = getfile(pointsDir, ['offsets_ch' num2str(c)], 'match');
        offsets{c} = offsetsT(offsetsT.ch == c,:);
    end
    
    % set offsets equal to the first hyb, and align to ch1. as chromatic
    % aberration and physical shift are fixed and aligned to ch1.
    for c = chArray
        offsets{c}.row(:) = offsets{c}.row(:) - offsets{1}.row(1);
        offsets{c}.col(:) = offsets{c}.col(:) - offsets{1}.col(1);
        offsets{c}.z(:) = offsets{c}.z(:) - offsets{1}.z(1);
    end
    
    
    %[points, intensity] = csv2hybxcell(pointsPath, chArray);
    load(pointsPath, 'points', 'intensity')

    for ch = chArray    
        pointsch{ch} = alignpoints(folderArray, points(:,ch), intensity(:,ch), offsets{ch}, chaTform{ch}, numRounds, physicalTforms{ch}, usechabboffsets, usephysicaloffsets); % load the chabberation transformations
    end

end

% combine this piece of code to the matlab function fro
% testNormalizeIntenisty as the full pipeline to get the points, aligne the
% images using the bead images, adjust the threshold ane decode the images