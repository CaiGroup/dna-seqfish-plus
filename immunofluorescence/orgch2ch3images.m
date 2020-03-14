function Iorg = orgch2ch3images(I, numHybs, numFolders)
% organize images

    numChannels = size(I, 2);
    
    if iscell(I)
        Iorg = cell(numHybs, 1);
    else
        Iorg = zeros(numHybs, 1);
    end
    for i = 1:numChannels
        startidx = numFolders*(i-1)+1;
        endidx = startidx + numFolders - 1;
        switchIdx = [1 3 2];
        Iorg(startidx:endidx) = I(1:numFolders,switchIdx(i));
    end
end