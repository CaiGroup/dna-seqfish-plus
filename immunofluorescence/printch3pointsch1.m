function saveName = printch3pointsch1(refPointsPath, experimentDir, experimentLabel, ...
    position)
% outputput channel 3 as channel 1 for alignment

    offsetsT = readtable(refPointsPath);
    channelSelect = 3;
    channelChange = 1;
    T = offsetsT(offsetsT.ch == channelSelect,:);
    T.ch(:) = channelChange;
    saveName = ['ref-ch3-points-pos' num2str(position) '-' experimentLabel '.csv'];
    savePath = fullfile(experimentDir, 'analysis', experimentLabel, 'points', 'pre_formated', saveName);
    writetable(T, savePath);

end