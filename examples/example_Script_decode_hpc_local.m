% script to run decode all cells hpc on Linus NIH3T3 08092018 experiment
%
% folder for new threshold
% E:\Mike_2\Linus_10k_cleared_080918_NIH3T3\Analysis_BackCorrect_new_thresholds_2\488_Analysis\Analysis_Details_NO_FISH_RCE_1.0\extractedData

experimentDir = 'E:\Mike_2\Linus_10k_cleared_080918_NIH3T3';
experimentName = 'NIH3T3_080918';
position = 0;
numRounds = 4;
numChannels = 20;
cell = 5;
channel = 488;
saveDir = experimentDir;
iters = 1:3;
searchradius = sqrt(4);
err = 1;
minseeds = [2 3 4];
offpercent = 100;


% the full barcode is different for each iteration, why is this - still
% when barcodekeys are the same
%
% change randi to randperm
% change 'header to no header for decodeallcellshpc




    %% get the points
    pointsDir = fullfile(experimentDir, 'Analysis_BackCorrect_new_thresholds_2', [num2str(channel) '_Analysis'], 'Analysis_Details_NO_FISH_RCE_1.0', 'extractedData');
    pointsPath = fullfile(pointsDir, ['FISH_only_Pos' num2str(position) '_' num2str(channel) 'nm_results.mat']);
    p = load(pointsPath, 'FISH_only');
    points = organizeFISH_only2points(p.FISH_only);

for iter = iters
    decodeallcellshpc(experimentDir, experimentName, position, ...
            numRounds, numChannels, points, cell, saveDir, iter, offpercent, searchradius, err, channel, minseeds);
end