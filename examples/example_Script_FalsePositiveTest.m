% test all the parameters of the false positive for the NIH3T3 data
% 1. radius of colocalizaiton
% 2. stringency of minimum numbner of seeds
% 3. which off-target barcodes
%
% Linus' idea of removing all the dots/spots with 1 or 2 seeds
%
% 1. need to load the points and reorganize the points:
% the function was already written
%
% start with the 647 channel
% name of file: FISH_only_Pos0_647nm_results
%
% Run on the cluster instead because it is too slow


% Initial Variables
experimentDir = 'E:\Mike_2\Linus_10k_cleared_080918_NIH3T3\';
experimentName = 'Linus_10k_cleared_080918_NIH3T3';
fovArray = 0;%[17, 18, 19, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];%13:16;
channel = 647;

% Variables for Decoding
numRounds = 4;
numChannels = 20;
superres = 'radial';
alloweddiff = 1;
sqrtradius = 6;

for position = fovArray
    
    %% Retrieve the Points
    pointsFileName = ['FISH_only_Pos' num2str(position) '_' num2str(channel) 'nm_results'];
    pointsPath = fullfile(experimentDir, 'Analysis', [num2str(channel) '_Analysis'], 'Analysis_Details_NO_FISH_RCE_1.0', 'extractedData', pointsFileName);
    p = load(pointsPath, 'FISH_only');
    
    %% Reorganize the points
    points = organizeFISH_only2points(p.FISH_only);
    
    

    %% Read the Barcode excel file
    barcodeFolder = 'barcodekey';
    if isempty(channel)
        barcodekeyPath = getfile(fullfile(experimentDir,barcodeFolder), 'barcodekey');
    else
        channelFolder = ['ch' num2str(channel)];
        barcodekeyPath = getfile(fullfile(experimentDir,barcodeFolder, channelFolder), 'barcodekey');
    end
    barcodekey = readbarcode(barcodekeyPath, 'header'); % barcodekey.barcode has another column in the beginning
    barcodekey.barcode(:,1) = [];
    
    %% Decode Points for each ROI or labeled cell
    % Get the path for segmentation - need to make useful for 3d
    segmentPath = fullfile(experimentDir, ['RoiSet_Pos' num2str(position) '.zip']);
    segment = 'roi'; % can use ['roi', '3d']
    saveoption = 'list'; % can use ['matrix', 'list'] % does both, decide what to do with the options
    if isempty(channel)
        explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius)];
    else
        explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius) '-ch' num2str(channel)];
    end
    
    %% Decode
    [finalPosList, dotlocations, numpointconsensus, numdotlocations, ...
        numfinalpoints,  numpointspercell, seeds] = processallcells(experimentDir, experimentName, ...
        points, position, numRounds, numChannels, barcodekey, segmentPath, ...
        segment, sqrtradius, alloweddiff, saveoption, explabel);
    
    
end