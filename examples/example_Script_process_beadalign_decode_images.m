% script for processing yodai dna fish rep2 for each channel
%
% Dependences: 
% 1. Packages: preprocessing, processing, AlignImages
% 2. threshold: round by pseudochannel
% 3. barcodekey
% 4. (if necessary) roi segmentations
%
% Date: 8/14/2019



% Variables for processing
experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
experimentName = '2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH';
fovArray = 4;
folderArray = 0:80;
numHybs = length(folderArray) - 1; % last is a repeat
backgroundFolder = 'HybCycle_initial_background';
useBackground = false;
saveProcessIms = false;
backgroundSubtract = false;

% Variables for Decoding
numRounds = 5;
numChannels = 16;
superres = 'radial';
channelArrayBarcode = 1:2;
alloweddiff = 2;
typedots = 'exons';
sqrtradius = 6;
savePoints = true;

% Variables fo Bead Alignment
targetNumPoints = 600;
fixedThreshold = [500,500,500];
hybThreshold = [];
channelArrayBeads = 1:2;

for position = fovArray
    fprintf('Processing Position: %.0f\n\n', position);
    
    % Process the images
    [I, hybIms, tformDapi] = preprocessimages(experimentName, experimentDir, position, ...
        folderArray, useBackground, backgroundFolder, backgroundSubtract, saveProcessIms);

    %{
    % Align by beads
    [hybImsBeadAligned, hybThreshold, fixedThreshold, beadMatchPoints] = imalignbeads(experimentName, experimentDir, position, ...
        channelArrayBeads, I, hybIms, hybThreshold, fixedThreshold, tformDapi, targetNumPoints);


    % Process the images for each channel
    % Add Bead points and remove any of these points
    for ch = channelArrayBarcode
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, ...
            numfinalpoints, numpointspercell, seeds, points] = ...
            processimages(experimentDir, experimentName, position, ...
            numRounds, numChannels, hybImsBeadAligned(1:numHybs,ch), sqrtradius, typedots, ...
            superres, savePoints, alloweddiff, ch, beadMatchPoints);
    end
        %}
    
end