% Script to decode the images for 
% final output will be saved in 'analysis\2error-sqrt6'

% Set up Variables
experimentDir = 'I:\2019-09-09-brain-rep2-2-DNAFISH';
position = 0;
numRounds = 5;
numChannels = 16;
sqrtradius = 6;
segment = 'whole';
typedots = 'exons';
superres = 'radial';
experimentName = '2019-09-09-brain-rep2-2-DNAFISH'; 
minseeds = 3;
alloweddiff = 2;

% To load the processed images
%dataDir = fullfile(experimentDir, 'processedimages\pos0\preProcessedData-pos0-2019-09-09-brain-rep2-2-DNAFISH-2019-09-28.mat');
%load(dataDir, 'I');


for channel = 2%1:2
    [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
    ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
    position, numRounds, numChannels, pointsch{channel}, segment, sqrtradius, typedots, superres, alloweddiff, ...
    channel, minseeds);
end

