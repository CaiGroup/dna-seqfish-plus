% main script to run DNA seqFISH+ analysis.
% ch1-2 and ch3 analysis are integrated.
% edit by Yodai Takei
% Last updated on 01/25/20.

%% Experiment dependent Variables
experimentDir = 'F:\Yodai\DNA+\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH';
experimentName = '2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH';
experimentLabel = '2020-01-29-full-decoding';
experimentLabelNew = '2020-02-04-pipeline-full-decoding-falsepos';
%posArray = 0:4; % for decoding. needs Pos0 anyway. Don't start from Pos1 or later.
posArray = 5;

%% barcoding parameters. may need to be changed and compared.
sqrtradius = 3;
minseeds = 3;
alloweddiff = 2;

%% save global variables
dateStart = datetime;
formatDate = 'yyyy-mm-dd';
dateString = datestr(dateStart, formatDate);

%% decoding
for position = posArray
    
    loadDecodeDir = fullfile(experimentDir, 'analysis', experimentLabel, 'decoding_variables');
    loadDecodeFile = fullfile(loadDecodeDir,['decoding-variables-pos' num2str(position) '-' experimentName '.mat']);
    % remove variables that should be changed and defined above.
    load(loadDecodeFile, 'experimentDir', 'experimentName', 'experimentLabel', 'superres', ...
    'thresholdadjust', 'posArray', 'position','folderArray', 'numChAll','numCh', 'chArray', 'sizeZ', 'numRounds', ...
    'numChannels', 'usechabboffsets', 'usephysicaloffsets', ...
    'filtersigma', 'numRefPoints', 'medRefIntensity', 'longestEdges', 'minEdgeMatch', ...
    'minDotMatch', 'finalpoints', 'segment', 'minseeds', 'numCh3pointshyb', ...
    'chaTform', 'physicalTform' ,'offsets');
   
    
    % save all the variables for decoding to change the decoding paramters
    % at the next run if necessary.
    saveDecodeDir = fullfile(experimentDir, 'analysis', experimentLabelNew, 'decoding_variables');
    if exist(saveDecodeDir, 'dir') ~= 7
        mkdir(saveDecodeDir);
    end
    saveDecodePath = fullfile(saveDecodeDir,['decoding-variables-pos' num2str(position) '-' experimentName '.mat']);
    save(saveDecodePath, 'experimentDir', 'experimentName', 'experimentLabel', 'superres', ...
    'thresholdadjust', 'posArray', 'position','folderArray', 'numChAll','numCh', 'chArray', 'sizeZ', 'numRounds', ...
    'numChannels', 'sqrtradius', 'alloweddiff', 'usechabboffsets', 'usephysicaloffsets', ...
    'filtersigma', 'numRefPoints', 'medRefIntensity', 'longestEdges', 'minEdgeMatch', ...
    'minDotMatch', 'finalpoints', 'segment', 'minseeds', 'numCh3pointshyb', ...
    'chaTform', 'physicalTform' ,'offsets');
    
    for channel = chArray
        [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimagespoints(experimentDir, experimentName, ...
        position, numRounds, numChannels, finalpoints{channel}, segment, sqrtradius, alloweddiff, ...
        channel, minseeds, experimentLabelNew);
    end
    
end