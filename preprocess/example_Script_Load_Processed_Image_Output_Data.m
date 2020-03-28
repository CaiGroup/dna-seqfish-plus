%script to load the processed Images and run the processimages.m function
%to output the final positions


experimentDir = 'I:\Limb\071919_543genes_E13.5';
experimentName = 'Limb-071919_543genes_E13_5';
channel = 3; % 488 channel
numFolders = 16;
numRounds = 4;
numChannels = 12;
sqrtradius = 4;
typedots = 'exons';
superres = 'radial';
fovArray = 0:3;
dateExp = '-2019-08-10';
totalHybs = numChannels * numRounds;

% Mike Experiment
% each channel is separate experiment
% 5 pseudochannels
% 4 rounds


% Limb Experiment
% 48 total
% 4 rounds 12 pseudochannels


% organize differently for experiments that use each channel
% 1. try without chTforms, then use if not giving good results
for position = fovArray
    
    processDataFileName = ['preProcessedData-pos' num2str(position) '-' experimentName]; %dateExp];
    processDataPath = fullfile(experimentDir, 'organizehybs', ['pos' num2str(position)], processDataFileName);
    load(processDataPath, 'I')
    Iorg = cell(totalHybs, 1);
    for f = 1:numFolders
        for ch = 1:3
            hybIndex = (f-1) * 3 + ch;
            Iorg{hybIndex} = I{f,ch};
        end
    end
    
    
    [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimages(experimentDir, experimentName, ...
        position, numRounds, numChannels, Iorg, sqrtradius, typedots, superres);
end
