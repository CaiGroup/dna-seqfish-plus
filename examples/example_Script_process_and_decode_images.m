% script to run all processing images and decode for each of the positions
% in channel 488
%
% Dependencies:
% preprocessing package
%
% Date: 8/9/2019

experimentDir = 'C:\Users\alignell\Desktop\Mike\Chicken\Chick_PSM\2019\76';
experimentName = 'Chick_PSM-76';
channel = 3; % 488 channel
numFolders = 20;
numRounds = 4;
numChannels = 5;
sqrtradius = 4;
typedots = 'introns';
superres = 'radial';
%fovArray = [0,1,5];
fovArray = [2,3,4,6,7,8];
folderArray = 0:20;
useBackgroundImages = false;
divideImsTransformation = false;



for position = fovArray
    % process images
    
    I = preprocessimages(experimentName, experimentDir, position, folderArray, ...
        divideImsTransformation, useBackgroundImages);

    % Get channel 3
    I488 = I(1:numFolders, channel);
    
    % decode data
    [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
        ,numpointspercell, seeds, points] = processimages(experimentDir, experimentName, ...
        position, numRounds, numChannels, I488, sqrtradius, typedots, superres);
end
