% create function to look at all the thresholds in each barcode round and
% channel, get a matrix of threshold values.
%
% Things to add: Save the figures for the threshold to view??
%
% Date: 2/20/2019

position = 2; % choose position starting from 0
channel = 1; % choose channel starting from 1
pathName = ['H:' filesep 'CZI' filesep 'hybs']; % Example: 'H:\CZI\VISp'
% filename
strCompare = ['MMStack_Pos' num2str(position) '.ome.tif'];
numPseudoChannels = 12;
numBarcodes = 4;

disp('Start Thresholding');
tic
threshold = thresholdbypos(pathName, position, channel, numBarcodes, numPseudoChannels);
toc
sdisp('Done');