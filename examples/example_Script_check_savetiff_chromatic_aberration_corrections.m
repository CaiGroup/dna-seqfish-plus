% script to check the chromatic aberration corrections for initial and
% final beads


% add paths
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');


%% Initialize Date for saving files
dateStart = datetime;
formatDate = 'yyyy-mm-dd';
endingDateString = datestr(dateStart, formatDate);


% variables
runLabel = 'global_chaTform_fiduciary_bead_check'; % name of folder experiment run
baseDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
finalfidDir = fullfile(baseDir, 'final_fiducial_markers');
initialfidDir = fullfile(baseDir, 'initial_fiducial_markers');
savefidDir = fullfile(baseDir, 'initial_fiducial_markers', runLabel);
pointsDirectory = fullfile(baseDir, 'points');
posArray = 0:4;

if exist(savefidDir, 'dir') ~= 7
    mkdir(savefidDir);
end

% load the chaTforms for each channel - only 2 channels
load('I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points\pre_formated\refpoints-chromaticaberrations-initial-final.mat');

for position = posArray


%% use initial beads
% load chromatic aberration corrections
%filename = ['beadInitialRefPoints_chaTforms-pos' num2str(position) '.mat'];
%chaTformFile = getfile(pointsDirectory, filename, 'match');
%load(chaTformFile, 'chaTform')

% get images
finalImPath = fullfile(finalfidDir, ['MMStack_Pos' num2str(position), '.ome.tif']);
[finalI, sizeC, sizeZ, ~,~] = grabimseries(finalImPath, position);
numCh = sizeC - 2;
% apply chaTforms
I = cell(numCh,1);
for i = 1:numCh
    I{i} = imwarp(finalI{i}, chaTformInitial{i}, 'OutputView', imref3d(size(finalI{i})));
end
% save images
savefolchimage(position, I, savefidDir, 'checkInitial2FinalBeads-Avg-chromaticaberrationcorrections', endingDateString);


%% use final beads
% load chromatic aberration corrections
filename = ['beadFinalRefPoints_chaTforms-pos' num2str(position) '.mat'];
chaTformFile = getfile(pointsDirectory, filename, 'match');
%load(chaTformFile, 'chaTform')

% get images
initialImPath = fullfile(initialfidDir, ['MMStack_Pos' num2str(position), '.ome.tif']);
[initialI, sizeC, sizeZ, ~,~] = grabimseries(initialImPath, position);
%numCh = sizeC - 1;
% apply chaTforms
I = cell(numCh,1);
for i = 1:numCh
    I{i} = imwarp(initialI{i}, chaTformFinal{i}, 'OutputView', imref3d(size(initialI{i})));
end
% save images
savefolchimage(position, I, savefidDir, 'checkFinal2InitialBeads-Avg-chromaticaberrationcorrections', endingDateString);

end
