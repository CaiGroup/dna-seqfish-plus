% script to check the chromatic aberration corrections for initial and
% final beads


% add paths
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');

% variables

baseDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';
finalfidDir = fullfile(baseDir, 'final_fiducial_markers');
initialfidDir = fullfile(baseDir, 'initial_fiducial_markers');
pointsDirectory = fullfile(baseDir, 'points');
posArray = 0:4;

for position = posArray


%% use initial beads
% load chromatic aberration corrections
filename = ['beadInitialRefPoints_chaTforms-pos' num2str(position) '.mat'];
chaTformFile = getfile(pointsDirectory, filename, 'match');
load(chaTformFile, 'chaTform')

% get images
finalImPath = fullfile(finalfidDir, ['MMStack_Pos' num2str(position), '.ome.tif']);
[finalI, sizeC, sizeZ, ~,~] = grabimseries(finalImPath, position);
numCh = sizeC - 1;
% apply chaTforms
I = cell(numCh,1);
for i = 1:numCh
I{i} = imwarp(finalI{i}, chaTform{i}, 'OutputView', imref3d(size(finalI{i})));
end
% save images
savefolchimage(position, I, finalfidDir, 'checkInitial2FinalBeads-Avg-chromaticaberrationcorrections', '11-20-2019');


%% use final beads
% load chromatic aberration corrections
filename = ['beadFinalRefPoints_chaTforms-pos' num2str(position) '.mat'];
chaTformFile = getfile(pointsDirectory, filename, 'match');
load(chaTformFile, 'chaTform')

% get images
initialImPath = fullfile(initialfidDir, ['MMStack_Pos' num2str(position), '.ome.tif']);
[initialI, sizeC, sizeZ, ~,~] = grabimseries(initialImPath, position);
numCh = sizeC - 1;
% apply chaTforms
I = cell(numCh,1);
for i = 1:numCh
I{i} = imwarp(initialI{i}, chaTform{i}, 'OutputView', imref3d(size(initialI{i})));
end
% save images
savefolchimage(position, I, initialfidDir, 'checkFinal2InitialBeads-Avg-chromaticaberrationcorrections', '11-20-2019');

end
