% example script alignFinalRNAimages.m
% Date: 10/31/2019

% add correct packages - need to add
addpath('C:\github\streamline-seqFISH\src\AlignImages', '-end'); % to align images
addpath('C:\github\streamline-seqFISH\src\AlignImages\bfmatlab\', '-end'); % to load images
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end'); % Fiji Package to use ImageJ: Miji(false)

% variables
position = 1;
experimentName = '2019-09-04-brain-rep2-2-RNAFISH';
% ref images
fixedDir = 'I:\2019-09-09-brain-rep2-2-DNAFISH\HybCycle_0';
% moving images
movingDir = 'I:\2019-09-04-brain-rep2-2-RNAFISH\final_alignment';
saveDir = movingDir;

[polya, movingI, fixedI, tform] = alignfinalRNAimages(movingDir, fixedDir, experimentName, saveDir, position);