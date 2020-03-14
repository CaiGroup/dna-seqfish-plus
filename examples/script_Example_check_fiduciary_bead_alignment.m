%script to check the alignment by getting the dapi and aligning it
addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');


%load('points-pos0-radial3d-alignedbeadsbyJonathan-chrabb.mat')
experimentDir = 'I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped';

% variables
folderArray = 0:9;
numChannels = length(folderArray);
position = 3;
chArray = 1:2;

allI = cell(numChannels, 1);
dapich1I = allI;
dapich2I = allI;
hybch1I = allI;
hybch2I = allI;
beadOffsetsCh1 = allI;
beadOffsetsCh2 = allI;

parfor folder = folderArray

    % get the images for each channel
    imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
    imagePath = fullfile(experimentDir, ['HybCycle_' num2str(folder)], imageName);
    [allI{folder+1}, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
    
end

offsetsDir = fullfile(experimentDir, 'points', 'positions');
offsetsName = ['20191125_pos' num2str(position) '_offsets_initial_roithresholdfiltered.csv'];
offsetsPath = getfile(offsetsDir, offsetsName, 'match');
offsets = cell(length(chArray),1);
offsetsT = readtable(offsetsPath);
for c = chArray
    offsets{c} = offsetsT(offsetsT.ch == c,:);
end


%% ch 1
ch = 1;
offsetsch1 = offsets{ch};
for folder = folderArray
    
    beadOffsetsCh1{folder+1} = maketform(offsetsch1.col(folder+1), offsetsch1.row(folder+1), offsetsch1.z(folder+1));
    tempI = applydapitform(allI{folder+1}, beadOffsetsCh1{folder+1});
    
    dapich1I{folder+1} = tempI{4};
    hybch1I{folder+1} = tempI{1};
end
saveDir = fullfile(experimentDir, 'points');
% save dapi images
startString = ['dapi-align-check-ch' num2str(ch)];
endingString = '11-6-2019';
savefolchimage(position, dapich1I, saveDir, startString, endingString);
% save hyb images
startString = ['hyb1-align-check-ch' num2str(ch)];
savefolchimage(position, hybch1I, saveDir, startString, endingString);



%% ch 2
ch = 2;
offsetsch2 = offsets{ch};

idx = 1;
for folder = folderArray
    
    beadOffsetsCh2{folder+1} = maketform(offsetsch2.col(folder+1), offsetsch2.row(folder+1), offsetsch2.z(folder+1));
    tempI = applydapitform(allI{folder+1}, beadOffsetsCh1{folder+1});
    
    dapich2I{folder+1} = tempI{4};
    hybch2I{folder+1} = tempI{1};
end


startString = ['dapi-align-check-ch' num2str(ch)];
savefolchimage(position, dapich2I, saveDir, startString, endingString);

startString = ['hyb2-align-check-ch' num2str(ch)];
savefolchimage(position, hybch2I, saveDir, startString, endingString);


disp('done');
