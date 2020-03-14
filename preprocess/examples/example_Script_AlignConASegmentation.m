% script to align final conA images

% add bfmatlab package and AlignImages package
addpath('C:\github\streamline-seqFISH\src\AlignImages\', '-end');
addpath('C:\github\streamline-seqFISH\src\AlignImages\bfmatlab', '-end');


refDir = 'I:\2019-09-09-brain-rep2-2-DNAFISH\HybCycle_0';
conADir = 'I:\2019-09-09-brain-rep2-2-DNAFISH\final_conA_segmentation';

posArray = 0:4;
tform = cell(length(posArray), 1);
dapiImages = cell(1, length(posArray) * 2);
dapiIdx = 1;


for pos = posArray
    idx = pos + 1;
    refPath = fullfile(refDir, ['MMStack_Pos' num2str(pos) '.ome.tif']);
    conAPath = fullfile(conADir, ['MMStack_Pos' num2str(pos) '.ome.tif']);
    [refIms, refSizeC, ~, ~, ~] = grabimseries(refPath, pos);
    [conAIms, conASizeC, ~, ~, ~] = grabimseries(conAPath, pos);
    
    % Save dapi images
    dapiImages{dapiIdx} = refIms{refSizeC};
    dapiImages{dapiIdx+1} = conAIms{conASizeC};
    dapiIdx = dapiIdx + 2;
    
    tform{idx} = grabtform(conAIms{conASizeC}, refIms{refSizeC});
    tform{idx}.T
    conAI = applydapitform(conAIms, tform{idx});
    dapiImages{dapiIdx+1} = conAI{conASizeC};
    savePath = fullfile(conADir, ['aligned-Pos' num2str(pos) '.tif']);
    savechannelsimagej(conAI, savePath);

end

% alignment check
saveDapiPath = fullfile(conADir, 'allDapiAlignmentCheck.tif');
savechannelsimagej(dapiImages, saveDapiPath);

% print tfroms
fileNameTform = fullfile(conADir, 'conA-offsets.csv');
printtform(tform, fileNameTform); % will print hyb as positon