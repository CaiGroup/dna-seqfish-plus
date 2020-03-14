conADir = 'I:\2019-09-09-brain-rep2-2-DNAFISH\final_conA_segmentation';
posArray = 0:4;
dapiIdx=1;
dapiImages = cell(1, length(posArray) * 2);
for pos = posArray
    idx = pos + 1;
    conAPath = fullfile(conADir, ['MMStack_Pos' num2str(pos) '.ome.tif']);
    [conAIms, conASizeC, ~, ~, ~] = grabimseries(conAPath, pos);
    conAI = applydapitform(conAIms, tform{idx});
    dapiImages{dapiIdx+1} = conAI{conASizeC};
    dapiIdx = dapiIdx + 2;
end

saveDapiPath = fullfile(conADir, 'allDapiAlignmentCheck.tif');
savechannelsimagej(dapiImages, saveDapiPath);