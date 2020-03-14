% script to get dots from the raw images

%fileName = 'imagesHybDapi-pos0-E14-DNA-seqFISH+rep3-1-2019-09-11.mat';
expDir = 'E:\Jonathan\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH';
%dataPath = fullfile(expDir, fileName);
%load(dataPath, 'hybIms');
%load('E:\Jonathan\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\threshold\ch1\threshold-Ch1-pos0-numHybCycles85-E14-DNA-seqFISH+rep3-1-DNAFISH.mat')

%points = cell(85,2);
typedots = 'exons';
for ch = 1:3
    listSavePath = fullfile(pwd, ['rawPointsForBeadAlignment-ch' num2str(ch) '.csv']);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID, '%s,%s,%s,%s,%s\n', 'hyb', 'x', 'y', 'z', 'int');
    for hyb = 1:85
        t = threshold(hyb,ch);
        [points, intensity, ~, ~] = detectdotsv2(hybIms{hyb, ch}, t, typedots);
        pointsSize = length(points);
        for i = 1:pointsSize
            x = points(i,1);
            y = points(i,2);
            z = points(i,3);
            int = intensity(i);
            fprintf(fileID, '%.0f,%.0f,%.0f,%.0f,%.0f\n', hyb, x, y, z, int);
        end
    end
    fclose(fileID);
end

disp('end');

% hyb, x, y, z, int