function medianIntensity = outputmedianintensity(pointsDir, pointsName, chArray, numRounds, numChannels)
% get median intensity for thresholded reference points

    %% reference median intensity for capture points in each hybridization
    %pointsRefPath = 'I:\2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH\points\hyb-points-0_80-ForBeadAlignment-pos0-radialcenter3d-2019-07-29-E14-DNA-seqFISH+rep3-1-DNAFISH.csv';
    
    %[refPoints, refIntensity] = csv2hybxcell(pointsRefPath, chArray); %
    %for csv values

    % for mat file with points
    [refPoints, refIntensity, threshold] = mat2hybxcell(pointsDir, pointsName, chArray, numRounds, numChannels);


    % get the median ref intensity
    medianIntensity = refIntensity;
    for i = 1:size(medianIntensity,1)
        for ch = chArray
            medianIntensity{i,ch} = median(medianIntensity{i,ch});
        end
    end

end