function [pointsorg, intensity, thresh] = mat2hybxcell(pointsDir, pointsName, chArray, numRounds, numChannels)
% converts the mat to a hyb by cell points structure
%
% assumes the points.mat file is in the 'analysis\2error-sqrt6-ch1'
% directory
%
% output threshold folder by channel matrix
% output intensity and points = hyb by channel cell array

    %T = readtable(pointsPath);
    %numHybs = max(T.hyb);
    thresh = [];
    numHybs = numRounds * numChannels;
    pointsorg = cell(numHybs, length(chArray));
    intensity = cell(numHybs, length(chArray));
    
    for c = chArray
        pointFilePath = getfile(pointsDir, pointsName, 'match');
        
        load(pointFilePath, 'points', 'threshold');
        
        for h = 1:numHybs
            r = ceil(h / numChannels);
            ch = mod(h, numChannels);
            if ch == 0
                ch = numChannels;
            end
            pointsorg{h,c} = points{r}(ch).channels;
            intensity{h,c} = points{r}(ch).intensity;
        end
        
        % convert threshold into hyb by 1 matrix
        thresh(:,c) = ch2wholethreshold(threshold);
    end


end