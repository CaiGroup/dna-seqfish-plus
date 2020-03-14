function physicaltform = getphysicaloffsets(posArray, refPoints, intensity, ...
    fiduciaryType, numCh, chaTform)
% wrapper function to grab offsets from the shift in each positin and
% channel after chromatic aberration.
%
% Date: 12/13/2019

%I:\2019-07-21-E14-DNA-seqFISH+rep2-2-DNAFISH - Swapped\points\pre_formated
    physicaltform = cell(length(posArray), 1);
    for position = posArray
        % Variables
        %initialPath = fullfile(fiduciaryDir, sprintf(filename, position));
        initialPoints = cell(numCh, 1);
        initialIntensity = cell(numCh, 1);
        
        % get the file
        %initial = readtable(initialPath);

        % Grab ref points and calculate the difference in distance
        % filter out a distance more than 1
        
        %{
        idxinitial = cell(numCh, 1);
        idxinitial{1} = find(initial.ch == 1);
        x = initial.x(idxinitial{1});
        y = initial.y(idxinitial{1});
        z = initial.z(idxinitial{1});
        %}
        
        % for each channel calculate the difference and tform
        idxinitial = cell(numCh, 1);
        for ch = 1:numCh
            % apply chromatic aberration corrections
            initialPoints{ch} = transformPointsForward(chaTform{ch}, refPoints{position+1}(ch).channels);
        end
        
        indexPairs = knnsearch(initialPoints{2},initialPoints{1},'NSMethod','Exhaustive');
        matchedPoints1 = initialPoints{2}(indexPairs,:); % matching to ref
        matchedPoints2 = initialPoints{1}; % ref
        dist = sum((matchedPoints1 - matchedPoints2).^2,2);
        dist = dist< 1.5;
        matchedPoints1 = matchedPoints1(dist,:);
        matchedPoints2 = matchedPoints2(dist,:);
        for ch = 1:numCh
            if ch > 1
                initialPoints{ch} = initialPoints{ch}(indexPairs,:);
                initialIntensity{ch} = intensity{position+1}{ch}(indexPairs);
                idxinitial{ch} = indexPairs;
            else
                initialPoints{ch} = initialPoints{ch};
                initialIntensity{ch} = intensity{position+1}{ch};
                idxinitial{ch} = 1:length(indexPairs);
            end
            initialPoints{ch}(dist,:) = [];
            initialIntensity{ch}(dist) = [];
            idxinitial{ch}(dist) = [];
        end
        
        % apply the tforms to the points to the specific channels and apply
        % to the points
        physicaltform{position+1} = cell(numCh, 1);
        physicaltform{position+1}{1} = affine3d(eye(4));
        initialPointsCorrected = cell(numCh, 1);
        for ch = 2:numCh
            % get the global tform
            physicaltform{position+1}{ch} = getglobaltform(matchedPoints2,matchedPoints1);
            
            % apply the tform
            x = initialPoints{ch}(:,1);
            y = initialPoints{ch}(:,2);
            z = initialPoints{ch}(:,3);
            [x2, y2, z2] = transformPointsForward(physicaltform{position+1}{ch}, x, y, z);
            initialPointsCorrected{ch} = [x2 y2 z2];
        end
        
        %{
        % loop through table 
        for ch = 2:numCh
            for i = 1:length(idxinitial{ch})
                idx = idxinitial{ch}(i);
                pointsTemp = initialPointsCorrected{ch};
                initial.x(idx) = pointsTemp(i,1);
                initial.y(idx) = pointsTemp(i,2);
                initial.z(idx) = pointsTemp(i,3);
            end
        end
        %}
        
        % output the points
        %filesavename = fullfile(fiduciaryDir, ['ref-points-pos' num2str(position) '-' fiduciaryType '_fiducial_markers-physicaloffsets.csv']);
        %writetable(initial, filesavename);

        fprintf([fiduciaryType ' Fiduciary Physical Shift\n']);
        for ch = 2:numCh
            fprintf('physical tform pos %.0f ch %.0f:\n', position, ch);
            physicaltform{position+1}{ch}.T
        end
        
    end

end