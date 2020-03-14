function I = useshadingsubtractback(experimentDir, hybIms, backIms, folderArray, numCh, shadingcorr, varargin)
% applies the given shading corrections and background subtraction if
% boolean variable is on

    %% Set up optional Parameters
    numvarargs = length(varargin);
    argsLimit = 1;
    if numvarargs > argsLimit
        error('myfuns:useshadingsubtractback:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    % set defaults for optional inputs
    optargs = {false}; 
    optargs(1:numvarargs) = varargin;
    % Place optional args in memorable variable names
    [subtractBackground] = optargs{:};
    
    
    %% Start backsubracting and apply shading corrections
    fprintf('Applying Shade Corrections for HybCycle:');
    I = cell(length(folderArray), numCh);
    for f = 1:length(folderArray)
        fprintf(' %.0f', f-1);
        for ch = 1:numCh
            
            if subtractBackground % option to subract background
                imageTemp = backsubtract(hybIms{f,ch}, backIms{ch});
            else
                imageTemp = hybIms{f,ch};
            end
            prct5 = prctile(backIms{ch}(:),5);
            % Remove Inf double values and set to percentile 5
            infInd = find(imageTemp == Inf);
            zeroInd = find(imageTemp == 0);
            if ~isempty(infInd)
                imageTemp(ind2sub(size(imageTemp),infInd)) = prct5;
            end
            if ~isempty(zeroInd)
                imageTemp(ind2sub(size(imageTemp),zeroInd)) = prct5;
            end

            % Apply the shading correctionsmean
            imageTemp = double(imageTemp) ./ shadingcorr{ch};
            I{f,ch} = uint16(imageTemp);

            % ImageJ Rolling Ball Back Subtract to remove noise using rad 3
            % replace with deconvolution - need to test first
            uniqueString = 'imageTempProcess-90jf0fa323j';
            I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString,...
                experimentDir);

        end
    end
    fprintf('\n');


end