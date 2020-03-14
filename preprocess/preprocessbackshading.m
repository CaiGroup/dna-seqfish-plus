function [I, shadingcorr] = preprocessbackshading(folderArray, subtractBackground, hybIms, backIms, numCh)

    I = cell(length(folderArray), numCh);
    shadingcorr = shadingcorrection(backIms(1:numCh));
    for f = 1:length(folderArray)
        for ch = 1:numCh
            if subtractBackground % option to subract background
                imageTemp = backsubtract(hybIms{f}{ch}, backIms{ch});
            else
                imageTemp = hybIms{f}{ch};
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
            imageTemp = double(imageTemp) ./ double(shadingcorr{ch});
            I{f,ch} = uint16(imageTemp);

            if imageJBackSubtract
                % ImageJ Rolling Ball Back Subtract to remove noise using rad 3
                % replace with deconvolution - need to test first
                uniqueString = 'imageTempProcess-90jf03j';
                I{f,ch} = imagejbackgroundsubtraction(I{f,ch}, uniqueString,...
                    experimentDir);
            end

        end
    end

end