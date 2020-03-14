function shadingcorr = shadingcorrection(backIms)
% Get shading corrections for uneven illumination
%
% Input: background images with a cell for each channel
%
% output backcorrections using morphological operator
%
% shading corrections are divided by the image
% ex: image ./ shadingcorr
%
% Author: Sheel Shah
% Date: 2018
% Modified By: Nico Pierson
% Date: 5/7/2019

    %% Declare Variables
    numCh = length(backIms);
    numZSlices = size(backIms{1}, 3);
    shadingcorr = cell(1, numCh);
    shadingBack = cell(1, numCh);
    shadingBackMed = cell(1, numCh);

    %% Open the images with 100 radius morphological disk
    for c = 1:numCh
        for j = 1:numZSlices
            tempSlice = imopen(backIms{c}(:,:,j),strel('disk',100));
            % add values to zeros
            zeroInd = find(tempSlice == 0);
            lowMed = prctile(prctile(tempSlice,25),25);
            tempSlice(ind2sub(size(tempSlice),zeroInd)) = lowMed;
            shadingBack{c}(:,:,j) = tempSlice;
        end
    end

    % Get median
    for c = 1:numCh
        shadingBackMed{c} = median(shadingBack{c}, 3);
        shadingcorr{c} = double(shadingBackMed{c})/double(max(max(shadingBackMed{c})));
        
        figure;
        surf(double(shadingcorr{c}(1:16:end,1:16:end))),zlim([0 1]);
        ax = gca;
        ax.YDir = 'reverse';
    end
    

        
    % close all the images
    close all;

end