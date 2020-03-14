function backcorrections = getbackcorrections(imagePath, imageFileName)
% Get background for an image.
%
% Ideally use dapi to align the background to the image, and apply to the
% images. Figure out wheter it is better to correct for uneven illumination
% or to subtract the image from the other images????
%
% Assumes Dapi is the last channel
%
% Author: Sheel Shah
% Date: 2018
% Modified By: Nico Pierson
% Date: 4/1/6/2019

    % need to make a function for grabim to return the channel and the
    % zSlice number, while grabbing the whole image
    [numAllChannels, numZslices] = getimageinfo(imagePath, imageInfoName);
    numChannels = numAllChannels - 1; % subtract dapi channel
    chArray = 1:numChannels;

    I = grabim(imagePath, imageFileName);
    image = cell(1, numChannels);
    backcorrections = cell(1, numChannels);

    for ch = chArray
        startIndex = (ch - 1) * numZslices + ch;
        endIndex = startIndex + numZslices - 1;
        zIndices = startIndex:endIndex;
        image{ch} = I(:,:,zIndices);
        
        imageBack = imopen(image{ch}, strel('disk', 100));
        imageMed = median(imageBack, 3);
        normImage = double(imageMed)/double(max(max(imageMed)));
        backcorrections{ch} = normImage;
    end

end