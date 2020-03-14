function I = backsubtract(hybimage, backimage)
% function divides the image by the background image and multiplies the
% median of the background image.
%
% Used for removing autofluorescent background
%
% Date: 8/9/2019

    medianBack = median(backimage(:));
    imageTemp = (double(hybimage) ./ double(backimage));
    I = imageTemp .* double(medianBack);
end