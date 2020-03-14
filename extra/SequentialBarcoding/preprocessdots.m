function [regMax, logFish] = preprocessdots(image)
% function preprocessdots preprocesses the images with a laplacian filter
% and regional maxima logical matrix
%
% Author: Nico Pierson
% Date: 1/22/2019
% nicogpt@caltech.edu


    %% Preprocess image and get regional Maxima
    logFish = zeros(size(image));
    fish = double(image);
    for i=1:size(fish,3)
        logFish(:,:,i)=logMask(fish(:,:,i));
    end
    regMax = imregionalmax(logFish); % use variables again only once
end
