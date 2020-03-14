function [I, divZSlices] = imdivideby4(image, varargin)
% function to divide the image into 4 pieces and make into a stack
%
% default division factor is 4
%
% Date 6/5/2019

    %% Set up optional Parameters for z-slice index
    argsLimit = 1;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('myfun:grabtform:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs >= 1 
        if ~isnumeric(varargin{1})
            error('myfun:imdivide:WrongInput', ...
                'initRadius requires type double');  
        end
    end
    
    % set defaults for optional inputs: 2d
    optargs = {4};
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [divFactor] = optargs{:};
    
    %% Main Function
    I1 = [];
    I2 = [];
    I3 = [];
    I4 = [];
    zSlices = size(image, 3);
    divZSlices = zSlices * divFactor;
    pixelDiv = fix(size(image,1) / (divFactor/2)); % for each axis
    % How to do this for more than one piece of image?
    for z = 1:zSlices
        image1 = image(1:pixelDiv, 1:pixelDiv, z);
        image2 = image(1:pixelDiv, pixelDiv+1:end, z);
        image3 = image(pixelDiv+1:end, 1:pixelDiv, z);
        image4 = image(pixelDiv+1:end, pixelDiv+1:end, z);
        I1 = cat(3, I1, image1);
        I2 = cat(3, I2, image2);
        I3 = cat(3, I3, image3);
        I4 = cat(3, I4, image4);
    end
    I12 = I1;
    I12 = cat(3, I12, I2);
    I34 = I3;
    I34 = cat(3, I34, I4);
    I = cat(3, I12, I34);
    clearvars I1 I2 I3 I4 I12 I34
    
end