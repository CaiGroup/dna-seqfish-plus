function [regMax, logImage] = preprocessdots(image, varargin)
% function preprocessdots preprocesses the images with a laplacian filter
% and regional maxima logical matrix
%
% Updated 1/29/2019: Set up option to use find dots for introns or exons.
%
% Author: Nico Pierson
% Date: 1/22/2019
% nicogpt@caltech.edu
% Modified: 1/29/2019

    %% Set up optional Parameters for z-slice index
    numvarargs = length(varargin);
    if numvarargs > 1
        error('src:preprocess:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs == 1 
        if ~ischar(varargin{1}) 
            error('myfun:grabim:WrongInput', ...
                'image organization requires type string');
        elseif ~strcmp(varargin{1}, 'exons') && ~strcmp(varargin{1}, 'introns') && ~strcmp(varargin{1}, 'exons2d') && ~strcmp(varargin{1}, 'log')
            error('myfun:grabim:WrongInput', ...
                'image organization requires type string: "exons" or "introns" or "exons2d"');
        end
    end

    % set defaults for optional inputs
    optargs = {'exons'}; % options are 'introns' or 'exons'
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [typeDetectDots] = optargs{:};


    %% Preprocess image and get regional Maxima
    logImage = [];
    regMax = [];
    switch typeDetectDots
        case 'exons'
            %% Find Dots using a laplacian filter
            logImage = zeros(size(image));
            fish = double(image);
            for i=1:size(fish,3)
                logImage(:,:,i)=logMask(fish(:,:,i));
            end
            regMax = imregionalmax(logImage); % use variables again only once
            
        case 'log'

            %% Log Filter developed by Sheel 9/2019
            rawImage = double(image);
            h = fspecial3('log',7,1);
            logImage = imfilter(rawImage, -20*h, 'replicate');
            regMax = [];
                

        case 'exons2d'
            %% Find Dots using Laplacian filter
            rawImage = double(image);
            logImage = zeros(size(rawImage)); % initialize logFish
            for i=1:size(rawImage, 3) % for each z-slice, apply laplacian filter
                logImage(:,:,i) = logMask(rawImage(:,:,i));
                regMax(:,:,i) = imregionalmax(logImage(:,:,i)); % get regional maxima
            end
            
        case 'introns'
            %% Find Dots using 3 x 3 x 3 logical mask
            logImage = double(image); % initialize maxFish
            regMax = [];
            
        otherwise
            error 'Invalid optional argument used';
    end
end
