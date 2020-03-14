function [dots, intensity, dotsLogical, logImage] = detectdotsv2(image, threshold, varargin)
% Process the image and check the threshold. Find dots using a
% laplacian filter. Code is usually used for detecting dots for exons as
% they are not as bright as introns. For 'exons', threshold values are
% usually 3000 to 20,000, and the number of dots detected range from 5,000 to
% 20,000. For 'introns', threshold values are usually in the range of 200
% to 600, and there are usually 12,000 to 20,000 dots.
%
% inputs: processed image, threshold (integer), debug flag (set to 1
% for true and 0 for false).
%
% outputs: dots.channels are the x y and z coordinates, dots.intensity
% are the intensity values, and dotsFinalLogical is a logical matrix for
% each found dot.
%
% Update 1/29/2019: Add optional arguments for debug and choosing the type
% of method to detect dots. Default will be false for debug, and 'exons'
% for type of dot detection.
%
% Update 2/6/2019: Add option for more than one channel - go into for loop
% for each channels and intensity of the dots structure.
%
% Author: Sheel Shah
% Date: March 2018
% Modified By: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 2/6/2019



    %% Set up optional Parameters for z-slice index
    numvarargs = length(varargin);
    argsLimit = 4;
    if numvarargs > argsLimit
        error('src:detectdots:TooManyInputs', ...
            'requires at most 4 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs > 0 
        if ~ischar(varargin{1}) 
            error('src:detectdots:WrongInput', ...
                'detectdots var typedetectdots requires type string');
        elseif ~strcmp(varargin{1}, 'exons') && ~strcmp(varargin{1}, 'introns') ...
                && ~strcmp(varargin{1}, 'exons2d') && ~strcmp(varargin{1}, 'log') ...
                 && ~strcmp(varargin{1}, 'log2d')
            error('myfun:detectdots:WrongInput', ...
                'detectdots var typedetectdots requires type string: "exons" or "introns" or "exons2d" or "log2d"');
        end
    end
    if numvarargs > 1
        if varargin{2} ~= 0 && varargin{2} ~= 1
            error('src:detectdots:WrongInput', ...
                'detectdots var debug requires type boolean');
        end
    end
    if numvarargs > 2
        if ~ischar(varargin{3})
            error('src:detectdots:WrongInput', ...
                'detectdots var saveFigPath requires type string');
        end
    end


    % set defaults for optional inputs
    optargs = {'introns', false, '', 1};
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [typeDetectDots, debug, saveFigPath, multiplier] = optargs{:};
    
    % Initialize Variables
    threshold = threshold * multiplier;


        switch typeDetectDots
            case 'exons'
 
                %% Find Dots using Laplacian filter
                rawImage = double(image);
                logImage = zeros(size(rawImage)); % initialize logFish
                for i=1:size(rawImage, 3) % for each z-slice, apply laplacian filter
                    logImage(:,:,i) = logMask(rawImage(:,:,i));
                end
                regMax = imregionalmax(logImage); % get regional maxima
                % create logical matrix of dots that are a regional maxima and above
                % the treshold value
                dotsLogical = regMax & logImage > threshold; 
                
            case 'log'
                
                %% Log Filter developed by Sheel 9/2019
                rawImage = double(image);
                h = fspecial3('log',7,1);
                logImage = imfilter(rawImage, -20*h, 'replicate');
                msk = true(3,3,3);
                msk(2,2,2) = false;
                apply = logImage < threshold;
                logImage(apply) = 0;
                s_dil = imdilate(logImage,msk);
                dotsLogical = logImage > s_dil; 
            case 'log2d'
                
                %% Log Filter developed by Sheel 9/2019
                rawImage = double(image);
                [ys,xs,zs] = size(rawImage);
                h = fspecial('log',7,1);
                
                % initialize
                logImage = zeros(ys,xs,zs);
                dotsLogical = false(ys,xs,zs);
                
                for z = 1:zs
                    logTemp = imfilter(rawImage(:,:,z), -20*h, 'replicate');
                    msk = true(3,3);
                    msk(2,2) = false;
                    apply = logTemp < threshold;
                    logTemp(apply) = 0;
                    s_dil = imdilate(logTemp,msk);
                    dotsLogical(:,:,z) = logTemp > s_dil; 
                    logImage(:,:,z) = logTemp;
                end
                
            case 'introns'
                %% Find Dots using 3 x 3 x 3 logical mask
                rawImage = double(image); % initialize maxFish
                % set mask to 3 by 3 by 3 logical matrix, and set middle to zero
                regMask = true(3,3,3);
                regMask(2,2,2) = false;
                % assign, to every voxel, the maximum of its neighbors
                belowThresh = rawImage < threshold;
                logImage = rawImage;
                logImage(belowThresh) = 0;
                dilate = imdilate(logImage,regMask);

                % create logical matrix of dotsLogical where 1 is the voxel's value is
                % greater than its neighbors
                dotsLogical = logImage > dilate; 
            case 'exons2d'

                %% Find Dots using Laplacian filter
                rawImage = double(image);
                logImage = zeros(size(rawImage)); % initialize logFish
                for i=1:size(rawImage, 3) % for each z-slice, apply laplacian filter
                    logImageSlice = logMask(rawImage(:,:,i));

                    regMax = imregionalmax(logImageSlice); % get regional maxima
                    % create logical matrix of dots that are a regional maxima and above
                    % the treshold value
                    dotsLogical(:,:,i) = regMax & logImageSlice > threshold; 
                end

            otherwise
                error 'detectdots var typedetectdots invalid argument';
        end

        %% Remove border dots
        bordSize = 5;
        bord = ones(size(dotsLogical));
        bord(1:bordSize,:,:) = 0;
        bord(end-bordSize:end,:,:) = 0;
        bord(:,end-bordSize:end,:) = 0;
        bord(:,1:bordSize,:) = 0;
        dotsLogical = dotsLogical.*logical(bord);

        %% Debug flag to visualize image and dots
        if debug
            f = figure();
            % View image with dots: can use logFish or fish image
            imshow(max(rawImage, [], 3), 'DisplayRange', [0 2000], 'InitialMagnification', 'fit');
            hold on;
            [v2,v1] = find(max(dotsLogical, [], 3) == 1);
            scatter(v1(:), v2(:), 75);
            hold off;
            savefig(f, saveFigPath);
            close all;
        end

        %% Get x, y, z coordinates of dots
        [y,x,z] = ind2sub(size(dotsLogical), find(dotsLogical == 1));
        dots = [x y z]; % structure is location or channels

        %% Get max intensity of each dot
        %im = max(rawImage, [], 3); % max of all images: fish is not found if introns are called
        intensity = zeros(length(y), 1);
        removeDots = [];
        for i = 1:length(y)
            intensity(i,1) = rawImage(y(i), x(i), z(i));
            if intensity(i,1) > 65000 % remove very bright dots - max is 65000
                removeDots = cat(2, removeDots, i);
            end
        end
        dots(removeDots,:) = [];
        intensity(removeDots) = [];
    end
    



function lapFrame = logMask(im)   %laplacian filter 4x4 matrix

k = [-4 -1  0 -1 -4;...
     -1  2  3  2 -1;...
      0  3  4  3  0;...
     -1  2  3  2 -1;...
     -4 -1  0 -1 -4];

lapFrame = imfilter(im, k, 'repl');

end