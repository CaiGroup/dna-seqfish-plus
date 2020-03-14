function [dots, intensity, dotsLogical] = detectdotsfindthreshold(image, regMax, logFish, threshold, varargin)
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

    rawImage = double(image); % initialize maxFish
    switch typeDetectDots
        %% Find Dots for exons or introns
        case 'exons'
        % Find Dots for exons
        dotsLogical = regMax & logFish > threshold;

        case 'exons2d'
            dotsLogical = regMax & logFish > threshold;

        case 'log'
            rawImage = logFish; % keep the raw Image
            % Find Dots for exons
            msk = true(3,3,3);
            msk(2,2,2) = false;
            apply = logFish < threshold;
            logFish(apply) = 0;
            s_dil = imdilate(logFish,msk);
            dotsLogical = logFish > s_dil; 

        case 'introns'
            rawImage = logFish; % keep the raw Image
            % Find Dots for introns: doesn't use regMax
            regMask = true(3,3,3);
            regMask(2,2,2) = false;
            % assign, to every voxel, the maximum of its neighbors
            belowThresh = logFish < threshold;
            logFish(belowThresh) = 0; % logFish needs to stay the same
            dilate = imdilate(logFish,regMask);
            % create logical matrix of dotsLogical where 1 is the voxel's value is
            % greater than its neighbors
            dotsLogical = logFish > dilate;
    end
        
                %{
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
                %}
 

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
            if intensity(i,1) > 30000 % remove very bright dots - max is 65000
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