function numPoints = numdetectdotsv3(image, threshold, varargin)
% Same as detectdotsv2.m except function is meant to return number of
% Points for minimizing the threshold
%
% inputs: processed image, threshold (integer), debug flag (set to 1
% for true and 0 for false).
%
% outputs: number of points
%
% To do: Look to make faster because the whole auto-threshold takes 10-12
% minutes for looking at 5 iterations, then 7 iterations at the number of
% points generated
%
% Author: Sheel Shah
% Date: March 2018
% Modified By: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 2/6/2019


    %% Set up optional Parameters for z-slice index
    numvarargs = length(varargin);
    if numvarargs > 4
        error('src:detectdots:TooManyInputs', ...
            'requires at most 4 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs == 1 
        if ~ischar(varargin{1}) 
            error('src:detectdots:WrongInput', ...
                'detectdots var typedetectdots requires type string');
        elseif ~strcmp(varargin{1}, 'exons') && ~strcmp(varargin{1}, 'introns') && ~strcmp(varargin{1}, 'exons2d')
            error('myfun:detectdots:WrongInput', ...
                'detectdots var typedetectdots requires type string: "exons" or "introns" or "exons2d"');
        end
    end
    if numvarargs == 2
        if varargin{2} ~= 0 && varargin{2} ~= 1
            error('src:detectdots:WrongInput', ...
                'detectdots var debug requires type boolean');
        end
    end
    % Set up rows or cols for images with more than one cell
    useRows = false;
    if numvarargs == 3
        rows = size(image, 1);
        cols = size(image, 2);
        if ~isnumeric(varargin{3}) || ~isscalar(varargin{3})
            error(['src:detectdots:WrongInput', ...
                'detectdots var numPseudoChannels requires type int']);
        end
        if rows > cols
            useRows = true;
        end
        if useRows
            if rows ~= varargin{3}
                error(['src:seqFISH-ImagePipeline:detectdots', ...
                    'detectdots numPseudoChannels does not match to number of cells in image']);
            end
        end
    end

    % set defaults for optional inputs
    optargs = {'introns', false, 1, ''};
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [typeDetectDots, debug, pseudoIndex, saveFigPath] = optargs{:};
    
    % Initialize Variables
    


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
            imshow(max(logImage, [], 3), [min(min(max(logImage,[],3))) mean(mean(max(logImage, [], 3))) + 5000]);
            hold on;
            [v2,v1] = find(max(dotsLogical, [], 3) == 1);
            scatter(v1(:), v2(:), 75);
            hold off;
            savefig(f, fullfile(saveFigPath, ['Dots' num2str(pseudoIndex)]));
            close all;
        end

        %% Get x, y, z coordinates of dots
        numPoints = size(find(dotsLogical == 1),1);
       
    end
    



function lapFrame = logMask(im)   %laplacian filter 4x4 matrix

k = [-4 -1  0 -1 -4;...
     -1  2  3  2 -1;...
      0  3  4  3  0;...
     -1  2  3  2 -1;...
     -4 -1  0 -1 -4];

lapFrame = imfilter(im, k, 'repl');

end