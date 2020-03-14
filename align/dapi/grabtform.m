function tform = grabtform(moving, fixed, varargin)
% grabtform returns the tform between the two images and the dapi image.
%
% Inputs: image path and part of the unique string of image
%
% All tforms should be the same.
%
% To Do List:
%
% Outputs: tform in 2d for each zSlice or in 3d for all slices
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 2/26/2019
% Modified: 

    %% Set up optional Parameters for z-slice index
    argsLimit = 2;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('myfun:grabtform:TooManyInputs', ...
            'requires at most 2 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs >= 1 
        if ~isnumeric(varargin{1})
            error('myfun:grabtform:WrongInput', ...
                'initRadius requires type double');  
        end
    end
    
    if numvarargs >= 2
        if ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
            error('myfun:grabtform:WrongInput', ...
                'numIterations requires type int');  
        end
    end
    
    % set defaults for optional inputs: 2d
    % use 0.0003 - 0.0001 for very stringent results, ball radius for
    % gradient descent; use 200-500 for stringent results for maximum
    % Iterations
    optargs = {0.0063, 100};
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [initRadius, maxIterations] = optargs{:};

% add path to CompareDotsError: addpath(['..' filesep 'CompareDotsError'']);
% need the getdirectory, and the get2dtform
  
    
    %% Check file existance
    if ~exist('moving', 'var') && ~exist('fixed', 'var')
        error 'image variables do not exist';
    end
    
    zSlices = size(moving, 3);

    if zSlices < 16
        
        moving4 = imdivideby4(moving);
        fixed4 = imdivideby4(fixed);
        
        if numvarargs <= (argsLimit - 1)
            initRadius = 0.0625;
        end
        [optimizer, metric] = imregconfig('monomodal'); % use monomodal for 3d tforms
        optimizer.MaximumStepLength = initRadius; % 0.0625 is the default for the regular step gradient descent
        optimizer.MaximumIterations = maxIterations;
        tform = imregtform(moving4, fixed4, 'translation', optimizer, metric);
        if abs(tform.T(4,3)) < 0.7
            tform.T(4,3) = 0;
        end
        
    else
        if numvarargs <= (argsLimit - 1)
            initRadius = 0.0625;
        end
        [optimizer, metric] = imregconfig('monomodal'); % use monomodal for 3d tforms
        optimizer.MaximumStepLength = initRadius; % 0.0625 is the default for the regular step gradient descent
        optimizer.MaximumIterations = maxIterations;
        tform = imregtform(moving, fixed, 'translation', optimizer, metric);
    end

end