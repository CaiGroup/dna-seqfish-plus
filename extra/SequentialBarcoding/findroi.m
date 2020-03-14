function roiNumber = findroi(roiPath, vertex, varargin)
% Find if the corresponding roi centroid inside the polygon points vertex(X
% and Y), specifically used to match the cytoplasm roi with the nucleus
% roi. If not found, set the default value to null.
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/23/2019

    %% Set up optional Parameters
    numvarargs = length(varargin);
    if numvarargs > 1
        error('myfun:findroi:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs == 1
        if ~isnumeric(varargin{1}) || ~isscalar(varargin{1})
            error('myfun:findroi:WrongInput', ...
                'total pixels is not numeric or scalar: requires type int');
        elseif ~(mod(varargin{1}, 1) == 0)
            error('myfun:findroi:WrongInput', ...
                'total pixels is not type integer: requires type int');
        end
    end
    
    % set defaults for optional inputs
    optargs = {2048};
    
    % now put these defaults, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [totPixels] = optargs{:};
    
    %% Find the Centroid of the roi with the PathName given
    vertexCentroid = selfseg(roiPath);
    for roi = 1:length(vertexCentroid)
        BWnucl = poly2mask(vertexCentroid(roi).x, vertexCentroid(roi).y, totPixels, totPixels);
        %area_cyto = sum(sum(BWnucl));
        nucl = regionprops(BWnucl,'centroid');
        xCentroid = nucl(1).Centroid(1);
        yCentroid = nucl(1).Centroid(2);

        %% Return number of roi if found
        include = inpolygon(xCentroid, yCentroid, vertex.x, vertex.y);
        if include
            roiNumber = roi;
            return; 
        end
    end
    
    % If not found return null
    roiNumber = [];
    return;
    
end