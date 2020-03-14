function [numRefPoints, medIntensity] = loadrefpoints(pointsDir, pointsName, chArray, numRounds, numChannels, varargin)
% loads reference points and reference threshold used 
%
% assumes the points.mat file is in the 'analysis\points\hyb80 pos0'
% directory
%
% output threshold folder by channel matrix
% output intensity and points = hyb by channel cell array


    %% Set up optional Parameters
    argsLimit = 1;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:adjustthreshold:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end   
    % set defaults for optional inputs
    optargs = {[]};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [roimask] = optargs{:};

    numHybs = numRounds * numChannels;
    numRefPoints = cell(numHybs, length(chArray));
    medIntensity = cell(numHybs, length(chArray));
    sizeZ = size(roimask,3);
    max = 2048;
    min = 1;
    
    pointFilePath = getfile(pointsDir, pointsName, 'match');
        
    load(pointFilePath, 'points', 'intensity');
        
    for c = chArray
        for h = 1:numHybs
            removeInd = [];
            if ~isempty(roimask)
                pointsdata = round(points{h,c});
                numPoints = size(pointsdata,1);
                for i = 1:numPoints
                    x = pointsdata(i,1);
                    y = pointsdata(i,2);
                    z = pointsdata(i,3);
                    if x < min
                        x = min;
                    elseif x > max
                        x = max;
                    end
                    if y < min
                        y = 1;
                    elseif y > max
                        y = max;
                    end
                    if z < min
                        z = min;
                    elseif z > sizeZ
                        z = sizeZ;
                    elseif isnan(z)
                        z = round(pointsdata(i-1,3));
                    end
                    if ~roimask(y, x, z)
                        removeInd = cat(1, removeInd, i);
                    end
                end
                points{h,c}(removeInd,:) = [];
                intensity{h,c}(removeInd) = [];
            end
            
            numRefPoints{h,c} = size(points{h,c},1);
            medIntensity{h,c} = median(intensity{h,c});
        end

    end


end