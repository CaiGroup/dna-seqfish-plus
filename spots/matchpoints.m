function [matchRef, mismatchRef, matchPoints, mismatchPoints] = matchpoints(pointsRef, pointsMatch, radius, intref, intch, varargin)
% Get all the points that match within a certain radius

% limiter is the maximum distance a point can be associated with another
% point

% returns structure match and mismatch with nested matrix of 'ref' and
% 'points'

%% Set up optional Parameters
    numvarargs = length(varargin);
    if numvarargs > 2
        error('myfuns:matchpoints:TooManyInputs', ...
            'requires at most 2 optional inputs');
    end
    optargs = {3}; % default of using 7 x 7 pixel grid for gaussian function
    optargs(1:numvarargs) = varargin;
    [sigmaref, sigma] = optargs{:};

    [idx,D]= rangesearch(pointsRef,pointsMatch,radius+.00001);
    [idx2,D2] = rangesearch(pointsMatch,pointsRef,radius+.00001);
    
    %limiter = 3; % matching point has to be less than 3; could be the same as the radius parameter
    [matchRef, mismatchRef] = colocalize2points(pointsRef, pointsMatch, D2, idx2, radius, intref, intch, sigma);
    [matchPoints, mismatchPoints] = colocalize2points(pointsMatch, pointsRef, D, idx, radius, intch, intref, sigmaref);
    % matchPoints.points gives the reference points without repeats and
    % matchPoints.ref gives the matching points
end

