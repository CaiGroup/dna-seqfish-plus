function [] = frequencyhist(points, varargin)
% function makes histogram with the points and number of Intervals
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 3/20/2019

    %% Set up optional Parameters
    
    numvarargs = length(varargin);
    if numvarargs > 1
        error('myfuns:getallbeads:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end

    % set defaults for optional inputs
    optargs = {100}; 
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Place optional args in memorable variable names
    [numIntervals] = optargs{:};

    intervalWidth = (max(points) - min(points))/numIntervals;
    x = 0:intervalWidth:max(points);
    ncount = histc(points,x);
    relativefreq = ncount/length(points);
    figure
    bar(x-intervalWidth/2, relativefreq,1)
    xlim([min(x) max(x)])
    %set(gca, 'xtick', x) % too many ticks

end
