function chTforms = getchtforms(pathName, numOfColors, varargin)
% getchtforms grabs all the chromatic aberration tforms from the tetraspec
% beads.
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/22/2019

    %% Set up optional Parameters for z-slice index
    numvarargs = length(varargin);
    if numvarargs > 3
        error('myfun:getchtforms:TooManyInputs', ...
            'requires at most 3 optional inputs');
    end
    
    % Error for type of arguments
    for i = 1:numvarargs
        if ~isnumeric(varargin{i}) || ~isscalar(varargin{i})
            error('myfun:getchtforms:WrongInput', ...
                'optional input %.0f is not numeric or scalar: requires type int', i);
        end
    end

    % set defaults for optional inputs
    optargs = {1, 0, 0}; %chTformPosition is looking at the Pos0 Image; there are 0-2 Position Images
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [reference, chTformPosition, debug] = optargs{:};

    %% Get Chromatic abberation tforms
    chTformChannels = 1:numOfColors;
    chAbFolderName = 'TetraSpec_beads_need_HyperSwap_1';
    tformPath = [pathName filesep chAbFolderName];
    tformDirectory = [pathName filesep chAbFolderName filesep 'channelsWBG'];
    mkdir(tformDirectory);
    chTforms = alltforms3D(tformPath,numOfColors,chTformPosition,chTformChannels,reference); % Reference is 1 by default; Only need one position...use other positions as a check
    if debug
        check1Position = 1;
        check1Tforms = alltforms3D(tformPath,numOfColors,check1Position,chTformChannels,reference);
        check2Position = 2;
        check2Tforms = alltforms3D(tformPath,numOfColors,check2Position,chTformChannels,reference);
        pause;
    end
    % position2 chTforms didn't match with position 0 or 1: is this a problem?
end