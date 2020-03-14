function [] = organizesequentialimage(pathName, folderArray, position, varargin)
% organizes the images used for sequential barcoding.
% Make sure the images are in the correct order of xyczt; if not hyperswap
% the images.
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/23/2019

% Examples:
% folderArray = 22:28;
% pathName = 'D:\MATLAB\CaiLab\Development\YodaiIntronRep4\Hybs';
% folderName = 'TestOrganized';
% position = 0;

    %% Set up optional Parameters
    numvarargs = length(varargin);
    if numvarargs > 1
        error('myfun:organizesequentialimage:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs == 1
        if ~ischar(varargin{1})
            error('myfun:organizesequentialimage:WrongInput', ...
                'imput is not a string: requires type char');
        end
    end
    
    % set defaults for optional inputs
    optargs = {'Organized'};
    
    % now put these defaults, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [folderName] = optargs{:};

    %% Start Function
    fprintf('Starting organizesequentialimage.m\n\n');
    tic
    
    %% Make the new directory for the new images
    basePath = [pathName filesep '..'];
    newDirectory = [basePath filesep folderName filesep 'Pos' num2str(position)];
    mkdir(newDirectory);
    
    %% Move Images running through Folders
    for folder = folderArray
        fprintf('Looking in folder %.0f\n', folder);
        source = [pathName filesep num2str(folder)];
        destination = newDirectory;
        newFileName = [num2str(folder) '.tif'];
        moveimage(source, destination, position, newFileName)
    end
    
    toc
    fprintf('Finished copyting images\n');

end