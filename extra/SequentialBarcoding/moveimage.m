function [] = moveimage(source, destination, position, varargin)
% moveimage moves an image from its source to the new destination
% with optional arguments of a new name and a compare string to find the
% correct image. Default for the compare string is ['Pos' (position) '.ome']
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/23/2019

    %% Set up optional Parameters
    numvarargs = length(varargin);
    if numvarargs > 2
        error('myfun:moveimage:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs == 1
        if ~ischar(varargin{1})
            error('myfun:moveimage:WrongInput', ...
                'imput is not a string: requires type char');
        end
    end
    
    if numvarargs == 2
        if ~ischar(varargin{2})
            error('myfun:moveimage:WrongInput', ...
                'imput is not a string: requires type char');
        end
    end
    
    % set defaults for optional inputs
    defaultStrCompare = ['Pos' num2str(position) '.ome'];
    optargs = {[], defaultStrCompare};
    
    % now put these defaults,
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [newName, strCompare] = optargs{:};

    %% Find Image using compare string
    allFiles = dir(source);
    % Finds file with position number
    for i = 1:length(allFiles)
        if ~isempty(strfind(allFiles(i).name, strCompare)) % Only need this if loop because strCompare is passed through with position number
            myFile = allFiles(i).name;
            break
        end
        % send error message if file is not found
        if i == length(allFiles)
            error('file: %s not found', source);
        end
    end
    
    %% Move Image to destination
    sourcePath = [source filesep myFile];
    if isempty(newName) % set new file name to same if null
        newName = myFile;
    end
    % Set up destination path
    destinationPath = [destination filesep newName];
    status = copyfile(sourcePath, destinationPath);
  
    % Check if successful
    if status
        fprintf('File: %s copied to \n%s\n\n', myFile, destinationPath);
    end

end