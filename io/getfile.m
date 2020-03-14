function filedirectory = getfile(directory, filename, varargin)
% function finds files in the directory and lists them for the user to
% choose
%
% Update: 2/6/2019
% - Added option to use 'select' (default) to choose a list or 'match' to 
% get a specific filename.
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 2/6/2019

    %% Set up optional Parameters for z-slice index
    argsLimit = 1;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:getfile:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs == 1 
        if ~ischar(varargin{1}) 
            error('myfun:getfile:WrongInput', ...
                'image organization requires type string');
        elseif ~strcmp(varargin{1}, 'select') && ~strcmp(varargin{1}, 'match') && ~strcmp(varargin{1}, 'strict')
            error('myfun:getfile:WrongInput', ...
                'image organization requires type string: "select", "match", or "strict"');
        end
    end

    % set defaults for optional inputs
    optargs = {'select'}; % options 'select' or 'match'
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [searchType] = optargs{:};

    
    
    %% Get file based on selection or by filename
    
    % Variables
    dirList = dir(directory);
    dirLength = size(dirList, 1);
    dirCount = 1;
    filelist = [];
    
    switch searchType
        case 'select'
            for i = 1:dirLength
                if ~strcmp(dirList(i).name, '.') && ~strcmp(dirList(i).name, '..') && isfile([dirList(i).folder filesep dirList(i).name])
                    filelist{dirCount} = dirList(i).name;
                    dirCount = dirCount + 1;
                end
            end

            if ~isempty(filelist)
                if size(filelist, 2) == 1
                    filedirectory = [directory filesep filelist{1}];
                else
                    [indx,~] = listdlg('ListString', filelist, 'SelectionMode', ...
                        'single', 'ListSize', [400 150], 'PromptString', ['Choose one file for the ' filename ':']);
                    filedirectory = [directory filesep filelist{indx}];
                end
            else
                % set filedirectory to null
                filedirectory = [];
            end
        case 'match'
            for i = 1:dirLength
                if ~isempty(strfind(dirList(i).name, filename)) && exist([dirList(i).folder filesep dirList(i).name], 'file') == 2
                    filelist{dirCount} = dirList(i).name;
                    dirCount = dirCount + 1;
                end
            end
            if ~isempty(filelist)
                if size(filelist, 2) == 1
                    filedirectory = [directory filesep filelist{1}];
                else
                    % Select if multiple files
                    %[indx,~] = listdlg('ListString', filelist, 'SelectionMode', ...
                    %    'single', 'ListSize', [400 150], 'PromptString', ['Choose one file for the ' filename ':']);
                    %filedirectory = [directory filesep filelist{indx}];
                    filedirectory = [directory filesep filelist{1}]; % get the first one or last?
                end
            else
                % set filedirectory to null
                filedirectory = [];
            end
        case 'strict' % get the exact filename
            lengthStr = length(filename);
            for i = 1:dirLength
                if strncmp(dirList(i).name, filename, lengthStr) && isfile([dirList(i).folder filesep dirList(i).name])
                    filelist{dirCount} = dirList(i).name;
                    dirCount = dirCount + 1;
                end
            end
            if ~isempty(filelist)
                if size(filelist, 2) == 1
                    filedirectory = [directory filesep filelist{1}];
                else
                    filedirectory = [directory filesep filelist{1}]; % get the first one or last?
                end
            else
                % set filedirectory to null
                filedirectory = [];
            end
            
        otherwise
            error 'error in getfile.m';
    end


end