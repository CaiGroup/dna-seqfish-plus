function directory = getdirectory(folderName, varargin)
% function finds the directory in the search space of a windows by
% searching for the main foldername provided, and uses the second
% foldername as a second check, otherwise the user provides the specified
% directory.
%
% Update: 2/7/2019
% - Fixed bug: Only finding one file, returns the file name with the
% directory. Added variable directoryOnly to only use directory.
%
% To Do:
% 1. Add options to search for a file or a directory.
% 2. So far the way to find a file is getdirectory('ij.jar', []); is there
% a cleaner way to write code for this?
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 2/5/2019



    %% Set up optional Parameters for z-slice index
    numvarargs = length(varargin);
    if numvarargs > 2
        error(['src:seqFISH-ImagePipeline:getdirectory:TooManyInputs', ...
            'requires at most 2 optional inputs']);
    end
    
    % Error for type of arguments
    if numvarargs == 1 
        if ~isempty(varargin{1})
            if ~ischar(varargin{1})
                error(['src:seqFISH-ImagePipeline:getdirectory:WrongInput', ...
                    'getdirectory varargin requires type string']);
            end
        end
    end
    if numvarargs == 2
        if ~isempty(varargin{2})
            if ~ischar(varargin{2})
                error(['src:seqFISH-ImagePipeline:getdirectory:WrongInput', ...
                    'getdirectory varargin requires type string']);
            end
        end
    end

    % set defaults for optional inputs
    optargs = {[], 'C:'};
    
    % now put these defaults into the valuesToUse cell array, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [subFolderName, startDirectory] = optargs{:};
    
    
    
        %% Check the presence of the Fiji folder in Windows - look into for
        % macs
        directory = [];
        chooseFile = false;
        searchDirectory = [];
        dirList = dir(fullfile(searchDirectory, '**', folderName));
        folderFound = false;
        if isempty(dirList)
            searchDirectory = [startDirectory filesep 'Program Files'];
            dirList = dir(fullfile(searchDirectory, '**', folderName));
        end
        if isempty(dirList)
            searchDirectory = [startDirectory filesep 'Users'];
            dirList = dir(fullfile(searchDirectory, '**', folderName));
        end
        if isempty(dirList)
            searchDirectory = [startDirectory filesep 'Program Files (86x)'];
            dirList = dir(fullfile(searchDirectory, '**', folderName));
        end
        if isempty(dirList)
            searchDirectory = startDirectory;
            dirList = dir(fullfile(searchDirectory, '**', folderName));
        end
        % Confirm existence of folder
        if ~isempty(subFolderName)
            for i = 1:length(dirList)
                if strcmp(dirList(i).name, subFolderName) % folder is scripts inside Fiji.app
                    folderFound = true;
                    break;
                end
            end
        else
            % else return directory: usually for files
            dirCount = 1;
            if ~isempty(dirList)
                % find multiple multiple directories
                numDir = size(dirList, 1);
                %if numDir > 1
                for dirname = 1:numDir
                    %'..' is the unique identifier for a directory
                    if strcmp(dirList(dirname).name, '..') 
                        % store all directories with name
                        directory{dirCount} = dirList(dirname).folder;
                        dirCount = dirCount + 1;
                    end
                end
                
                % For using folderName to find files
                if isempty(directory)
                    % [folderName] is the unique identifier for a directory
                    for dirname = 1:numDir
                        if strcmp(dirList(dirname).name, folderName) 
                            % store all directories with name
                            directory{dirCount} = [dirList(dirname).folder filesep dirList(dirname).name];
                            directoryOnly{dirCount} = dirList(dirname).folder; % use for returning only directory
                            dirCount = dirCount + 1;
                        end
                    end
                    chooseFile = true;
                end
                
                % have user choose list of directories if there are
                % multiple
                if size(directory, 2) > 1
                    list = directory;
                    directory = directory'; % transpost
                    if ~chooseFile
                        % choose directory
                        [indx,~] = listdlg('ListString', list, 'SelectionMode', ...
                            'single', 'ListSize', [400 150], 'PromptString', 'Choose one directory:');
                        directory = directory{indx};
                    else
                        % choose filename
                        %[indx,~] = listdlg('ListString', list, 'SelectionMode', ...
                        %    'single', 'ListSize', [400 150], 'PromptString', 'Choose one file:');
                        %directory = dirList(indx).folder;
                        % Choose only last filename
                        %directory = directoryOnly{length(directoryOnly)}; % get last directory in list
                        directory = directoryOnly{1}; % get first directory in list
                    end
                    
                else
                    if ~chooseFile
                        directory = directory{1}; % set equal to first directory
                    else
                        % return only director and not filename
                        directory = directoryOnly{1};
                    end
                end

                fprintf('%s file found\n\n', [directory filesep folderName]);
                return;
            else
                % user selects manually
                fprintf('%s folder directory not found\n', folderName);
                disp('Select folder directory');
                directory = uigetdir(['C:' filesep], 'Select a directory, ex: C:\Program Files\Folder1\Folder2\');
                return;
            end
            
        end
        
        
        % Checks the second folder directory
        if ~folderFound
            fprintf('Directory not found\nPlease Select a new folder directory\n');
            directory = uigetdir([dirList(1).folder filesep ], 'Select script folder directory, ex: C:\Program Files\Folder1\Folder2\');
        else
            directory = [dirList(1).folder filesep subFolderName];
        end
        
            
        

        if exist(directory, 'dir') || exist(directory, 'file')
            fprintf('%s directory found\n\n', directory);
        else
            error('%s directory was not found', directory);
        end

end