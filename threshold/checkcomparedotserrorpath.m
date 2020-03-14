function pass = checkcomparedotserrorpath()
% checks if CompareDotsError package is in the MATLAB path
%
% Date: 7/9/2019
% Author: Nico Pierson
% Email: nicogpt@caltech.edu

    %% Check if CompareDotsError is in the path
    try
        packagePath = fullfile(pwd, '..');
        pathCell = regexp(path, pathsep, 'split');
        compareDotsFolder = 'CompareDotsError';
        onPathCompareDots = false;
        breakCompareDots = false;
        for i = 1:size(pathCell, 2)
            pathCell2 = regexp(pathCell{i}, filesep, 'split');
            if strcmp(pathCell2{end}, compareDotsFolder)
                onPathCompareDots = true;
                breakCompareDots = true;
                fprintf('CompareDotsError directory already in path\n');
            end
            
            if breakCompareDots
                break;
            end
        end
        % if not on path, add to MATLAB path
        if ~onPathCompareDots
            compareDotsDirectory = getdirectory(compareDotsFolder, [], packagePath);
            addpath(compareDotsDirectory, '-end');
        end
    catch
        error(['CompareDotsError package not found or unable to add to path' ...
            'Download package from CaiLab github or the Caltech Box account']);
    end
end