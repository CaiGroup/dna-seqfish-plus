function pass = checkbfmatlabpath()
% checks if bfmatlab is in the MATLAB path
%
% Date: 8/6/2019
% Author: Nico Pierson
% Email: nicogpt@caltech.edu

    %% Check if Fiji is in the path
    pass = 0;
    
    try
        foundBFPath = false;
        pathCell = regexp(path, pathsep, 'split');
        for i = 1:size(pathCell, 2)
            pathCell2 = regexp(pathCell{i}, filesep, 'split');
            if strcmp(pathCell2{end}, 'bfmatlab')
                foundBFPath = true;
                fprintf('bfmatlab directory already in path\n');
                break;
            end
        end
        
        if ~foundBFPath
            % flag for now
            fprintf('Finding bfmatlab directory\n');
            bfDirectory = getdirectory('bfmatlab');
            addpath(bfDirectory, '-end'); % add to MATLAB path
        end
    catch
        error(['bfmatlab directory not found or unable to add to path' ...
            'Download package at https://downloads.openmicroscopy.org/bio-formats/5.3.4/']);
    end
    
    pass = 1;


end