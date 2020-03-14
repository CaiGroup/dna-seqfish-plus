function fijiDirectory = checkfijipath()
% checks if Fiji.app/scripts is in the MATLAB path
%
% Date: 7/9/2019
% Author: Nico Pierson
% Email: nicogpt@caltech.edu

    %% Check if Fiji is in the path
    
    try
        foundFijiPath = false;
        pathCell = regexp(path, pathsep, 'split');
        for i = 1:size(pathCell, 2)
            pathCell2 = regexp(pathCell{i}, filesep, 'split');
            if strcmp(pathCell2{end}, 'scripts') && strcmp(pathCell2{end-1}, 'Fiji.app')
                foundFijiPath = true;
                fijiDirectory = strjoin(pathCell2, filesep);
                fprintf('Fiji.app\\scripts directory already in path\n');
                break;
            end
        end
        
        if ~foundFijiPath
            % flag for now
            useOnlyFiji = true;
            if useOnlyFiji
            
                fprintf('Finding Fiji.app\\script folder\n');
                fijiDirectory = getdirectory('Fiji.app', 'scripts');
                addpath(fijiDirectory, '-end'); % add to MATLAB path
            else
                % need to download these two jar files to run imageJ in Matlab
                % download here:  http://bigwww.epfl.ch/sage/soft/mij/
                % Need both ij.jar and mij.jar for MIJ java class
                fprintf('Finding mij.jar file\n');
                %ijDirectory = [getdirectory('ij.jar', []) filesep 'ij.jar']; % choose a start directory to make it faster
                mijDirectory = [getdirectory('mij.jar', []) filesep 'mij.jar'];
                %javaaddpath(ijDirectory, '-end'); % add to java path
                javaaddpath(mijDirectory, '-end');
            end
        end
    catch
        error(['mij.jar file not found or unable to add to path' ...
            'Download package at https://downloads.openmicroscopy.org/bio-formats/4.4.9/']);
        % later add download function:
        % https://downloads.openmicroscopy.org/bio-formats/4.4.9/artifacts/bfmatlab.zip
        % or /loci_tools.jar  file
    end


end