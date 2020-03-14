function bfmatlabDirectory = checkbfmatlabpath()
% checks if FindThreshold\bfmatlab is in the MATLAB path
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
            if strcmp(pathCell2{end}, 'bfmatlab') && strcmp(pathCell2{end-1}, 'FindThreshold')
                foundFijiPath = true;
                bfmatlabDirectory = strjoin(pathCell2, filesep);
                fprintf('FindThreshold\\bfmatlab directory already in path\n');
                break;
            end
        end
        
        if ~foundFijiPath
            
            lastwarn('');    
            bfPathAdd = fullfile('C:', 'github', 'streamline-seqFISH', 'src', 'FindThreshold', 'bfmatlab');
            addpath(bfPathAdd, '-end');
            [warnMsg, warnId] = lastwarn;
            
            % flag for now
            if ~isempty(warnMsg) % if warning not empty then find directory
            
                fprintf('Finding FindThreshold\\bfmatlab folder\n');
                bfmatlabDirectory = getdirectory('FindThreshold', 'bfmatlab');
                addpath(bfmatlabDirectory, '-end'); % add to MATLAB path
            end
        end
    catch
        error(['bfmatlab package not found or unable to add to path' ...
            'https://downloads.openmicroscopy.org/bio-formats/4.4.9/artifacts/bfmatlab.zip']);
        % later add download function:
        % https://downloads.openmicroscopy.org/bio-formats/4.4.9/artifacts/bfmatlab.zip
        % or /loci_tools.jar  file
    end


end