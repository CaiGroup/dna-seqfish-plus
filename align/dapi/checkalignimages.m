function [movingImage, fixedImage] = checkalignimages(moving, fixed, dapi1, dapi2, zSlices, savePath)    
% checkalignimages aligns the images with the tforms and return the images,
% while also saving the fixed and moving images together to compare
% alignment.
%
% Return the tformImage, and the refImage as cells with each channel
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 2/26/2019
% Modified:

    %% Check if Fiji is in the path
    try
        foundFijiPath = false;
        pathCell = regexp(path, pathsep, 'split');
        for i = 1:size(pathCell, 2)
            pathCell2 = regexp(pathCell{i}, filesep, 'split');
            if strcmp(pathCell2{end}, 'scripts') && strcmp(pathCell2{end-1}, 'Fiji.app')
                foundFijiPath = true;
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
                fprintf('Finding ij.jar, and mij.jar files\n\n');
                ijDirectory = [getdirectory('ij.jar', []) filesep 'ij.jar']; % choose a start directory to make it faster
                mijDirectory = [getdirectory('mij.jar', []) filesep 'mij.jar'];
                javaaddpath(ijDirectory, '-end'); % add to java path
                javaaddpath(mijDirectory, '-end');
            end
        end
    catch
        error(['bfmatlab package not found or unable to add to path' ...
            'Download package at https://downloads.openmicroscopy.org/bio-formats/4.4.9/']);
        % later add download function:
        % https://downloads.openmicroscopy.org/bio-formats/4.4.9/artifacts/bfmatlab.zip
        % or /loci_tools.jar  file
    end

    % Variables
    numZslices = zSlices;
    numChannels = ceil(size(fixed, 3) / numZslices);

    % How to include the number of channels as well???
    movingImage = cell(1, numChannels);
    fixedImage = cell(1, numChannels);
    for i = 1:size(fixed, 3) % loop over each zSlice
        imageIndex = ceil(i / numZslices);
        movingImage{imageIndex} = cat(3, movingImage{imageIndex}, moving(:,:,i));
        fixedImage{imageIndex} = cat(3, fixedImage{imageIndex}, fixed(:,:,i));
    end


    %% Save Images using tformImage and refImage
    Miji(false);
    
    for channel = 1:numChannels
        
        loop = 1;
        name1 = ['C' num2str(1) '-'  num2str(loop) '.tif'];
        name2 = ['C' num2str(2) '-'  num2str(loop) '.tif'];
        MIJ.createImage(name1, movingImage{channel}, true);
        MIJ.createImage(name2, fixedImage{channel}, true);

        str = [];
        str = [str ' image' num2str(1) '=C' num2str(1) '-' num2str(loop) '.tif'];
        str = [str ' image' num2str(2) '=C' num2str(2) '-' num2str(loop) '.tif'];


        saveFileName = ['FirstHybToLastHyb_AlignCheck_Ch_' num2str(channel) '.tif'];
        saveFinalPath = fullfile(savePath, saveFileName);
        try
            MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(2) ' slices=' num2str(numZslices) ' frames=1 display=Grayscale']);
            MIJ.run('Save', ['save=[' saveFinalPath ']']);
            MIJ.run('Close All')
        catch
            MIJ.closeAllWindows;
            MIJ.exit;
            error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
        end
        
    end
    
    
    %% save the dapi images as two channels
    loop = 1;
    name1 = ['C' num2str(1) '-'  num2str(loop) '.tif'];
    name2 = ['C' num2str(2) '-'  num2str(loop) '.tif'];
    MIJ.createImage(name1, dapi1, true);
    MIJ.createImage(name2, dapi2, true);

    str = [];
    str = [str ' image' num2str(1) '=C' num2str(1) '-' num2str(loop) '.tif'];
    str = [str ' image' num2str(2) '=C' num2str(2) '-' num2str(loop) '.tif'];


    saveDapiFileName = 'FirstHybToLastHyb_DapiAlignCheck.tif';
    saveDapiFinalPath = fullfile(savePath, saveDapiFileName);
    try
        MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
        MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(2) ' slices=' num2str(numZslices) ' frames=1 display=Grayscale']);
        MIJ.run('Save', ['save=[' saveDapiFinalPath ']']);
        MIJ.run('Close All')
    catch
        error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
    end
    
    
    
    
    
    %% debug to view images
    viewImage = false;
    if viewImage
        % get the max z-projection and view the image
        close all;
        addHigherPixelRange = 1500;
        figure;
        imshow(max(moving, [], 3), 'DisplayRange', [min(min(max(moving,[],3))) mean(mean(max(moving,[],3))) + addHigherPixelRange], 'InitialMagnification', 'fit');
        pause;
        close all;
    end
    
    %warning ('on','all');
    
end