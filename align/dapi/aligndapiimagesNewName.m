function [hybIms, dapiIms, tformDapi, numHybCh, numZSlice] = aligndapiimagesNewName(experimentName, ...
    experimentDir, folderArray, position, saveDir, varargin)
% Aligns all dapi images from each set of folders to the first one, finding 
% the transformation, transformaing the image, then saving the dapi images.
%
% Inputs: splitFactorHybImages is used to split the hybIms if they are
% saved using the optional argument saveHybIms (default: false);
%
% Outputs: Saves dapi images as tiff file as 'AllHybRegistration.tif', and
% saves the data of the hybIms, dapiIms, and tform.
% Note: if folderArray is 4:6, then a 3 x channel cellarray will be output.
%
% % Dependencies: bfmatlab and Fiji.app is in the path
%
% Date: 6/2019
% Author: Nico Pierson
% Updated: 8/6/2019




    %% Set up optional Parameters
    argsLimit = 5;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:aligndapiimages:TooManyInputs', ...
            'requires at most 4 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 0
        if ~ischar(varargin{1}) 
            error('src:aligndapiimages:WrongInput', ...
                'aligndapiimages var dapiRefPath requires type string');
        end
    end
    if numvarargs > 1
        if varargin{2} ~= 0 && varargin{2} ~= 1
            error('src:aligndapiimages:WrongInput', ...
                'aligndapiimages var saveHybIms requires type int');
        end
    end
    if numvarargs > 2
        if varargin{3} ~= 0 && varargin{3} ~= 1
            error('src:aligndapiimages:WrongInput', ...
                'aligndapiimages var saveHybIms requires type int');
        end
    end
    if numvarargs > 3
        if varargin{4} ~= 0 && varargin{4} ~= 1
            error('src:aligndapiimages:WrongInput', ...
                'aligndapiimages var divideIms requires type int');
        end
    end
    % set defaults for optional inputs
    optargs = {[], true, false, false, '3d'};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [dapiRefPath, saveData, saveHybIms, divideIms, dim] = optargs{:};
    
    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = [datestr(dateStart, formatDate)];


    
    %% initialze variables
    numZSlice = [];
    dapiref = [];
    dapiIms = cell(length(folderArray), 1);
    hybIms = cell(length(folderArray), 1);
    tformDapi = cell(length(folderArray), 1);
    
    % get dapi reference if another path is specified
    if ~isempty(dapiRefPath)
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(dapiRefPath, imageName);
        [allIms, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
        dapiref = allIms{sizeC};
    end
    

    for folder = 1:length(folderArray)
        fprintf('Retrieving Position %.0f Folder %.0f images\n', position, folderArray(folder));



        % get the images for each channel
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(experimentDir, [num2str(folderArray(folder)+1)], imageName);
        [allIms, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
        if sizeZ < 4
            error 'Images need 4 or more z-slices';
        end
        numDapiCh = sizeC;
        numHybCh = sizeC - 1;
        numZSlice = sizeZ;
        divideFactor = 4;
        numZSliceDivide = sizeZ * divideFactor;
        if folder == 1
            % initialize hybIms
            hybIms = cell(length(folderArray), numHybCh);
        end
        hybIms(folder,:) = allIms(1, 1:numHybCh);
        dapiIms{folder} = allIms{numDapiCh};
        
        
        %% Get the tform and apply the transformation
        % Use the dapi transformations and the chromatic aberrations (from barcoded
        % experiments)
        if folder ~= 1 || ~isempty(dapiRefPath)% if not the first folder
            
            if sizeZ < 16 || divideIms
                % divide image into 4 pieces from 1 image if zslices < 16
                [newIm, ~] = imdivideby4(dapiIms{folder});
                [newImRef, numZSliceDivide] = imdivideby4(dapiref);
            else
                newIm = dapiIms{folder};
                newImRef = dapiref;
            end
            
            initialRadius = 0.0625; %0.0625 for 3d is default
            numIterations = 100; % 100 is default
            % grabtform uses imregtform 3d for zslice>16 and imregtform 2d
            % otherwise
            tformDapi{folder} = grabtform(newIm, newImRef, initialRadius, numIterations);
            fprintf('Tform fov %.0f folder %.0f\n', position, folder-1);
            [tformDapi{folder}.T]
            fprintf('\n');

            % get median of the tform
            if numZSliceDivide < 16
                for c = 1:numHybCh
                    hybIms{folder, c} = imwarp(hybIms{folder, c}, tformDapi{folder}, 'OutputView', imref2d(size(hybIms{folder, c})));
                end
                dapiIms{folder} = imwarp(dapiIms{folder}, tformDapi{folder}, 'OutputView', imref2d(size(dapiIms{folder})));
            else
                tformDapiUse = tformDapi{folder};
                if strcmp(dim, '2d')
                    % remove 3rd dim
                    if abs(mod(tformDapiUse.T(4,3),1)) <= 0.5
                        tformDapiUse.T(4,3) = round(tformDapiUse.T(4,3));
                    end
                end
                    
                for c = 1:numHybCh
                    hybIms{folder, c} = imwarp(hybIms{folder, c}, tformDapiUse, 'OutputView', imref3d(size(hybIms{folder, c})));
                end
                dapiIms{folder} = imwarp(dapiIms{folder}, tformDapiUse, 'OutputView', imref3d(size(dapiIms{folder})));
            end
        else
            if numZSliceDivide < 16
                tformDapi{folder} = affine2d(eye(3));
            else
                tformDapi{folder} = affine3d(eye(4));
            end
            if isempty(dapiRefPath)
                dapiref = dapiIms{1};
            end
        end      

    end
    
    
    
    %% save .mat files
    if saveData
        fprintf('Saving images for position %.0f\n', position);
        saveDirImages = fullfile(saveDir, ['imagesHybDapi-pos' num2str(position) ...
            '-' experimentName '-' endingDateString '.mat']);
        save(saveDirImages, 'hybIms', 'dapiIms', 'tformDapi', '-v7.3');
    end
    
    
    
    %% Save Dapi Images for all hyb registraction
    startStringDapiReg = 'AllHybRegistrationCheck';

    savefolchimage(position, dapiIms, saveDir, startStringDapiReg, endingDateString);
    


    %% Save the Hyb images as tif images
    if saveHybIms
        startString = 'hybImages';
        savefolchimage(position, hybIms, saveDir, startString, endingDateString)
    end
    
    
end