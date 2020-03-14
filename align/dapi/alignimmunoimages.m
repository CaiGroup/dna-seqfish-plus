function [hybIms, dapiIms, tformFinal, numHybCh, numZSlice] = alignimmunoimages(experimentName, ...
    experimentDir, folderArray, position, saveDir, tformRNA2DNA, varargin)
% Aligns all images to the reference image
%
% reference image will first be aligned to reference image using the tform,
% then desample the reference image and align all hybs to this reference
% image.
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
% Date: 10/7/2019
% Author: Nico Pierson




    %% Set up optional Parameters
    argsLimit = 6;
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
    optargs = {[], true, false, false, '3d', []};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [refPath, saveData, saveHybIms, divideIms, dim, alignCh] = optargs{:};
    
    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = [datestr(dateStart, formatDate)];


    
    %% initialze variables
    numZSlice = [];
    ref = [];
    dapiIms = cell(length(folderArray), 1);
    hybIms = cell(length(folderArray), 1);
    tformFinal = cell(length(folderArray), 1);
    tformDapi = cell(length(folderArray), 1);
    
    % align RNA HybCycle0 to RNA final alignment (undersample)
    if ~isempty(refPath)
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(refPath, imageName);
        [allIms, sizeC, sizeRefZ, ~, ~] = grabimseries(imagePath, position);
        if ~isempty(alignCh)
            sizeC = alignCh;
        end
        IRNAImmuno = allIms{sizeC};
        
        % apply tform to reference
        %ref = applydapitform(ref, reftform);
        
        % desample reference
        % need to downsample by taking out 2 for each three
        idx1 = 2:3:sizeRefZ;
        idx2 = 3:3:sizeRefZ;
        idx = sort(cat(2,idx1,idx2));
        IRNAImmuno(:, :, idx) = []; % downsample
    end
    

    for folder = 1:length(folderArray)
        fprintf('Retrieving Position %.0f Folder %.0f images\n', position, folderArray(folder));



        % get the images for each channel
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(experimentDir, ['HybCycle_' num2str(folderArray(folder))], imageName);
        [allIms, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
        if sizeZ < 4
            error 'Images need 4 or more z-slices';
        end
        numDapiCh = sizeC;
        numHybCh = numDapiCh - 1;
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
        if folder ~= 1
            
            if sizeZ < 16 || divideIms
                % divide image into 4 pieces from 1 image if zslices < 16
                [newIm, ~] = imdivideby4(dapiIms{folder});
                [newImRef, numZSliceDivide] = imdivideby4(ref);
            else
                newIm = dapiIms{folder};
                newImRef = ref;
            end
            
            initialRadius = 0.0625; %0.0625 for 3d is default
            numIterations = 100; % 100 is default
            % grabtform uses imregtform 3d for zslice>16 and imregtform 2d
            % otherwise
            tformDapi{folder} = grabtform(newIm, newImRef, initialRadius, numIterations);
            tformFinal{folder} = addtform(tformRNAFinal, tformDapi{folder});
            tformFinal{folder} = addtform(tformRNA2DNA, tformFinal{folder});
            
            fprintf('Tform fov %.0f folder %.0f\n', position, folder-1);
            [tformFinal{folder}.T]
            fprintf('\n');

            % get median of the tform
            if numZSliceDivide < 16
                for c = 1:numHybCh
                    hybIms{folder, c} = imwarp(hybIms{folder, c}, tformFinal{folder}, 'OutputView', imref2d(size(hybIms{folder, c})));
                end
                dapiIms{folder} = imwarp(dapiIms{folder}, tformFinal{folder}, 'OutputView', imref2d(size(dapiIms{folder})));
            else
                tformDapiUse = tformFinal{folder};
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
                tformFinal{folder} = affine2d(eye(3));
            else
                tformFinal{folder} = affine3d(eye(4));
                tformDapi{folder} = affine3d(eye(4));
                
            end
            ref = dapiIms{1};
            tformRNAFinal = grabtform(ref, IRNAImmuno);
            tformFinal{folder} = addtform(tformRNAFinal, tformDapi{folder});% add to BELWOAOWO ELSE STATMET
            tformFinal{folder} = addtform(tformRNA2DNA, tformFinal{folder});
            tformDapiUse = tformFinal{folder};
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

    end
    
    
    
    %% save .mat files
    if saveData
        fprintf('Saving images for position %.0f\n', position);
        saveDirImages = fullfile(saveDir, ['imagesHybDapi-pos' num2str(position) ...
            '-' experimentName '-' endingDateString '.mat']);
        save(saveDirImages, 'hybIms', 'dapiIms', 'tformFinal', 'tformRNAFinal', 'tformRNA2DNA', 'tformDapi', '-v7.3');
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