function [hybIms, dapiIms, tformDapiRNADNA, tformDapiRot, numZSlice] = aligndapiimagesRotglobal_v2(experimentName, ...
    folderArray, position, saveDir, hybIms, dapiIms, dapiRefPath, varargin)
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
% Date: 01/07/2020
% Author: Yodai Takei
% ytakei@caltech.edu




    %% Set up optional Parameters
    argsLimit = 2;
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
    optargs = {true, false};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [saveData, saveHybIms] = optargs{:};
    
    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    endingDateString = [datestr(dateStart, formatDate)];

    numHybCh = 2;
    
    %% initialze variables
    numZSlice = [];
    dapiref = [];

    % get dapi reference if another path is specified
    if ~isempty(dapiRefPath)
        %imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        listing = dir([dapiRefPath '\*MMStack_Pos' num2str(position) '.ome.tif']);
        if length(listing) == 1
            imageName = listing(1).name;
        else
            error 'Image names are ambiguous';
        end
        imagePath = fullfile(dapiRefPath, imageName);
        [allIms, sizeC, sizeZ, ~, ~] = grabimseries(imagePath, position);
        dapiref = allIms{sizeC};
    else
        error('set dapiRefPath');
    end
    

    for folder = 1:length(folderArray)
            
            if folder == 1
                newImRef = dapiref;%dna hyb1 dapi
                newIm = dapiIms{1};%rna hyb1 dapi
                initialRadius = 0.0625; %0.0625 for 3d is default
                numIterations = 100; % 100 is default
                tformDapiRNADNA = grabtform(newIm, newImRef, initialRadius, numIterations);
                newIm_align = imwarp(newIm, tformDapiRNADNA, 'OutputView', imref3d(size(newIm)));% align rna hyb1 dapi to dna hyb1 dapi.
                tformDapiRot = grabtformRot(newIm_align, newImRef, initialRadius, numIterations);
            end

            for c = 1:numHybCh
                hybIms{folder, c} = imwarp(hybIms{folder, c}, tformDapiRNADNA, 'OutputView', imref3d(size(hybIms{folder, c})));
                hybIms{folder, c} = imwarp(hybIms{folder, c}, tformDapiRot, 'OutputView', imref3d(size(hybIms{folder, c})));
            end
            dapiIms{folder} = imwarp(dapiIms{folder}, tformDapiRNADNA, 'OutputView', imref3d(size(dapiIms{folder})));
            dapiIms{folder} = imwarp(dapiIms{folder}, tformDapiRot, 'OutputView', imref3d(size(dapiIms{folder})));

    end
    
    
    %% save .mat files
    if saveData
        fprintf('Saving images for position %.0f\n', position);
        saveDirImages = fullfile(saveDir, ['imagesHybDapiRotCorr-pos' num2str(position) ...
            '-' experimentName '-' endingDateString '.mat']);
        save(saveDirImages, 'hybIms', 'dapiIms', 'tformDapiRot', '-v7.3');
    end
    
    dapiImsCheck{1} = dapiref;
    dapiImsCheck{2} = newIm; % before rotation correction
    dapiImsCheck{3} = dapiIms{1}; % after rotation correction
    
    %% Save Dapi Images for all hyb registraction
    startStringDapiReg = 'RotationHyb1RegistrationCheck';
    savefolchimage(position, dapiImsCheck, saveDir, startStringDapiReg, endingDateString);
    
    startStringDapiReg2 = 'RotationAllHybRegistrationCheck';
    savefolchimage(position, dapiIms, saveDir, startStringDapiReg2, endingDateString);

    %% Save the Hyb images as tif images
    if saveHybIms
        startString = 'hybImagesRotCorr';
        savefolchimage(position, hybIms, saveDir, startString, endingDateString)
    end
    
    
end