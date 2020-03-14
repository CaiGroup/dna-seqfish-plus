function threshold = thresholdbych(experimentDir, experimentName, position, ...
    folderArray, typedots, varargin)
% threshold for sequential imaging for each channel
%
% Outputs: threshold in a hyb by 1 matrix for each channel
%
% Assumptions:
% 1. Order is assumed to be in seqFISH+, where each channel is a separate
% experiment. 
%
% 2. Image names assumed to be in format: 'MMStack_Pos{number}.ome.tif'.
%
% 3. Folder names assumed to be in format: 'HybCycle_{number}'.
%
% folderCompare gets the folders with this string, strCompare gets the
% filename
%
% Dependencies: ImageJ (mij.jar file in Fiji.app/scripts path) 
%
% Author: Nico Pierson
% Date: 8/16/2019
% nicogpt@caltech.edu


    %% Set up optional Parameters
    numvarargs = length(varargin);
    argsLimit = 4;
    if numvarargs > argsLimit
        error('myfuns:thresholdbych:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    % set defaults for optional inputs
    optargs = {false, 'initial_background', [], []}; 
    optargs(1:numvarargs) = varargin;
    % Place optional args in memorable variable names
    [saveImages, backgroundFolderName, numCh, images] = optargs{:};

    
    %% Set up paths in Matlab: CompareDotsError and Fiji packages
    checkcomparedotserrorpath();
    checkfijipath();


    %% Set up variables
    options.typedetectdots = typedots; % 'exons' or 'introns'
    options.processimages = false; % do it before
    HIGH_THRESHOLD_VALUE = 99999;
    numHybCycles = length(folderArray);


    try
        if isempty(images)
            % get number of channles, and z
            imageInitPath = fullfile(experimentDir, 'HybCycle_0', 'MMStack_Pos0.ome.tif');
            positionInit = 0;
            [~, sizeC, sizeZ, ~, ~] = grabimseries(imageInitPath, positionInit);
            if isempty(numCh)
                numCh = sizeC - 1;
            end
            threshold = ones(numHybCycles, numCh) * HIGH_THRESHOLD_VALUE;
            I = cell(numHybCycles, numCh);
            points = cell(numHybCycles, numCh);
            intensity = cell(numHybCycles, numCh);
            fprintf('Preprocessing Images\n'); 
            
            image = cell(length(folderArray), 1);
            parfor f = folderArray
                fprintf('Retrieving Position %.0f Folder %.0f images\n', position, f);
                imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
                folderName = ['HybCycle_' num2str(f)]; 
                imagePath = fullfile(experimentDir, folderName, imageName);
                [image{f+1}, numDapi, sizeZ, ~, ~] = grabimseries(imagePath, position);
            end
            

            %% Uneven Illumination Correction for Background Images
            
            fprintf('Getting Background Images...\n');
            
            backImBasePath = fullfile(experimentDir, backgroundFolderName);
            backImPath = fullfile(backImBasePath, ['MMStack_Pos' num2str(position) '.ome.tif']);
            if exist(backImPath, 'file') == 2
                [backIms, numDapi, maxZ, ~, ~] = grabimseries(backImPath, position);
                % Get Shading Corrections
                shadingcorr = shadingcorrection(backIms(1:numCh));
            else
                shadingcorr = [];
            end

            %% Apply Background Subtraction and Uneven Illumination Corrections
            for f = folderArray
                fprintf('Background Subtraction HybCycle %.0f...\n', f);
                uniqueString = 'imageTemp-2903754kjfh';
                for c = 1:numCh
                    if ~isempty(shadingcorr)
                        I{f+1,c} = uint16(double(image{f+1}{c}) ./ double(shadingcorr{c}));
                    else
                        I{f+1,c} = uint16(double(image{f+1}{c}));
                    end
                    I{f+1, c} = imagejbackgroundsubtraction(I{f+1,c}, uniqueString, experimentDir);
                    
                end
                image{f+1} = []; % delete images from work space
            end
            
            
            if saveImages
                savethresholdimages(experimentDir, experimentName, I, numCh, position)
            end
            
        else
            I = images;
            numHybCycles = size(images, 1);
            if isempty(folderArray)
                folderArray = 1:numHybCycles;
            end
            if isempty(numCh)
                numCh = size(images, 2);
            end
            threshold = ones(numHybCycles, numCh) * HIGH_THRESHOLD_VALUE;
            points = cell(numHybCycles, numCh);
            intensity = cell(numHybCycles, numCh);
        end

        fprintf('Choose Threshold:\n');
        
        for c = 1:numCh
            fIter = 1;
            for f = folderArray
                fprintf('folder %.0f channel %0.f, ', f, c);
                options.channelLabel = c;
                options.folderLabel = f;
                [threshold(fIter, c), points{fIter,c}, intensity{fIter,c}] = getthreshold(I{fIter, c}, options);
                fIter = fIter + 1;
            end
            fprintf('\n');
        end
        
        % Save in the threshold directory
        thresholdDir = fullfile(experimentDir, 'threshold');
        if ~exist(thresholdDir, 'dir')
           mkdir(thresholdDir); 
        end
        savePath = fullfile(thresholdDir, ['thresholdAllCh-pos' num2str(position) '-numHybCycles' num2str(numHybCycles) '-numCh-' num2str(numCh) '-' experimentName '.mat']);
        save(savePath, 'threshold');

        % Save for each channel
        thresholdRef = threshold;
        for c = 1:numCh
            thresholdDirCh = fullfile(thresholdDir, ['ch' num2str(c)]);
            if exist(thresholdDirCh, 'dir') ~= 7
                mkdir(thresholdDirCh);
            end
            threshold = thresholdRef(:,c);
            thresholdChannel = c;
            savePathCh = fullfile(thresholdDirCh, ['threshold-Ch' num2str(c) '-pos' num2str(position) '-numHybCycles' num2str(numHybCycles) '-' experimentName '.mat']);
            save(savePathCh, 'threshold', 'thresholdChannel');
        end
        
        threshold = thresholdRef;
        
        %% save points as a csv file
        pointDir = fullfile(experimentDir, 'points');
        if exist(pointDir, 'dir') ~= 7
            mkdir(pointDir);
        end
        savePointsPath = fullfile(pointDir, ['threshold-points-intensity-pos' num2str(position) '-' experimentName '.mat']);
        save(savePointsPath, 'threshold', 'points', 'intensity'); 
        listSavePath = fullfile(pointDir, ['hyb-points-0_80-rawpoints-pos' num2str(position) '.csv']);
        fileID = fopen(listSavePath,'w');
        fprintf(fileID, '%s,%s,%s,%s,%s,%s\n', 'ch', 'hyb', 'x', 'y', 'z', 'int');
        for hyb = 1:length(folderArray)
            for ch = 1:numCh
                pointsSize = length(points{hyb,ch});
                for i = 1:pointsSize
                    x = points{hyb,ch}(i,1);
                    y = points{hyb,ch}(i,2);
                    z = points{hyb,ch}(i,3);
                    int = intensity{hyb,ch}(i);
                    fprintf(fileID, '%.0f,%.0f,%.3f,%.3f,%.3f,%.0f\n', ch, hyb, x, y, z, int);
                end
            end
        end
    
    catch err
        rethrow(err);
    end

end