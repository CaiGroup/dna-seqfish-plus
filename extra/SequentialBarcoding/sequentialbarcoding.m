function [] = sequentialbarcoding(pathName, folderArray, numOfChannels, fieldOfView, varargin)
% function used for single gene barcoding for each position or field of
% view. After processing all images in specified position or field of view,
% threshold values will be assigned via dialog prompts.
%
% Notice: Make sure the Tetraspec Bead photos, Blanks, and Hyb photos are
% Hyperswapped or in the order xyczt.
%
% Default folders are 'Hybs' for the raw images, 'Organized 
%
% Things to Do: 1. Combine all the variables into a few mat files
% 2. Turn this script into a simple function
% 3. I believe I have to process all images for every position before
% assigning thresholding
% 4. Change ImageJ directories for Fiji.app/scripts folder - addpath
%
% Example Variables:     
%    startFolder = 22;
%    endFolder = 28;
%    folderArray = startFolder:endFolder;
%    numOfChannels = 3;
%    fieldOfView = 0;
%    currDirectory = pwd;
%    pathName = [currDirectory filesep '..' filesep '..' filesep 'YodaiIntronRep4'];
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/23/2019

    %% Set up optional Parameters
    numvarargs = length(varargin);
    if numvarargs > 3
        error('myfun:sequentialbarcoding:TooManyInputs', ...
            'requires at most 3 optional inputs');
    end
    
    % Error for type of arguments
    if numvarargs > 0
        for i = 1:numvarargs
            if ~ischar(varargin{i})
                error('myfun:sequentialbarcoding:WrongInput', ...
                    'imput is not a string: requires type char');
            end
        end
    end
    
    % set defaults for optional inputs
    optargs = {'Hybs', 'Organized', 'SequentialData'};
    
    % now put these defaults, 
    % and overwrite the ones specified in varargin.
    optargs(1:numvarargs) = varargin;
    
    % Default Value of ref image is 1
    [rawImageFolder, organizedFolder, seqDataFolder] = optargs{:};

    
    %% Set DIARY FileName and Initialize Date
    % Initialize Date
    dateStart = datetime;
    formatDate = 'yyyymmdd';
    dateSaveString = datestr(dateStart, formatDate);
    dateStartString = datestr(dateStart);
    % Set diaryFileName to same
    diary 'testSequentialBarcodingPos0_IntronsRep4_Jan23_2019.txt' % comment out diary later
    functionName = 'scriptSequentialBarcoding.m';
    fprintf('**************************************************************************\n\n'); 
    fprintf('Starting function %s on %s\n\n', functionName, dateStartString); 

    
    %% Set up path for ImageJ
    % find mij.jar and ij.jar files and add to the java path
    javaaddpath 'C:\Program Files\MATLAB\R2018a\java\mij.jar'
    javaaddpath 'C:\Program Files\MATLAB\R2018a\java\ij.jar'

    
    %% Set up Directories from the main path
    imagePath = [pathName filesep rawImageFolder];
    compilePath = [pathName filesep organizedFolder];
    organizedImageDir = [compilePath filesep 'Pos' num2str(fieldOfView)];
    savePath = [organizedImageDir filesep seqDataFolder];
    mkdir(savePath);


    %% Move the Images (Organize) to Different Folders and Rename
    organizeImage = 0; % true
    for i = folderArray % check to see if files exist
        source = [organizedImageDir filesep num2str(i)];
        if ~exist(source, 'file')
            organizeImage = 1; % set to false if file doesn't exist
            break;
        end
    end
    if organizeImage
        organizesequentialimage(imagePath, folderArray, fieldOfView)
    end

    
    %% Load preprocess Data Or Calculate Chromatic Aberration Tforms and Background Corrections
    % Load data
    processDataPath = [organizedImageDir filesep 'Pos' num2str(fieldOfView) '_FirstCorrections'];
    if exist(processDataPath, 'file')
        processdata = load(processDataPath);
    else
        processdata = [];
    end

    if isempty(processdata)
        % Get chromatic aberration tforms
        chTforms = getchtforms(pathName, numOfChannels);
        % Get the back corrections
        backcorrections = getbackcorrections(pathName, numOfChannels);
        % add later:  dateSaveString
        save([savePath filesep dateSaveString 'Pos' num2str(fieldOfView) '_FirstCorrections'], 'backcorrections', 'chTforms');
    else
        chTforms = processdata.chTforms;
        backcorrections = processdata.backcorrections;
    end


    %% Compile the Images 
    channelsall = cell(1, length(folderArray));
    registration = ones(1,length(folderArray)) * 4;
    for i = 1:length(folderArray)
       channelsall{i} = 1:numOfChannels; 
    end
    [prehybnum, regis] = CompileImages(compilePath, fieldOfView, folderArray, channelsall, registration, backcorrections, chTforms);


    %% Register the images to the first Hyb1 Dapi Image
    % convert regimage to uint16 - loss of precision converting from 32 to 16
    % bits
    regimage = 'FirstHybDAPI.tif'; % add as an optional argument
    [hybnum, tformsAll] = RegisterImages(compilePath, fieldOfView, regimage, channelsall, prehybnum, regis, seqDataFolder);
    

    %% Get the Threshold
    debug = 0;
    threshold = [];
    % Get Threshold of each set of Images if there are no threshold values
    if isempty(threshold)
        threshold = findthreshold_sequential(folderArray, hybnum, fieldOfView);
        % add saveDateString later in beginning of save file
        processDataPath = [savePath filesep dateSaveString '_Pos' num2str(fieldOfView) 'processImageData_YodaiIntronRep4'];
        save(processDataPath, 'backcorrections', 'chTforms', 'threshold');
    end


    %% Find Dots
    roiCytoPath = [organizedImageDir filesep 'RoiSetCyto']; % path to roiset;maybe test something elle with the coded
    roiNuclPath = [organizedImageDir filesep 'RoiSet'];
    channelArray = 1:numOfChannels; % number of channels to find dots
    debug = 0;
    [color, regMax, dots, Data] = FindDotsSeq(roiCytoPath, length(folderArray), channelArray, threshold,'roi', debug, hybnum);

    % save the dots
    save([savePath filesep dateSaveString '_Pos' num2str(fieldOfView) '-hyb-output_YodaiHybnumRep4.mat'],'Data','regMax','dots','color');


    %% Create Final Data Set Using ROIs 
    outputCyto = seqcopycompiler(fieldOfView, Data, roiCytoPath, roiNuclPath, length(folderArray), numOfChannels);
    save([savePath filesep dateSaveString '_Pos' num2str(fieldOfView) 'SequentialHybOutputCyto'], 'outputCyto');

    
    %% Finishing Diary Output
    % Save time completed
    dateEnd = datetime;
    dateStartToEnd = [dateStart dateEnd]; 
    dateTotal = diff(dateStartToEnd);
    [h,min,s] = hms(dateTotal);
    dateTotalString = sprintf('Total Elapsed Time: %.0f hours %.0f minutes and %.0f seconds...\n\n', h, min, s);
    % Save the Data and Time the script started
    fprintf(['Completed ' functionName ' Script Successfully\n']);
    fprintf(['Date: ' datestr(dateEnd) '\n']);
    fprintf(dateTotalString);
    fprintf('**************************************************************************\n\n'); 

    % Turn off the diary 
    diary off
end