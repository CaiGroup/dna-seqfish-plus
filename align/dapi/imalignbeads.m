function [hybImsBeadAligned] = imalignbeads(experimentName, experimentDir, ...
    position, channelArray, hybImsDapiAligned, rawHybIms, tformDapi, varargin)
% Take the images from the dapi alignment, and align the images with
% fiduciary beads.
%
% Inputs: rawImages with dapi, hyb aligned images with dapi, tform for dapi
% (hybCycel by 1 cell array)
%
% Outputs: ???
%
% Paths:
% 1. addpath('C:\github\streamline-seqFISH\src\AlignImages\bfmatlab',
% '-end'); add the check function for this, and update the bfreader like
% the other packages. 
%
% Dependencies:
% 1. FindThreshold package:
% addpath('C:\github\streamline-seqFISH\src\FindThreshold', '-end');
% 2. preprocessing package: 
% addpath('C:\github\streamline-seqFISH\src\preprocessing', '-end');
%
% To Do:
% 1. Add packages to path function or check if in path
%
% Update:
% 1. 9/11/2019: sped up autothreshold for beads; threshold beads first; and
% added option to provide dapiRefPath to align to different dapi images
%
% Assumptions: threshold for hybs and fixed beads are the same
%
% Assumptions: hybIms are assumed to be aligned by tform and background
% subtracted
% 
% Check the transformations for the dots and compare the difference between
% the dapi alignment and the bead alignment


    %% Set up optional Parameters
    argsLimit = 6;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:imalignbeads:TooManyInputs', ...
            'requires at most 5 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 0
        if ~ischar(varargin{1}) 
            error('src:imalignbeads:WrongInput', ...
                'imalignbeads beadfoldername requires type string');
        end
    end
    if numvarargs > 1
        if ~ischar(varargin{2}) 
            error('src:imalignbeads:WrongInput', ...
                'imalignbeads dapiRefPath requires type string');
        end
    end
    if numvarargs > 2
        if ~isnumeric(varargin{3}) && ~isscalar(varargin{3})
            error('src:imalignbeads:WrongInput', ...
                'imalignbeads var searchradius requires type int');
        end
    end
    if numvarargs > 3
        if varargin{4} ~= 1 && varargin{4} ~= 0
            error('src:imalignbeads:WrongInput', ...
                'imalignbeads process requires type boolean');
        end
    end
    if numvarargs > 4
        if ~isnumeric(varargin{5}) && ~isscalar(varargin{5})
            error('src:imalignbeads:WrongInput', ...
                'imalignbeads backradius requires type int');
        end
    end
    if numvarargs >= argsLimit
        if ~ischar(varargin{6}) 
            error('src:imalignbeads:WrongInput', ...
                'imalignbeads var typedots requires type string');
        elseif ~strcmp(varargin{6}, 'exons') && ~strcmp(varargin{6}, 'introns') && ~strcmp(varargin{6}, 'exons2d')
            error('src:imalignbeads:WrongInput', ...
                'imalignbeads var typedots requires type string: "exons" or "introns" or "exons2d"');
        end
    end
    
    % set defaults for optional inputs
    optargs = {'initial_fiducial_markers', [], sqrt(4), false, 3, 'exons'};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [beadfoldername, dapiRefPath, searchradius, process, backradius, typedots] = optargs{:};
    


    %% add path to use Find Threshold package
    currentPath = pwd;
    addpath(fullfile(currentPath, '..', 'FindThreshold'), '-end');
    
    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    dateSaveString = datestr(dateStart, formatDate);

    

    %% Initialize Variables
    targetNumPoints = 400; % number of target points to threshold for beads...
    sliding = false;
    endingString = [experimentName '-Pos' num2str(position) '-' dateSaveString];
    saveDir = fullfile(experimentDir, 'beaderror', ['pos' num2str(position)]);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    processBeads = true;
    debug = 1;
    typedotsFixed = 'introns';

    
    
    %% threshold for beads in each channel first
    imBeadsPath = fullfile(experimentDir, beadfoldername, ['MMStack_Pos' num2str(position) '.ome.tif']);
    [beadImsTemp, numBeadDapi, numBeadZSlice, ~, ~] = grabimseries(imBeadsPath, position);
    numHyb = numBeadDapi -1;
    fixedThreshold = []; % Can add threshold if already set
    [fixedPoints, fixedThreshold] = matchpointsims(beadImsTemp(1:numHyb), ...
        fixedThreshold, typedotsFixed, processBeads, backradius, sliding, ...
        [], searchradius, 0);
    
    
    
    %% Grab Dapi Images
    numHybCycles = size(hybImsDapiAligned, 1);
    %numHybCycles = 80;
    %position = 0;
    % Get first dapi Image
    if isempty(dapiRefPath)
        dapiImsPath = fullfile(experimentDir, 'HybCycle_0', ['MMStack_Pos' num2str(position) '.ome.tif']);
    else
        dapiImsPath = fullfile(dapiRefPath, ['MMStack_Pos' num2str(position) '.ome.tif']);
    end
    [dapiIms, numDapi, numZSlice, ~, ~] = grabimseries(dapiImsPath, 0);
    %imBeadsPath = 'D:\MATLAB\CaiLab\Development\YodaiPointsAlignment\initial_fiducial_markers\MMStack_Pos0.ome.tif';
    
    
    
    
    %% Align Images to the first dapi images in HybCycle_0
    tformDapi4Beads = grabtform(beadImsTemp{numBeadDapi}, dapiIms{numDapi});
    beadIms = applydapitform(beadImsTemp, tformDapi4Beads);


    
    %% Get the fiduciary points across all channels
    % addpath('C:\Users\nicog\Documents\fiji-win64\Fiji.app\scripts', '-end');
    % need to choose the threshold, which will be used for all other hybs
    fprintf('Aligning Beads to Hyb0 Reference Image\n');
    
    % use test threshold - normally have to manually threshold
    %fixedThreshold = [500,450,500];
    %fixedThreshold = [];
    saveFigPathBeads = fullfile(saveDir, ['FixedPoints-Pos' num2str(position) '.fig']);
    [fixedPoints, fixedThreshold] = matchpointsims(beadIms(1:numHyb), ...
        fixedThreshold, typedotsFixed, processBeads, backradius, sliding, [], searchradius, debug, saveFigPathBeads);
    
    
    %% Get Points only in ROI - add as option later
    %{
    roiPath = fullfile(experimentDir, 'segmentation', ['Pos' num2str(position)], 'RoiSet.zip');
    fixedPointsRoi = getpointsinroi(fixedPoints(1).channels, roiPath);
    numRefPoints = length(fixedPointsRoi);
    %}
    fixedPointsRoi = fixedPoints;
    
    %% Declare Variables for Global Alignment
    saveFilePath = fullfile(saveDir, ['beadError-Pos' num2str(position) '-' endingString '.csv']);
    hybImsBeadAligned = cell(numHybCycles, numHyb);
    tformBeads = cell(numHybCycles, 1);
    fixedMatchPoints = tformBeads;
    movingMatchPoints = tformBeads;
    beadMatchPoints = tformBeads;
    rawThreshold = []; % use the same rawThreshold for getting the reference points in the beaderror.m
    %targetNumPoints = numRefPoints * 15;
    
    
    %% Parfor loop to get all threshold
    numWorkers = 10; % number of cores to use
    numIters = 25; % number of iterations to determin autothreshold
    typedots2 = 'exons';
    hybThreshold = autothresholdnumpoints(hybImsDapiAligned, targetNumPoints, numIters, numWorkers, channelArray, typedots2);

    
    
    
    %% Find the global alignment between images and the fiduciary points
    for h = 1:numHybCycles
        

        % Get the moving points (gaussian) across all channels using the aligned dapi Images
        % most likely the hybIms are processed
        % might need (hybImsDapiAligned(h,1:numHyb)
        %typedotshyb = 'exons';
         % images for rep2 yodai dna fish experiment don't match channel 3
        debug = 1;
        saveFigPathHyb = fullfile(saveDir, ['HybPoints-Pos' num2str(position) 'HybCycle-' num2str(h) '.fig']);
        [movingPoints, hybThreshold{h}] = matchpointsims(hybImsDapiAligned(h,channelArray), ...
            hybThreshold{h}, typedots, process, backradius, sliding, ...
            targetNumPoints, searchradius, debug, saveFigPathHyb); % option to back subtract because images are background subtracted sometimes
        
        

        % Match fixedPoints to the movingPoints, and calculate the tform
        % for beads: figure out if .ref or .points for fixed and moving
        % points
        % - use channel 1
        % - use the 3 option, grabbing more points - ref refers to the 
        [fixedMatchPoints{h}, ~, movingMatchPoints{h}, ~] = matchpoints(fixedPointsRoi(1).channels, movingPoints(1).channels, backradius);
        tformBeads{h} = getglobaltform(movingMatchPoints{h}.points, movingMatchPoints{h}.ref); % check if the fixedMatch points is actually the ref points
       
        

        %% Apply tform to the images - getting points form image is different
        %numRefBeadPoints = [];
        beadMatchPoints{h} = transformPointsForward(tformBeads{h}, movingMatchPoints{h}.ref);
        fprintf('Pos %.0f, HybCycle %.0f:\nTform for Hybs:\n', position, h);
        hybImsBeadAligned(h,:) = applydapitform(hybImsDapiAligned(h,:), tformBeads{h});
        %[beadAlignedPoints, ~] = matchpointsims(hybImsBeadAligned(h,:), hybThreshold, typedots, process, backradius, sliding, numRefBeadPoints);
        %[~, ~, beadMatchPoints{h}, ~] = matchpoints(fixedPoints(1).channels, beadAlignedPoints{h}(1).channels, backradius);
        
        %% Error for alignment and figures
        % use channel 1 for fixed points
        rawThreshold = beaderror(saveFilePath, saveDir, h, rawHybIms(h,channelArray), hybImsBeadAligned(h,:), beadMatchPoints{h}, ...
            fixedPointsRoi(1).channels, movingMatchPoints{h}, tformBeads{h}, tformDapi{h}, ...
            position, typedots, backradius, sliding, rawThreshold, targetNumPoints);
        if h == 1
            fprintf('Threshold for raw images:,%.2f,%.2f,%.2f\n', rawThreshold);
        end
        
        
    end


    %% Save the images
    saveImName = ['beadAlignedImages-' endingString];
    saveImPath = fullfile(saveDir, saveImName);
    save(saveImPath, 'hybImsBeadAligned', 'tformBeads', 'tformDapi', ...
        'fixedMatchPoints', 'movingMatchPoints', 'beadMatchPoints', ...
        'hybThreshold', 'fixedThreshold', 'numHybCycles', ...
        'targetNumPoints', 'fixedPointsRoi', 'fixedPoints', '-v7.3');
end



%imHyb0Path = 'D:\MATLAB\CaiLab\Development\YodaiPointsAlignment\HybCycle_0\MMStack_Pos0.ome.tif';

% first hyb for initial dapi transformations
% use the hyb ims from the image processing to find dots and find the
% transformation
%[hybImsTemp, numDapi, numHybZSlice, ~, ~] = grabimseries(imHyb0Path, position);