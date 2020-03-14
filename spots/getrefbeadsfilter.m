function [chaTform, refPointsAligned, intensity, firstThreshold, numRefPoints] ...
    = getrefbeadsfilter(experimentDir, position, beadfoldername, roimask, varargin) 
% gets bead ref points from bead images and transform them to the hyb0 ref
% dapi image
%
% Path Dependencies
%addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end');
%addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\', '-end');
%addpath('C:\github\streamline-seqFISH\src\AlignImages\', '-end');
%addpath('C:\github\streamline-seqFISH\src\FindThreshold\', '-end');

    %% Set up optional Parameters
    argsLimit = 6;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:adjustthreshold:TooManyInputs', ...
            'requires at most 6 optional inputs');
    end   
    % set defaults for optional inputs
    optargs = {'gaussian', false, ''};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [superres, manualThreshold, saveEnding, firstThreshold, numRefPoints] = optargs{:};
    filtersigma = true; %added


    saveDir = fullfile(experimentDir, 'points');
    if exist(saveDir, 'dir') ~= 7
        mkdir(saveDir);
    end

    backradius = 3;
    sliding = false;
    searchradius = 3;
    processBeads = true;
    typedotsFixed = 'log';
    


    % initial beads
    %beadfoldername = 'initial_fiducial_markers';
    listing = dir([experimentDir '\' beadfoldername '\*MMStack_Pos' num2str(position) '.ome.tif']);
    imBeadsPath = fullfile(experimentDir, beadfoldername, listing(1).name);
    %imBeadsPath = fullfile(experimentDir, beadfoldername, ['MMStack_Pos' num2str(position) '.ome.tif']);
    [beadImsTemp, numBeadDapi, numBeadZSlice, ~, ~] = grabimseries(imBeadsPath, position);
    numCh = numBeadDapi -1;
    fixedThreshold = ones(numCh, 1) * 10000;
    numpointError = ones(numCh, 1);
    if ~manualThreshold
        fixedThreshold = [15000 15000 12000];
        %{
        for ch = 1:numCh
            logFish = max(beadImsTemp{ch},[],3);
            thresh = multithresh(logFish,2);
            if thresh(2) < 1000 || thresh(2) > 15000
                fixedThreshold(ch) = 10000;%fixedThreshold(ch-1);
            else
                fixedThreshold(ch) = thresh(2);
            end
        end
        %}
    elseif ~isempty(firstThreshold)
        for ch = 1:numCh
            [~, ~, fixedThreshold(ch), numpointError(ch)] = ...
                autothresholdbynumberpoints(beadImsTemp{ch}, numRefPoints(ch), firstThreshold(ch), ...
                typedotsFixed, roimask);
        end
    else
        fixedThreshold = [];
    end
    
    [refPoints, fixedThreshold, matches] = matchpointsims(beadImsTemp(1:numCh), ...
    fixedThreshold, typedotsFixed, processBeads, backradius, sliding, ...
    [], searchradius, 0, pwd, superres, filtersigma);
    
    maxZ = size(beadImsTemp{1},3);
    maxXY = 2048;
    min = 1;

    % make save directory
    savePath = fullfile(experimentDir, 'analysis', saveEnding, beadfoldername);
    if exist(savePath, 'dir') ~= 7
        mkdir(savePath);
    end
    
    %% Save the intensity vs sigma
    f1 = figure('Name', 'int vs sigma - ch1');
    scatter(matches{1}.sigma, matches{1}.intensity);
    xlabel('sigma')
    ylabel('intensity')
    filenamefig = fullfile(savePath, ['int-vs-sigma-pos' num2str(position) '-ch1.fig']);
    savefig(f1, filenamefig);
    f2 = figure('Name', 'int vs sigma - ch2');
    scatter(matches{2}.sigma, matches{2}.intmatch);
    xlabel('sigma')
    ylabel('intensity')
    filenamefig = fullfile(savePath, ['int-vs-sigma-pos' num2str(position) '-ch2.fig']);
    savefig(f2, filenamefig);
    f3 = figure('Name', 'int vs sigma - ch3');
    scatter(matches{3}.sigma, matches{3}.intmatch);
    xlabel('sigma')
    ylabel('intensity')
    filenamefig = fullfile(savePath, ['int-vs-sigma-pos' num2str(position) '-ch3.fig']);
    savefig(f3, filenamefig);
    
    
    %{
    % Apply the tform of dapi from final beads to the dapi of Hyb0
    % get Dapi from Hyb0
    dapiPath = fullfile(experimentDir, 'HybCycle_0', ['MMStack_Pos' num2str(position) '.ome.tif']);
    [allhyb0, numhyb0Dapi, numhyb0ZSlice, ~, ~] = grabimseries(dapiPath, position);
    dapiRef = allhyb0{numhyb0Dapi};
    % find tform from the bead dapi to the dapiRef
    beadDapi = beadImsTemp{numBeadDapi};
    bead2dapiTform = grabtform(beadDapi, dapiRef);
%}


    intensity = cell(1,numCh);
    for ch = 1:numCh
        intensity{ch} = refPoints(ch).intensity;
    end
    
    %% filter the points
    for ch = 1:numCh
        removeInd = [];
        if ~isempty(roimask)
            adjustedPoints = refPoints(ch).channels;
            numPoints = size(refPoints(ch).channels,1);
            for i = 1:numPoints
                y = round(adjustedPoints(i,2));
                x = round(adjustedPoints(i,1));
                z = round(adjustedPoints(i,3));
                if z > maxZ
                    z = maxZ;
                elseif z < min
                    z = min;
                end
                if ~roimask(y, x, z)
                    removeInd = cat(1, removeInd, i);
                end
            end
            refPoints(ch).channels(removeInd,:) = [];
            refPoints(ch).intensity(removeInd,:) = [];
            intensity{ch}(removeInd) = [];
        end
    end
    if isempty(numRefPoints)
        numRefPoints = zeros(numCh, 1);
        for ch = 1:numCh
            numRefPoints(ch) = size(refPoints(ch).intensity,1);
        end
    end
    if isempty(firstThreshold)
        firstThreshold  = fixedThreshold;
    end
    
    
    
    %% print reference bead points
    %listSavePath = fullfile(savePath, ['ref-points-pos' num2str(position) '-' beadfoldername '-' endingString '.csv']);
    listSavePath = fullfile(savePath, ['ref-points-pos' num2str(position) '-' beadfoldername '.csv']);
    fileID = fopen(listSavePath,'w');
    fprintf(fileID, '%s,%s,%s,%s,%s\n', 'ch', 'x', 'y', 'z', 'int');
    
    save(fullfile(savePath, ['ref-points-pos' num2str(position) '-' beadfoldername '.mat']), ...
        'refPoints', 'position');
    
    for ch = 1:numCh
        pointsSize = length(refPoints(ch).channels);
        % apply transformation to the points
        pointsTemp = refPoints(ch).channels;
        for i = 1:pointsSize
            x = pointsTemp(i,1);
            y = pointsTemp(i,2);
            z = pointsTemp(i,3);
            int = intensity{ch}(i);
            fprintf(fileID, '%.0f,%.3f,%.3f,%.3f,%.0f\n', ch, x, y, z, int);
        end
    end
    fclose(fileID);
    
    refPointsAligned = refPoints;
    
    
    %% save images with points
    for ch = 1:numCh
        saveFigPath = fullfile(savePath, ['beadCheck-' beadfoldername '-pos' num2str(position) '-ch' num2str(ch) '.fig']);
        [dots, ~, dotsLogical, ~] = detectdotsv2(beadImsTemp{ch}, fixedThreshold(ch), typedotsFixed);
        % filter the points for Jonathan
        printfigure(beadImsTemp{ch}, dotsLogical, refPoints(ch).channels, saveFigPath);
    end
    
    %{
    for ch = 1:numCh
        % apply transformation to the points
        refPointsAligned(ch).channels = transformPointsForward(bead2dapiTform, refPoints(ch).channels);
    end
%}
    
    %% get chromatic abberations
    chaTform = loadchabberationtform(1:numCh, refPointsAligned);
end