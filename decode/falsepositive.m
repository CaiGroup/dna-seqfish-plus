function [finalPosListFP, dotlocationsFP, numtotalConsensusSpotsFP, seedsFP, barcodekeyfull] ...
    = falsepositive(finalPosList, numtotal, dotlocations, points, avgpointsperround, ...
    barcodekeyfull, savePath, projectName, numRounds, numChannels, varargin)
% function calculates the false positive rate by removing found points,
% then using the whole codebook to find barcodes.
%
% Current Method for calculating false positives:
% - decode barcodes with current barcode
% - remove points used in dotlocations
% - use the full barcode key to calculate false positive rate
% 
%
% To do:
% 1. Make function simpler - less input variables
% 2. Add different types of ways to calculate false positives 
%
% Date: 8/1/2019
% Author: Nico Pierson
%% Set up optional Parameters

    argsLimit = 2;
    numvarargs = length(varargin);
    if numvarargs > argsLimit
        error('src:decodeimages:TooManyInputs', ...
            'requires at most 2 optional inputs');
    end
    % Error for type of arguments
    if numvarargs > 1
        if ~ischar(varargin{1}) && ~isscalar(varargin{1})
            error('src:decodeimages:WrongInput', ...
                'decodeimages var typedots requires type int');
        end
    end
    if numvarargs >= argsLimit 
        if ~isnumeric(varargin{2}) && ~isscalar(varargin{2})
            error('src:decodeimages:WrongInput', ...
                'decodeimages superresolve requires type int');
        end
    end
    % set defaults for optional inputs
    optargs = {1, 6};
    % assign defaults
    optargs(1:numvarargs) = varargin;
    % Default Value of ref image is 1
    [alloweddiff, sqrtradius] = optargs{:};
    
    
    %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    dateSaveString = datestr(dateStart, formatDate);
    
    
    %% Variables
    minNumSeeds = numRounds - alloweddiff;
    % Directories for saving the csv files for each position
    roiPathName = [projectDir filesep 'RoiSet.zip'];
    saveFileNameEnd = [projectName '-Pos' num2str(position) '-' dateSaveString];
    
    
    
    %% False Positive Check
    % Remove the Points
    [pointsFP, numRemovedPoints] = removepoints(dotlocations, points);
    if exist(roiPathName, 'file')
        pointsFP = filterpointspercell(roiPathName, pointsFP);
    end

    % Find barcodes
    [dotlocationsFP,seedsFP] = BarcodeFinder(numChannels, pointsFP, numRounds, barcodekeyfull,sqrtradius,alloweddiff);

    % Get the Final Positiion
    [numtotaldotconsensusFP, dotlocationsFP, numtotaldotlocationsFP] = filterseedsv2(seedsFP, dotlocationsFP, minNumSeeds);
    [finalPosListFP, numCells] = assignspots2cells(roiPathName, dotlocationsFP(:,6));

    % Check the number of false positives per barcode
    offtargetbarcodegraph(finalPosList, finalPosListFP, savePath, saveFileNameEnd);
    
    
    
    %% Print Data for false positive rate
    % Calculate efficiency of dotlocations and the size of the barcode
    % space
    offtargetpercell = numtotaldotconsensusFP / numCells;
    numofftargetefficiencyrate = numtotaldotlocationsFP / avgpointsperround;
    numofftargetbarcode = size(barcodekeyfull, 1);
    fprintf('\nFalse Positive Data ---------\n');
    fprintf('Efficiency of on-target barcodes: %.0f\n', numofftargetefficiencyrate); 
    fprintf('Size of off-target barcodekey: %.0f\n', numofftargetbarcode); 
    fprintf('Number of Cells: %.0f\n', numCells);
    fprintf('Average number of off-target barcodes per cell: %.0f\n', offtargetpercell); 
    fprintf('Number of points removed: %.0f\n', numRemovedPoints);
    fprintf('Number of consnsensus points found: %.0f\n', numtotaldotconsensusFP);
    
    
    %% Save the Data
    saveFileName = ['dataFP-' saveFileNameEnd];
    saveDataPath = fullfile(savePath, saveFileName);
    save(saveDataPath, 'finalPosListFP', 'dotlocationsFP', 'numtotaldotconsensusFP', ...
        'seedsFP', 'pointsFP', 'fullbarcodekey', 'numRemovedPoints', ...
        'numtotaldotlocationsFP', 'numofftargetbarcode', 'offtargetpercell', ...
        'numofftargetefficiencyrate')

    
end
