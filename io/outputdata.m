function [] = outputdata(savePath, projectSaveName, barcodekey, finalPosListC, posListC, dotlocations, ...
    numpointconsensusC, numdotlocationsC, numpointspercellC, numfinalpointsC, seedsC, option)
% outputs the data into csv gene by cell matrices or into csv files for
% each position
%
%
% Date: 8/8/2019



    %% Check options
    if ~ischar(option) 
        error('src:outputdata:WrongInput', ...
            'outputdata var option requires type string');
    elseif ~strcmp(option, 'matrix') && ~strcmp(option, 'list')
        error('myfun:outputdata:WrongInput', ...
            'outputdata var option requires type string: "matrix" or "list"');
    end
    
    
    
    %% Declare variables
    numGenes = size(barcodekey.names,1);
    numCells = size(finalPosListC,1);
    saveEnd = projectSaveName;
    
    
    %% output the data
    switch option
        case 'matrix'
            
            finalPosList = cell(numGenes, numCells + 1);
            finalPosList(:,1) = barcodekey.names;
            posList = finalPosList;
            numpointconsensus = finalPosList;
            numdotlocations = finalPosList;
            seeds = finalPosList;
            numfinalpoints = finalPosList;
            numpointspercell = cell(1, numCells + 1);
            numpointspercell{1} = 'total points per cell';

            
            for c = 1:numCells % only for number of cells
                startIdx = 1;
                endIdx = size(finalPosListC{c}, 1);
                finalPosList(startIdx:endIdx, c+1) = finalPosListC{c};
                posList(startIdx:endIdx, c+1) = posListC{c};
                numpointconsensus(startIdx:endIdx, c+1) = numpointconsensusC{c};
                numdotlocations(startIdx:endIdx, c+1) = numdotlocationsC{c};
                seeds(startIdx:endIdx, c+1) = seedsC{c};
                numpointspercell{1, c+1} = numpointspercellC{c};
                numfinalpoints(startIdx:endIdx, c+1) = numfinalpointsC{c};
            end
            
            % make a function to print all the variables:
            % printall(...,...,saveBegin, saveEnd);
            printcsv(finalPosList, fullfile(savePath, ['finalPosList-' saveEnd]));
            printcsv(posList, fullfile(savePath, ['posList-' saveEnd]));
            printcsv(numpointconsensus, fullfile(savePath, ['numpointconsensus-' saveEnd]));
            printcsv(numdotlocations, fullfile(savePath, ['numdotlocations-' saveEnd]));
            printcsv(seeds, fullfile(savePath, ['seeds-' saveEnd]));
            printcsv(numpointspercell, fullfile(savePath, ['numpointspercell-' saveEnd]));
            printcsv(numfinalpoints, fullfile(savePath, ['numfinalpoints-' saveEnd]));

        case 'list'
            % Variblles
            listSavePath = fullfile(savePath, ['outputData-' saveEnd '.csv']);
            fileID = fopen(listSavePath,'w');
            fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                            'cellID', 'geneID', 'regionID', 'x', 'y', 'z', ...
                            'seeds', 'intensity');
            fprintf(fileID,'\n'); % first new line doesn't create a new line
            
            
            for c = 1:numCells
                % go through the dotlocations 1 = poslist, 2 = round, 3 =
                % channel, 4 = intnesity, 5 = scaled intensity, 6 = finalposlist, 7 =
                % numdotlocations, 8 = numpointconsensus, 9 = barcode
                numRows = size(dotlocations{c}, 1);
                for r = 1:numRows
                    if ~isempty(numfinalpointsC{c}{r}) && numfinalpointsC{c}{r} ~= 0
                        numPointsPerRow = numfinalpointsC{c}{r};
                        for p = 1:numPointsPerRow
                            cellID = c;
                            geneID = barcodekey.names{r};
                            regionID = p; % regionID is the same gene in the same cell
                            x = dotlocations{c}{r,6}(p,1);
                            y = dotlocations{c}{r,6}(p,2);
                            z = dotlocations{c}{r,6}(p,3);
                            seed = dotlocations{c}{r,10}(p); 
                            intensity = dotlocations{c}{r,4}(p);

                            % print each line
                            fprintf(fileID,'%.0f,%s,%d,%.3f,%.3f,%.3f,%d,%.0f\n', ...
                                cellID, geneID, regionID, x, y, z, seed, intensity);
                        end
                    end
                end
            end
            fclose(fileID);
            
            
            %% make matrices for saving
            finalPosList = cell(numGenes, numCells + 1);
            finalPosList(:,1) = barcodekey.names;
            posList = finalPosList;
            numpointconsensus = finalPosList;
            numdotlocations = finalPosList;
            seeds = finalPosList;
            numfinalpoints = finalPosList;
            numpointspercell = cell(1, numCells + 1);
            numpointspercell{1} = 'total points per cell';

            
            for c = 1:numCells % only for number of cells
                startIdx = 1;
                endIdx = size(finalPosListC{c}, 1);
                finalPosList(startIdx:endIdx, c+1) = finalPosListC{c};
                posList(startIdx:endIdx, c+1) = posListC{c};
                numpointconsensus(startIdx:endIdx, c+1) = numpointconsensusC{c};
                numdotlocations(startIdx:endIdx, c+1) = numdotlocationsC{c};
                seeds(startIdx:endIdx, c+1) = seedsC{c};
                numfinalpoints(startIdx:endIdx, c+1) = numfinalpointsC{c};
            end
            numpointspercell(1, 2:end) = numpointspercellC';
            % make a function to print all the variables:
            % printall(...,...,saveBegin, saveEnd);
            
            %printcsv(finalPosList, fullfile(savePath, ['finalPosList-' saveEnd]));
            %printcsv(posList, fullfile(savePath, ['posList-' saveEnd]));
            %printcsv(numpointconsensus, fullfile(savePath, ['numpointconsensus-' saveEnd]));
            %printcsv(numdotlocations, fullfile(savePath, ['numdotlocations-' saveEnd]));
            %printcsv(seeds, fullfile(savePath, ['seeds-' saveEnd]));
            %printcsv(numpointspercell, fullfile(savePath, ['numpointspercell-' saveEnd]));
            %printcsv(numfinalpoints, fullfile(savePath, ['numfinalpoints-' saveEnd]));

            
        otherwise
            error 'myfun:outputdata:WrongInput outputdata var option is not valid';
    end
    
    %% save the data
    dataSaveName = fullfile(savePath, ['data-' saveEnd '.mat']);
    save(dataSaveName, 'finalPosList', 'posList', 'dotlocations', ...
    'numpointconsensus',  'numdotlocations', 'numpointspercell', 'seeds', ...
    'barcodekey', 'dotlocations', 'numfinalpoints');


end


%{
        %% Initialize Date for saving files
    dateStart = datetime;
    formatDate = 'yyyy-mm-dd';
    dateSaveString = datestr(dateStart, formatDate);

    %% initialize the variables
    % Directories and Paths
    expType = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius)];
    minNumSeeds = numRounds - alloweddiff;
    % Directories for saving the csv files for each position
    analysisFolder = 'analysis';
    roiPathName = [projectDir filesep 'RoiSet.zip'];
    savePath = fullfile(projectDir, analysisFolder, expType);
    if ~exist(savePath, 'dir')
        mkdir(savePath);
    end
    saveFileNameEnd = [projectName '-Pos' num2str(position) '-' dateSaveString];
    finalPosListDir = [savePath filesep 'FinalPosList' saveFileNameEnd]; % Main data to sav


    %% Read the Barcode excel file
    if isempty(barcodekey)
        barcodeFolder = 'barcodekey';
        barcodekeyPath = getfile([projectDir filesep barcodeFolder], 'barcodekey');
        barcodekey = readbarcode(barcodekeyPath, 'header');
    end
%}