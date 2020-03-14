function [] = outputdatahpc(saveDir, experimentDir, experimentName, pos, ...
                        ch, err, sqrtrad, iter, offpercent, numCells, first)
% organizes the saved files from the output and turns them into cells that
% are converted into an output list and matrix
%
% Still needs work: need to output the false positive rate based on the
% on-target and off-target hits
%
% Date: 9/9/2019

    %% load the data for each cell in the 
    % initialize variables
    posList = cell(numCells, 1);
    finalPosList = posList;
    dotlocations = posList;
    numpointconsensus = posList;
    numdotlocations = posList;
    seeds = posList;
    numfinalpoints = posList;
    numpointspercell = posList;
    barcodekey = posList;
    % Combine data together
    for i = 1:numCells
        explabel = ['minSeeds' num2str(i) 'Pos' num2str(pos) '-Cell' ...
            num2str(i) '-' num2str(err) 'error-sqrt' num2str(sqrtrad) ...
            '-iter' num2str(iter) '-ch' num2str(ch) '-.' num2str(offpercent) 'offpercent'];
        dataPath = fullfile(saveDir, ['decodeData-' explabel '-' experimentName '.mat']);
        if exist(dataPath, 'file') == 2
            d = load(dataPath, 'finalPosList', 'posList', 'dotlocations', 'numpointconsensus', ...
                'numdotlocations', 'numfinalpoints', 'seeds', 'numpointspercell', 'barcodekey', 'numpointspercell');
            finalPosList{i} = d.finalPosList;
            posList{i} = d.posList;
            dotlocations{i} = d.dotlocations;
            numpointconsensus{i} = d.numpointconsensus;
            numdotlocations{i} = d.numdotlocations;
            seeds{i} = d.seeds;
            numfinalpoints{i} = d.numfinalpoints;
            numpointspercell{i} = d.numpointspercell;
            barcodekey{i} = d.barcodekey;
        end
    end
    
    
    

    %% Organize the cell arrays to output data
    %explabel = [num2str(alloweddiff) 'error-sqrt' num2str(sqrtradius) '-iter' num2str(iter) '-ch' num2str(channel)];
    projectSaveName = ['Pos' num2str(pos) '-' num2str(err) 'error-sqrt' num2str(sqrtrad) ...
            '-iter' num2str(iter) '-ch' num2str(ch) '-' experimentName];
    savePath = fullfile(experimentDir, 'output-data');
    if exist(savePath,'dir') ~= 7
        mkdir(savePath);
    end
    saveoption = 'list';
    %outputdata(savePath, projectSaveName, d.barcodekey, finalPosList, posList, dotlocations, ...
    %    numpointconsensus,  numdotlocations, numpointspercell, numfinalpoints, seeds, saveoption);
    
    
    
    %% Remove empty cells because some files will not be present
    numpointconsensus = numpointconsensus(~cellfun('isempty',numpointconsensus));
    numdotlocations = numdotlocations(~cellfun('isempty',numdotlocations));
    seeds = seeds(~cellfun('isempty',seeds));
    numfinalpoints = numfinalpoints(~cellfun('isempty',numfinalpoints));
    numpointspercell = numpointspercell(~cellfun('isempty',numpointspercell));
    barcodekey = barcodekey(~cellfun('isempty',barcodekey));
    
    
    %% calculate the false positive rate
    % details: use the dataForFP to calculate the false positive rate
    outputfp(savePath, barcodekey, seeds, ...
        numpointconsensus, numdotlocations, numpointspercell, numfinalpoints, ...
        pos, err, sqrtrad, ch, iter, offpercent, first)
            % What is the rate for this
        % the rate is determined by the number of false positives vs tru on
        % target barcodes...or per cell rate as well..on-target efficiency
        %
        % number of on-target / average points per round
        % off-target / average points per round
        % put the number of each

end