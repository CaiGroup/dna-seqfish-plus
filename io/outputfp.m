function first = outputfp(savePath, barcodekey, seeds,...
    numpointconsensus, numdotlocations, numpointspercell, numfinalpoints, ...
    position, error, sqrt, ch, iter, offpercent, first)
% outputs false positive data
%
%
% Date: 9/10/2019

    
    
    %% Declare variables
    % Variables
    listSavePath = fullfile(savePath, ['falsePositiveData-per-cell-Pos' num2str(position) '.csv']);
    fileID = fopen(listSavePath,'a+');
    if first
        fprintf(fileID,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
                        'fov', 'ch', 'error', 'sqrtradius', 'iter', 'off percent','true barcodes', ...
                        'false positives', 'consensus points', 'colocalized points', ...
                        'avg points', 'error rate', 'calling rate');%, ...
    end
    
    % find number of seeds for unfiltered points found
    numCells = size(seeds, 1);
    
    % Calculate total
    numpointconsC = zeros(1,numCells);
    numdotlocC = numpointconsC;
    truebarC = numpointconsC;
    falsepC = numpointconsC;
    for c = 1:numCells
        numpointconsC(c) = sum(cell2mat(numpointconsensus{c}), 1);
        numdotlocC(c) = sum(cell2mat(numdotlocations{c}), 1);
        %numfinal = mean(sum(cell2mat(numfinalpoints), 1));
        offidx = find(contains(barcodekey{c}.names,'offtarget'));
        onidx = find(~contains(barcodekey{c}.names,'offtarget'));
        while onidx(end) > size(numfinalpoints{c}, 1) % checks indices, but why are they different
            onidx(end) = [];
        end
        while offidx(end) > size(numfinalpoints{c}, 1)
            offidx(end) = [];
        end
        truebarC(c) = sum(cell2mat(numfinalpoints{c}(onidx, :)),1);
        falsepC(c) = sum(cell2mat(numfinalpoints{c}(offidx, :)), 1);
    end
    numpointspc = mean(cell2mat(numpointspercell{1}));
    numpointcons = mean(numpointconsC);
    numdotloc = mean(numdotlocC);
    truebar = mean(truebarC);
    falsep = mean(falsepC);
    
    
    % print each line
    fprintf(fileID,'%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.3f,%.3f\n', ...
        position, ch, error, sqrt, iter, offpercent, truebar, falsep, numpointcons, numdotloc, numpointspc, falsep/truebar, truebar/numpointcons);%, numpointspc);

    fclose(fileID);
    first = false;

end
