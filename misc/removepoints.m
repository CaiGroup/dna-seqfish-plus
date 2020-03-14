function points = removepoints(points, dotlocations, removeInd)
    % removes points using dotlocations, and indices if provided
    % 
    % Date: 8/29/2019

    %% Take away the points here if getting offtarget barcodes
    disp('Take away points from on-target barcodes');

    numRemovedPoints = 0;
    numCells = size(dotlocations,1);
    for cellIndex = 1:numCells
        fprintf('Removing Points in cell %.0f of %.0f...\n', cellIndex, numCells);
        numGenes = size(dotlocations{cellIndex},1);
        if isempty(removeInd)
            % remove all genes
            geneArray = 1:numGenes;
        else
            % only remove genes with specified indices
            geneArray = removeInd;
        end
        for geneIndex = geneArray
            numPoints = size(dotlocations{cellIndex}{geneIndex, 1},1);
            for pointIndex = 1:numPoints
                % grab all the variables for the point
                x = dotlocations{cellIndex}{geneIndex,1}(pointIndex,1);
                y = dotlocations{cellIndex}{geneIndex,1}(pointIndex,2);
                z = dotlocations{cellIndex}{geneIndex,1}(pointIndex,3);
                serialhyb = dotlocations{cellIndex}{geneIndex,2}(pointIndex);
                barcoderound = dotlocations{cellIndex}{geneIndex,3}(pointIndex);

                % find the index of the row
                [row, col] = find(points{barcoderound}(serialhyb).channels(:,1) == x & points{barcoderound}(serialhyb).channels(:,2) == y & points{barcoderound}(serialhyb).channels(:,3) == z);

                % delete the point
                if ~isempty(row)
                    points{barcoderound}(serialhyb).channels(row,:) = [];
                    points{barcoderound}(serialhyb).intensity(row,:) = [];
                    points{barcoderound}(serialhyb).scaledIntensity(row,:) = [];
                    numRemovedPoints = numRemovedPoints + 1;
                end
            end
        end
    end

end
                