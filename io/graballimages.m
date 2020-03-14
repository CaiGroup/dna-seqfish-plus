function I = graballimages(position, folderArray, experimentDir, numWorkers)

    % Initialize variables
    I = cell(length(folderArray), 1);
    % Initialize parfor pool
    p = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
    if poolsize <= 0
        parpool('local', numWorkers);
    end
    
    parfor folder = folderArray
        %get the images for each channel
        imageName = ['MMStack_Pos' num2str(position) '.ome.tif'];
        imagePath = fullfile(experimentDir, ['HybCycle_' num2str(folderArray(folder))], imageName);
        [I{folder+1}, ~, ~, ~, ~] = grabimseries(imagePath, position);
    end

end