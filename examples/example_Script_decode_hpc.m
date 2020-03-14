% script to output the hpc data for all the parameters used

% enter the directory for the data
%experimentDir = 'L:\FalsePositiveTest';
experimentDir = 'E:\Mike_2\Linus_10k_cleared_080918_NIH3T3';
experimentName = 'NIH3T3_080918';
saveDir = fullfile(experimentDir, 'hpc-data');
if exist(saveDir, 'dir') ~= 7
    mkdir(saveDir);
end
chs = 488;%[647 561 488];
positions = 0;%0:6; % only check position 0 
cells = {[1 5]};%{1:16, 1:14, 1:14, 1:18, 1:13, 1:15, 1:13};
errors = [1];%[0 1];%[0 1 2];
iters = 1:3;
sqrts = [2];%[1 2];%[1 2 3 4 6];
offpercents = [100];%[10, 25, 50, 75, 100];


for pos = positions
    first = true; % controls printing title in csv for false positive
    for ch = chs
        for error = errors
            for sqr = sqrts
                for iter = iters
                    for offpercent = offpercents
                        numCells = length(cells{pos+1});

                        % save the outputData, the point files, and the
                        % number of real barcodes vs off-target barcodes
                        % for each iteration....

                        outputdatahpc(saveDir, experimentDir, experimentName, pos, ...
                            ch, error, sqr, iter, offpercent, numCells, first);
                        first = false;
                    end
                end
            end
        end
    end
end