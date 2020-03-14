function [] = printtform(tform, fileName)
% function takes the tforms and prints the tforms:
% 
% take the data
% each x, y, and z - output into a csv file using fopen
% do the same for the colocalization to compare the efficiency
%
% Date: 8/7/2019
% Author: Nico Pierson


    %% Variables
    numTforms = length(tform);
    numDim = length(tform{1}.T);
    
    
    %% Open the file reader
    fileID = fopen(fileName,'w');
    fprintf(fileID,'Transformations for each folder named HybCycle\n');
    
    if numDim == 3
        fprintf(fileID,'%s,%s,%s\n', 'HybCycle', 'x', 'y');
    elseif numDim == 4
        fprintf(fileID,'%s,%s,%s,%s\n', 'HybCycle', 'x', 'y', 'z');
    end
    
    %% Go through each tform
    for i = 1:numTforms
        
        % get x, y, and z displacement
        x = tform{i}.T(numDim,1);
        y = tform{i}.T(numDim,2);
        if numDim == 3
            fprintf(fileID,'%.0f,%.3f,%.3f\n', i, x, y);
        elseif numDim == 4
            z = tform{i}.T(numDim,3);
            fprintf(fileID,'%.0f,%.3f,%.3f,%.3f\n', i, x, y, z);
        else
            error 'Wrong number of dimensions in printtform';
        end
        
    end
    fclose(fileID);
    

end