function stringOutput = mat2strFormat(matrix)
% format string output to mathematica format
%
% example:
% tm = [203.23, 894.23, 12.3; 432.4, 930.23, 4.3];
% ouputs to:
% tm = {{203.23, 894.23, 12.3}, {432.4, 930.23, 4.3}}

    numberElements = size(matrix, 1);
    numberCols = size(matrix, 2);
    stringOutput = '{';
    
    for r = 1:numberElements
        stringOutput = strcat(stringOutput, '{');
        for c = 1:numberCols
            numString = num2str(matrix(r,c));
            stringOutput = strcat(stringOutput, numString);
            if c ~= numberCols
                stringSepCol = ',';
            else
                stringSepCol = '}';
            end
            stringOutput = strcat(stringOutput, stringSepCol);
        end
        if  r ~= numberElements
            stringSepRow = ',';
        else
            stringSepRow = '}';
        end
        stringOutput = strcat(stringOutput, stringSepRow);
    end

end