function cellin = round2dec(cellin)

if ~ischar(cellin)

    [row col] = size(cellin);

    for i = 1:row
        for j = 1:col
            cellin(i,j) = round(cellin(i,j)*100)/100;
        end
    end
end

