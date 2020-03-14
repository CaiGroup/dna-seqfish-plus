function num = add2num(x,y)

	fprintf('first values x: %s, y: %s.\n', x,y);
	num = str2num(x) + str2num(y);
	curr = 'C:\github\streamline-seqFISH\src\';
	filePath = fullfile(curr, 'testnumbers.txt');
	fileID = fopen(filePath, 'w');
	fprintf(fileID, '%s,%s,%s\n', 'x','y','result');
	fprintf(fileID, '%d,%d,%d\n', str2num(x),str2num(y),num);
	fclose(fileID);
end
