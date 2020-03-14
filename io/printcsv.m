function [] = printcsv(finalPosList, fileName)
% printcsv rounds the results in each cell, turns the gene names
% into strins and sorts by the column and saves into the xls file name
%
% Inputs: 
%
% Outputs: 
%
% Update: 2/22/2019
% Fixed bug for switching column and row indices, whereby it was removig
% the first row of the results.
%
% Author: Nico Pierson
% Date: December 6, 2018

zeroFinalPosList = cellfun(@removezeros, finalPosList, 'UniformOutput', false);
roundFinalPosList = cellfun(@round2dec, zeroFinalPosList, 'UniformOutput', false); % round values

% Transpose: DON't Because xlswrite can't write so many columns, but csv
% files can handle it
%transFinalPosList = roundFinalPosList'; 

% Convert all matrices to string
strFinalPosList = cell(size(roundFinalPosList));
strFinalPosList(:,1) = roundFinalPosList(:,1);

% error when the gene name starts with a number
strFinalPosList(:, 2:end) = cellfun(@mat2str, roundFinalPosList(:, 2:end), 'UniformOutput', false);
% Sort Rows by the first column
%sortFinalPosList = sortrows(strFinalPosList, 1);

% Write to excel
%xlswrite([fileName '.xls'], strFinalPosList);
% Write to csv
cell2csv([fileName '.csv'], strFinalPosList);