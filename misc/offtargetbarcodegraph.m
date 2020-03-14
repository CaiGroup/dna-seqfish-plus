function [] = offtargetbarcodegraph(finalPosList, finalPosListFP, savePath, saveEnding)
% script to turn offtarget copynumfinalsum of barcodes to an array with the
% barcode name as the first column and the numbeer of barcodes called per
% gene as the second column

%load('H:\CZI\organizehybs\Pos1\imagedata\postionSeedsList_2019-03-04.mat')

% sum the barcodes for each gene in each row for finalPosList
num = zeros(size(finalPosList, 1),1);
for i = 1:size(finalPosList, 1)
    for j = 2:size(finalPosList, 2)
        if ~isempty(finalPosList{i,j})
            num(i) = num(i) + size(finalPosList{i,j}, 1);
        end
    end
end
% sort
ontarget = sort(num, 'descend');
% list
ontarget_list(:,1) = 1:1:size(finalPosList, 1);
ontarget_list(:,2) = ontarget;

% load positionSeedsList or just copynumfinalrevised
% sum the barcodes for each gene in each row for finalPosList
num = zeros(size(finalPosListFP, 1),1);
for i = 1:size(finalPosListFP, 1)
    for j = 2:size(finalPosListFP, 2)
        if ~isempty(finalPosListFP{i,j})
            num(i) = num(i) + size(finalPosListFP{i,j}, 1);
        end
    end
end
% sort
offtarget= sort(num, 'descend'); % sort
offtarget_list(:,1) = (1+length(ontarget_list)):1:(length(ontarget_list)+size(finalPosListFP, 1)); % add length of on target
offtarget_list(:,2) = offtarget; % make list


% plot the histogram for true positives, then false positives
f1 = figure();
m1Data  = line(ontarget_list(:, 1), ontarget_list(:, 2));
set(m1Data                         , ...
      'LineStyle'       , '--'      , ...
      'Marker'          , '.'         );
    set(m1Data                         , ...
      'Marker'          , 'o'         , ...
      'MarkerSize'      , 3           , ...
      'MarkerEdgeColor' , 'none'      , ...
      'LineWidth'       , 2.0         , ...
      'Color'       , [.75 .75 1] , ...
      'MarkerFaceColor' , 'b' );
  hold on
m2Data  = line(offtarget_list(:, 1), offtarget_list(:, 2));
set(m2Data                         , ...
      'LineStyle'       , '--'      , ...
      'Marker'          , '.'         );
    set(m2Data                         , ...
      'Marker'          , 'o'         , ...
      'MarkerSize'      , 3           , ...
      'MarkerEdgeColor' , 'none'      , ...
      'LineWidth'       , 2.0         , ...
      'Color'       , 'r' , ...
      'MarkerFaceColor' , 'r' );
  hold on
title('On/Off-Target Barcodes False Positive Histogram')
xlabel('Barcodes')
ylabel('# of barcodes found')
legend('on-target', 'off-target','Location','best');
axis([-250 2500 -1000 1000])
% per barcode data...divide by number of ROIs for per cell
mean(offtarget_list(:,2))
std(offtarget_list(:,2))

mean(ontarget_list(:,2))
std(ontarget_list(:,2))
fileName1 = ['HistogramOnOff-TargetBarcodesPerBarcode-' saveEnding];
save1FigPath = fullfile(savePath, fileName1);
savefig(f1, save1FigPath);

%% histogram for cell with number of total detected barcodes
num = zeros(size(finalPosList, 2)-1, 1);
for i = 2:size(finalPosList, 2)
    for j = 1:size(finalPosList, 1)
        if ~isempty(finalPosList{j,i})
            num(i-1) = num(i-1) + size(finalPosList{j,i}, 1);
        end
    end
end
numTotalPerCell_ontarget(:,1) = 1:(size(finalPosList, 2)-1);
numTotalPerCell_ontarget(:,2) = num';
% offtarget - load offtarget finalPosList
num = zeros(size(finalPosListFP, 2)-1, 1);
for i = 2:size(finalPosListFP, 2)
    for j = 1:size(finalPosListFP, 1)
        if ~isempty(finalPosListFP{j,i})
            num(i-1) = num(i-1) + size(finalPosListFP{j,i}, 1);
        end
    end
end
numTotalPerCell_offtarget(:,1) = 1:(size(finalPosListFP, 2)-1);
numTotalPerCell_offtarget(:,2) = num';


f2 = figure();
m1Data  = line(numTotalPerCell_ontarget(:, 1), numTotalPerCell_ontarget(:, 2));
set(m1Data                         , ...
      'LineStyle'       , '--'      , ...
      'Marker'          , '.'         );
    set(m1Data                         , ...
      'Marker'          , 'o'         , ...
      'MarkerSize'      , 8           , ...
      'MarkerEdgeColor' , 'none'      , ...
      'LineWidth'       , 2.0         , ...
      'Color'       , [.75 .75 1] , ...
      'MarkerFaceColor' , [.75 .75 1] );
hold on
m2Data  = line(numTotalPerCell_offtarget(:, 1), numTotalPerCell_offtarget(:, 2));
set(m2Data                         , ...
      'LineStyle'       , '--'      , ...
      'Marker'          , '.'         );
    set(m2Data                         , ...
      'Marker'          , 'o'         , ...
      'MarkerSize'      , 4           , ...
      'MarkerEdgeColor' , 'none'      , ...
      'LineWidth'       , 2.0         , ...
      'Color'       , [1 .5 .5] , ...
      'MarkerFaceColor' , [1 .5 .5] );
  xlabel('cells')
ylabel('number of barcodes')
  legend('on target', 'off target')
  title('Total Number of Decoded On and Off-Target Barcodes per Cells')
  
fileName2 = ['TotalNumberOnOff-TargetBarcodesPerCell-' saveEnding];
save2FigPath = fullfile(savePath, fileName2);
savefig(f2, save2FigPath);
  
% make csv for number of final positions
numTotalList = finalPosList;
for i = 1:size(finalPosList, 1)
    for j = 2:size(finalPosList, 2)
        if ~isempty(finalPosList{i,j})
            numTotalList{i,j} = size(finalPosList{i,j}, 1);
        else
            numTotalList{i,j} = 0;
        end
    end
end

% Matrix Form: getting mean per cell
numTotalListMat = zeros(size(finalPosList(:,2:end)));
for i = 1:size(finalPosList, 1)
    for j = 2:size(finalPosList, 2)
        if ~isempty(finalPosList{i,j})
            numTotalListMat(i,j) = size(finalPosList{i,j}, 1);
        else
            numTotalListMat(i,j) = 0;
        end
    end
end

numPerCell = sum(numTotalListMat);
medianPerCell = median(numPerCell(2:end));
meanPerCell = mean(numPerCell(2:end));
totalPerCell = sum(numPerCell);

%% other things
%number of dots called in one round = 12* 15,000 vs. how many total
%barcodes 40k
% need to consider points in the ROI


%{
%% Make script to find number of points in one barcode round
    vertex = selfseg(PathForROI);
    allcalled = [];
    for i = 1:length(vertex)
        for j = 1:hyb

            allcalled = [];
            for k = 1:channum
                include = inpolygon(points(j).dots(k).channels(:,1),points(j).dots(k).channels(:,2),vertex(i).x,vertex(i).y);
                allcalled = [allcalled; include];
            end
        end
    end
    numPointsOneBarcode = length(allcalled);
%}
