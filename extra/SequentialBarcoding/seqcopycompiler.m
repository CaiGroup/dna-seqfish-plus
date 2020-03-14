function output = seqcopycompiler(posnum, Data, roiMainPath, roiFindPath, hybNumber, channelNumber)
% seqcopycompiler produces an output file of the copy number of each roi in
% the cytoplasm, matches the roi's of the nucleus, caclulates area and
% centroids of the cytoplasm roi; and finally returns the output file as a
% cell array.
% Order of the output file is as follows below
% 1:posnum, 2: # of ROI Nucleus (DAPI) 3:# of ROI Cyto, 4: x centroid of ROI, 5: y centroid of ROI,
% 6:area, 7-: copy, positions, intensity; copy, intensity; (hyb1 channel1, 2,3, hyb2 channel 1,2,3...) 
%
% Author: Yodai Takei
% Email: ytakei@caltech.edu
% Modified By: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/23/2019

    %% Declare the variables
    tot_pixel = 2048;
    numberOfFirstColumns = 6;
    vertex = selfseg(roiMainPath);
    cols = numberOfFirstColumns + channelNumber * hybNumber;
    rows = length(vertex);
    output = cell(rows, cols);

    %% Start the function 
    for roi = 1:length(vertex) % ROI of cytoplasm
        BWcyto = poly2mask(vertex(roi).x, vertex(roi).y, tot_pixel,tot_pixel);
        area_cyto = sum(sum(BWcyto));
        cyto = regionprops(BWcyto,'centroid');
        output{roi,1} = posnum;

        % find the number of roi.....use a function here to find it
        roiNucl = findroi(roiFindPath, vertex(roi));
        output{roi,2} = roiNucl;

        output{roi,3} = roi;
        output{roi,4} = cyto(1).Centroid(1);
        output{roi,5} = cyto(1).Centroid(2);
        output{roi,6} = area_cyto;
    end


    for hyb = 1:hybNumber
        copymRNA = Data(hyb).copy;
        for roi = 1:length(vertex)
            for channel = 1:channelNumber
                copyValues = num2cell(copymRNA(channel, :)');
                colIndex = numberOfFirstColumns + (hyb-1)*channelNumber + channel;
                output(:, colIndex) = copyValues;
            end
        end
        clearvars copymRNA
    end
end