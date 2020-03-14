function threshold = thresholdbarcodeims(experimentDir, position, numPseudoChannels, numBarcodes, typeofdots)
% function allows manual thresholding for each channel in a barcoded 
% experiment
%
% Dependences: packages CompareDotsError and Fiji, mij.jar in
% Fiji.app/scripts folder
%
% Outputs: threshold in barcode by channel matrix
%
% Testing: Not tested on mac or linux OS; 'introns' has not been tested.
%
% Variables:
% experimentDir is the directory that contains HybCycle[1-9]
% position starts from 0
% numPseudoChannels is the number of pseudochannels in the experiment
% numBarcodes is the number of barcodes in the experiment
% typeofdots is the type of dotdetection: 'exons' (laplaction of gaussian 
% peaks) or 'introns' (peaks using masks)
%
% Assumptions: 
% Images are organized by HybCycle for each round of imaging
% 
% Date: 7/29/2019
% Author: Nico Pierson

    tic
    threshold = thresholdbypos(experimentDir, position, numBarcodes, numPseudoChannels, typeofdots);
    toc

end