function [dotlocations,seeds] = BarcodeFinder(channels, points, hyb, barcodekey,radius,alloweddiff)
%#codegen
% decodes the barcoded images with the number of channels, the points, the 
% nubmber of rounds, barcodekey, search radius and number of barcode rounds 
% needed to decode a gene (alloweddiff).
%
% Debugging: the false positive rate
%
% Author: Sheel Shah
% Date: 8/4/2019

[consensuscell,copynumfinal ] = BarcodeNoMiji_v8( channels, points, hyb, barcodekey,radius,alloweddiff);

[dotlocations] = PointLocations(hyb, channels, points, consensuscell,copynumfinal, radius);

% added 01/08/20 by Yodai Takei to fix formating bug. convert [] to 0 for consistency.
for j = 1:size(barcodekey,1)
    if j > size(dotlocations,1) || isempty(dotlocations{j,1})
        dotlocations{j,1} = 0;
        dotlocations{j,2} = 0;
        dotlocations{j,3} = 0;
        dotlocations{j,4} = 0;
        dotlocations{j,5} = 0;
    end
end

[seeds] = numseeds(dotlocations);

end

