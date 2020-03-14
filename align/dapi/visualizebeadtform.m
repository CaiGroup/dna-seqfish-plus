function I = visualizebeadtform(beadImDir, position, beadtform)
% visualizes bead tform using beaad images as a reference
%
% Assumes the images are named 'MMStack_Pos0.ome.tif' where position number
% varies
%
% Dependencies; 
% Fiji Path
% addpath('C:\Users\Long Cai - 2\Desktop\Fiji.app\scripts\', '-end'); 
% bfmatlab path to open images
% addpath('C:\github\streamline-seqFISH\src\process_with_beads\bfmatlab\');

    
    %% get the images
    beadImPath = fullfile(beadImDir, ['MMStack_Pos' num2str(position) '.ome.tif']);
    [I, sizeC, sizeZ, ~] = grabimseries(beadImPath, position);
    
    %% Use tform
    for ch = 1:sizeC
    	I{ch} = imwarp(I{ch}, beadtform{ch}, 'OutputView', imref3d(size(I{ch})));
    end


end