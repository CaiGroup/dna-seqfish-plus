function [] = savechannelsimagej(hybIms, savePath)
% saves the tif images from matlab to tiff images.
%
% Date: 10/22/2019
% Author: Nico Pierson
% Email: nicogpt@caltech.edu

    numChannels = length(hybIms);


    %% Need to save the images per barcode round for each channel
    Miji(false);
    numZSlice = size(hybIms{1}, 3);
    str = [];
    iter = 1;
    for ch = 1:numChannels
        namesh{iter} = ['C' num2str(1) '-'  num2str(ch) '.tif'];
        MIJ.createImage(namesh{iter}, hybIms{ch}, true);
        iter = iter + 1;
    end    
    iter = 1;
    for ch = 1:numChannels
        str = [str ' image' num2str(iter) '=C' num2str(1) '-' num2str(ch) '.tif'];
        iter = iter + 1;
    end

    
   
    try
        MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
        MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(numChannels) ' slices=' num2str(numZSlice) ' frames=1 display=Grayscale']);
        MIJ.run('Save', ['save=[' savePath ']']);
        MIJ.run('Close All')
    catch
        MIJ.exit;
        error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
    end

end