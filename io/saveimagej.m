function pass = saveimagej(I, savePath)
% function saveimagej saves the images from a cellarray and returns the
% pass value if function is completed
%
% Author: Nico Pierson
% Date: 4/12/2019

% add a way to check the path of Fiji.app

        Miji(false);
        pass = 0;
        zSlice = size(I{1}, 3);
        numChannels = length(I);
        for c = 1:numChannels
            namesh{c} = ['C' num2str(c) '-'  num2str(1) '.tif'];
            MIJ.createImage(namesh{c}, I{c}, true);
        end

        str = [];
        for c = 1:numChannels
                str = [str ' image' num2str(c) '=C' num2str(c) '-' num2str(1) '.tif'];
        end

        
        try
            MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
            MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(numChannels) ' slices=' num2str(zSlice) ' frames=1 display=Grayscale']);
            MIJ.run('Save', ['save=[' savePath ']']);
            MIJ.run('Close All');
            MIJ.exit;
        catch
            MIJ.exit;
            error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
        end
        pass = 1;
        
end