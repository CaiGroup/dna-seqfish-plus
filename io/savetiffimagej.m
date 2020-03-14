function pass = savetiffimagej()
% savetiffimagej saves the images in the dimension order xyzct to view 3D
% images properly
%
% Dependencies: run imageJ beforehand
%
% Author: Nico Pierson
% Date: 4/11/2019

% what is needed: numSerialHybs, barcodeNumber, channel for naming, images,
% zslices, savePath, image{numSerialHybs}

            Miji(false);
            for f = 1:numberOfSerialHybs
                    namesh{f} = ['C' num2str(f) '-'  num2str(barcode) '.tif'];
                    MIJ.createImage(namesh{f}, dapi{f}, true);
            end

            str = [];
            for f = 1:numberOfSerialHybs
                    str = [str ' image' num2str(f) '=C' num2str(f) '-' num2str(barcode) '.tif'];
            end


            try
                MIJ.run('Concatenate...', ['  title=[Concatenated Stacks] open' str]);
                MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(numberOfSerialHybs) ' slices=' num2str(zSlice) ' frames=1 display=Grayscale']);
                savePath = fullfile(saveDir, ['dapi_pos' num2str(position) '_bar' num2str(barcode)]);
                MIJ.run('Save', ['save=[' savePath '.tif' ']']);
                MIJ.run('Close All');
                MIJ.exit;
            catch
                MIJ.exit;
                error('MIJ exited incorrectly: most likely caused by out of memory in the java heap\n');
            end
            



end