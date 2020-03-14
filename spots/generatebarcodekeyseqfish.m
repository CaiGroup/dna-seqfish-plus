function barcodekey = generatebarcodekeyseqfish()
% function to generate all the barcodekeys for seqfish, using a one round
% barcode error correction.
%
% structure fields: 'names' as number, and 'barcode' with the barcodekey
%
% Date: 9/9/2019

    bcI = 0;
    for i = 1:20
        for j = 1:20
            for k = 1:20
                bcI = bcI + 1;
                bc(bcI,:) = [i j k rem(i+j-k+20,20)];
                if rem(i+j-k+20,20) == 0
                    bc(bcI,:) = [i j k 20];
                end
            end
        end
    end
    barcodekey.barcode = bc;
    for i = 1:length(bc)
        barcodekey.names{i} = sprintf('%.0f',i);
    end
    barcodekey.names = barcodekey.names';

end