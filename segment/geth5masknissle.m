function I = geth5masknissle(h5filepath)
% function geth5mask retrieves the mask from an .h5 file usually created
% from ilastik software

    infoh5File = h5info(h5filepath, 'TextEncoding', 'UTF-8');
    dataset = [infoh5File.Groups(2).Name '/' infoh5File.Groups(2).Datasets.Name];
    data = h5read(h5filepath, dataset);
    if size(data, 3) > 1
        imageSegment = permute(data, [2 1 3]);
    elseif size(data, 3) == 1
        imageSegment = permute(data, [2 1]);
    end
    I = imresize(imageSegment, 4, 'nearest');

end