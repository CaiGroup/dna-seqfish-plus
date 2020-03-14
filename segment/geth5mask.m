function I = geth5mask(h5filepath)
% function geth5mask retrieves the mask from an .h5 file usually created
% from ilastik software

    infoh5File = h5info(h5filepath, 'TextEncoding', 'UTF-8');
    dataset = ['/' infoh5File.Datasets.Name];
    data = h5read(h5filepath, dataset);
    sizeImage = size(data);
    if size(sizeImage, 2) == 4
        imageSegment = reshape(data, sizeImage(2:4));
        imageSegment = permute(imageSegment, [2 1 3]);
    elseif size(sizeImage, 2) == 3
        imageSegment = reshape(data, sizeImage(2:3));
        imageSegment = permute(imageSegment, [2 1]);
    end
    I = imresize(imageSegment, 4, 'nearest');

end