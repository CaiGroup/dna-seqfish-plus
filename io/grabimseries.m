function [image, sizeC, sizeZ, physicalsizeX, physicalsizeY] = grabimseries(imagePath, position)
% create an image reader for images that are connected by metadata:
% Anderson data and automation scripts. bfopen opens all of the images for
% every position
% This function will only open the images for one position
%
% Code adapted from bfmatlab package
% Date: 5/7/2019

    r = bfGetReader(imagePath);
    sizeC = [];
    sizeZ = [];
    sizeT = [];
    numSeries = r.getSeriesCount();
    if numSeries == 1
        position = 0;
    end
    r.setSeries(position)
    pixelType = r.getPixelType();
    bpp = javaMethod('getBytesPerPixel', 'loci.formats.FormatTools', ...
    pixelType);
    bppMax = power(2, bpp * 8);
    numImages = r.getImageCount();
    imageList = cell(numImages, 2);
    colorMaps = cell(numImages);
    id = imagePath;
    s = 1; % never need to change this
    globalMetadata = r.getGlobalMetadata();
    for i = 1:numImages
        %{ 
        take out print statements
        %if mod(i, 72) == 1
            %fprintf('\n    ');
        %end
        %fprintf('.');
        %}
        arr = bfGetPlane(r, i);
        % retrieve color map data
        if bpp == 1
            colorMaps{s, i} = r.get8BitLookupTable()';
        else
            colorMaps{s, i} = r.get16BitLookupTable()';
        end
        warning_state = warning ('off');
        if ~isempty(colorMaps{s, i})
            newMap = single(colorMaps{s, i});
            newMap(newMap < 0) = newMap(newMap < 0) + bppMax;
            colorMaps{s, i} = newMap / (bppMax - 1);
        end
        warning (warning_state);
        % build an informative title for our figure
        label = id;
        if numSeries > 1
            seriesName = char(r.getMetadataStore().getImageName(s - 1));
            if ~isempty(seriesName)
                label = [label, '; ', seriesName];
            else
                qs = int2str(s);
                label = [label, '; series ', qs, '/', int2str(numSeries)];
            end
        end
        if numImages > 1
            qi = int2str(i);
            label = [label, '; plane ', qi, '/', int2str(numImages)];
            if r.isOrderCertain()
                lz = 'Z';
                lc = 'C';
                lt = 'T';
            else
                lz = 'Z?';
                lc = 'C?';
                lt = 'T?';
            end
            zct = r.getZCTCoords(i - 1);
            sizeZ = r.getSizeZ();
            if sizeZ > 1
                qz = int2str(zct(1) + 1);
                label = [label, '; ', lz, '=', qz, '/', int2str(sizeZ)];
            end
            sizeC = r.getSizeC();
            if sizeC > 1
                qc = int2str(zct(2) + 1);
                label = [label, '; ', lc, '=', qc, '/', int2str(sizeC)];
            end
            sizeT = r.getSizeT();
            if sizeT > 1
                qt = int2str(zct(3) + 1);
                label = [label, '; ', lt, '=', qt, '/', int2str(sizeT)];
            end
        end
        % save image plane and label into the list
        imageList{i, 1} = arr;
        imageList{i, 2} = label;
    end

    result{s, 1} = imageList;
    % extract metadata table for this series
    seriesMetadata = r.getSeriesMetadata();
    javaMethod('merge', 'loci.formats.MetadataTools', ...
    globalMetadata, seriesMetadata, 'Global ');
    dimensionOrder = r.getDimensionOrder();
    result{s, 2} = seriesMetadata;
    result{s, 3} = colorMaps;
    result{s, 4} = r.getMetadataStore();
    
    % get physical pixel size
    omeMeta = r.getMetadataStore(); %omXML = char(omeMeta.dumpXML());
    if ~isempty(omeMeta.getPixelsPhysicalSizeX(0))
        physicalsizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER).doubleValue();
        physicalsizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER).doubleValue();
    else
        physicalsizeX = [];
        physicalsizeY = [];
    end

    image = cell(1,sizeC);
    for i = 1:sizeC
        if strcmp(dimensionOrder, 'XYCZT')
            indices = i:sizeC:numImages;
        elseif strcmp(dimensionOrder, 'XYZCT')
            startIndex = sizeZ * (i - 1) + 1;
            endIndex = startIndex + sizeZ - 1;
            indices = startIndex:endIndex;
        end
        for j = indices
            image{i} = cat(3,image{i}, result{1,1}{j,1});
        end
    end
    r.close();
end