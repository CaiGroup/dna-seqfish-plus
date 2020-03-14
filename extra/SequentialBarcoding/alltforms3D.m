function [tforms] = alltforms3D(PathName, colors, position, channels, reference)
% Switch channel 2 and 3 for chromatic aberrations

%   colors = # of channels in stack minus DAPI
%   channels = channels numbers that need to be aligned
%   reference = number of reference channel
    fld = pwd;
    cd(fld);
    %javaaddpath 'C:\Program Files\MATLAB\R2018a\java\mij.jar'
    %javaaddpath 'C:\Program Files\MATLAB\R2018a\java\ij.jar'
    MIJ.start;

    if isempty(PathName) 
        [FileName,PathName,FilterIndex] = uigetfile('.tif');
        path_to_fish = ['path=[' PathName FileName ']'];
    else
        FileName = [];
        listing = dir([PathName filesep '*.tif']);
        for i = 1:length(listing)
            if ~isempty(strfind(listing(i).name,num2str(position))) && ~isempty(strfind(listing(i).name,'ome')) 
                FileName = listing(i).name;
                break
            end
        end
        path_to_fish = ['path=[' PathName filesep FileName ']'];
    end

    MIJ.run('Open...', path_to_fish);
    MIJ.run('Split Channels');
    mkdir([PathName filesep 'channelsWBG']);
    for i = 1:colors
        name = ['C' num2str(i) '-' FileName];
        Images = uint16(MIJ.getImage(name));
        saveastiff(Images, [PathName filesep 'channelsWBG' filesep 'channel' num2str(i) '.tif'])
    end

    MIJ.run('Close All');
    MIJ.exit;
    for i = 1:colors
        tforms{i} = [];
    end
    for i = 1:length(channels)
        if i ~= 1
            [tforms{i}] = findtformV3([PathName filesep 'channelsWBG' filesep 'channel' num2str(channels(i)) '.tif'],[PathName filesep 'channelsWBG' filesep 'channel' num2str(reference) '.tif'],num2str(channels(i)));
        else
            tforms{i} = affine3d([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]);
        end
        close all;
    end
    
    % Switch channel 2 and 3 
    temp = tforms{2};
    tforms{2} = tforms{3};
    tforms{3} = temp;

    tforms(cellfun(@isempty,tforms)) = tforms(reference);
end


