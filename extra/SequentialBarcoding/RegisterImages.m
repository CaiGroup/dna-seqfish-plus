function [hybnum, tformsAll] = RegisterImages(PathName, posnum,regimage,channelsall, prehybnum, regis, seqDataFolder)
% function registers the images
%
% Author: Yodai Takei and Sheel Shah
% Email: ytakei@caltech.edu
% Date: 2018
% Modified By: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/23/2019

%% Get the Dapi and Save it in the path
MIJ.start;
path_to_fish = ['path=[' PathName filesep 'Pos' num2str(posnum) filesep '1.tif' ']'];
MIJ.run('Open...', path_to_fish);
MIJ.run('Duplicate...', 'duplicate channels=4'); % The number should be changed based on the total color
MIJ.run('Save', ['save=[' PathName filesep 'Pos' num2str(posnum) filesep 'FirstHybDAPI.tif' ']']);
MIJ.run('Close All');
MIJ.exit;

%% Register the Images
fld = pwd;
MIJ.start;
cd(fld)

rPath = [PathName filesep 'Pos' num2str(posnum) filesep regimage];
path_to_fish = ['path=[' rPath ']'];
MIJ.run('Open...', path_to_fish);
regimage = MIJ.getImage(regimage);%modified for single image
MIJ.run('Close All')
MIJ.exit
registered{1} = regimage;

for i = 1:length(prehybnum)
    channels = channelsall{i};
    [optimizer, metric] = imregconfig('monomodal');
    tform = imregtform(regis{i}, regimage, 'translation' ,optimizer, metric);
    tformsAll{i} = tform;
    registered{i+1} = imwarp(regis{i},tform,'OutputView',imref3d(size(regimage)));
 
    for j = 1:length(channels)
        hybnum(i).color{j} = imwarp(prehybnum(i).color{j}, tform, 'OutputView', imref3d(size(regimage)));
    end
end
clearvars prehybnum regis; % clear space

fld = pwd;
MIJ.start;
cd(fld);
za = cellfun(@(x) size(x,3),registered);
for i = 1:length(registered)
    MIJ.createImage(num2str(i),registered{i}(:,:,1:min(za)),true);
end

hordor = [];

for i = 1:length(registered)
    temp = ['image' num2str(i) '=' num2str(i)];
    hordor = [hordor ' ' temp];
end

str = ['title=[Concatenated Stacks] ' hordor];
MIJ.run('Concatenate...', str);
MIJ.run('Stack to Hyperstack...', ['order=xyzct channels=' num2str(length(registered)) ' slices=' num2str(min(za)) ' frames=1 display=Grayscale']);

MIJ.run('Save', ['save=[' PathName filesep 'Pos' num2str(posnum) filesep seqDataFolder filesep 'SequentialRegistrationCheck.tif' ']']);

MIJ.run('Close All')

MIJ.exit

save([PathName filesep 'Pos' num2str(posnum) filesep seqDataFolder filesep ...
    'Pos' num2str(posnum) 'SequentialHybsNewRegis.mat'], 'hybnum', 'tformsAll','registered','-v7.3')