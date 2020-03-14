% Script to decode from the hpc

experimentDir = 'I:\2019-09-09-brain-rep2-2-DNAFISH';
%imDir = fullfile(experimentDir, 'aligned-images');
%I = cell(93,1);
position = 0;
channel = 1;

%{
for i = 0:92
imPath = fullfile(imDir, ['processed-hyb-p' num2str(position) '-h' num2str(i) '-c' num2str(channel) '.tif']);
I{i+1} = loadtiff(imPath);
end
%}

numRounds = 5;
numChannels = 16;
sqrtradius = 6;
segment = 'whole';
typedots = 'exons';
superres = 'radial';
experimentName='2019-09-09-brain-rep2-2-DNAFISH'; 
mkdir(pointsSaveDir);
minseeds = 3;

alloweddiff = 2;
for ch = 1
    [finalPosList, dotlocations, numpointconsensus, numdotlocations, numfinalpoints...
    ,numpointspercell, seeds, points] = processimages(experimentDir, experimentName, ...
    position, numRounds, numChannels, I(:,ch), segment, sqrtradius, typedots, superres, alloweddiff, ...
    ch, minseeds);
end

%{

% Test if the images are null

for i = 1:93
test = max(max(I{i}(:,:,20)));
fprintf('Max for hyb %.0f is %.0f\n', i, test);
end


%}