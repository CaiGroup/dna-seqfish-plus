function backcorrections = backcorrintron(pathName, numOfChannels)
% Get background over all positions for each folder
%switch channels because yodai switches channel 2 with 3

    MIJ.start;
    backcorrections = cell(1, numOfChannels);
    listing = dir([pathName filesep '*.tif']);   
    cy7back = [];
    cy5back = [];
    a594back = [];


    for i = 1:length(listing)
        path_to_fish = ['path=[' pathName filesep listing(i).name ']'];
        MIJ.run('Open...', path_to_fish);
        MIJ.run('Split Channels');
        for loop = 1:numOfChannels
            name = ['C' num2str(loop) '-' listing(i).name];
            im = uint16(MIJ.getImage(name));
            hybnum.color{loop} = im;
        end
        MIJ.run('Close All');

        cy7 = max(hybnum.color{1},[],3);
        cy5 = max(hybnum.color{2},[],3);
        a594 = max(hybnum.color{3},[],3);

        cy7back(:,:,i) = imopen(cy7,strel('disk',100));
        cy5back(:,:,i) = imopen(cy5,strel('disk',100));
        a594back(:,:,i) = imopen(a594,strel('disk',100));

    end

    cy7med = median(cy7back,3);
    cy5med = median(cy5back,3);
    a594med = median(a594back,3);

    cy7med = double(cy7med)/double(max(max(cy7med)));
    cy5med = double(cy5med)/double(max(max(cy5med)));
    a594med = double(a594med)/double(max(max(a594med)));

    backcorrections{1} = cy7med;
    %switch channels because yodai switches channel 2 with 3
    backcorrections{3} = cy5med;
    backcorrections{2} = a594med;

    figure;
    surf(double(cy7med(1:16:end,1:16:end))),zlim([0 1]);
    ax = gca;
    ax.YDir = 'reverse';

    figure;
    surf(double(cy5med(1:16:end,1:16:end))),zlim([0 1]);
    ax = gca;
    ax.YDir = 'reverse';

    figure;
    surf(double(a594med(1:16:end,1:16:end))),zlim([0 1]);
    ax = gca;
    ax.YDir = 'reverse';
    
    close all
    
    MIJ.exit;

    % figure;
    % surf(double(cy3med(1:8:end,1:8:end))),zlim([0 1]);
    % ax = gca;
    % ax.YDir = 'reverse';
    % 
    % figure;
    % surf(double(a488med(1:8:end,1:8:end))),zlim([0 1]);
    % ax = gca;
    % ax.YDir = 'reverse';


end