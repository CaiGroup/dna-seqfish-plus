<<<<<<< HEAD
function pass = makefigvideo(numHybs, numCh, tempSavePath, saveVideoPathName)
% make a video from cell arra of figures


    pass = 0;
    % 10 different colors for different regions
    %colorR = {'m', [0.4660 0.6740 0.1880], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], 'c', [0.6350 0.0780 0.1840], 'b', 'r','y', [0 0.4470 0.7410], [0.3 0.4 0.1], [0.9 0.7 0.9], [0.9 1.0 0.7], [0.7 0.5 0.1], [0.5 0.7 0.8]};
    
    for ch = 2:numCh
        % create the video writer with 1 fps
        saveVideoPath = [saveVideoPathName '-ch' num2str(ch) '.avi'];

        writerObj = VideoWriter(saveVideoPath);
        writerObj.FrameRate = 10;% set the seconds per image
        % open the video writer
        open(writerObj);


        for hyb = 1:numHybs   
            filename = fullfile(tempSavePath, ['temp-hyb' num2str(hyb) 'ch' num2str(ch) '-figure.tif']);
            img = imread(filename);
            delete(filename);
            writeVideo(writerObj, img);
        end
        % close the writer object
        close(writerObj);
        close all;
    end
=======
function pass = makefigvideo(figures, saveVideoPathName)
% make a video from cell array of figures


    pass = 0;
    numFigs = length(figures);
    % 10 different colors for different regions
    %colorR = {'m', [0.4660 0.6740 0.1880], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], 'c', [0.6350 0.0780 0.1840], 'b', 'r','y', [0 0.4470 0.7410], [0.3 0.4 0.1], [0.9 0.7 0.9], [0.9 1.0 0.7], [0.7 0.5 0.1], [0.5 0.7 0.8]};
    
    
    for f = 1:numFigs

        %{
        lgd = legend(h, legend_str, 'Location', 'southeast');
        lgd.FontSize = 6;
        lgd.NumColumns = 3;
        %}
        %f.WindowState = 'maximized';
        hold off
        %F(z) = getframe(f);
        %set(gcf,'units','normalized','outerposition',[0 0 1 1])
        %truesize(f, [2048 2048]);
        %F(z) = getframe(f);
        filename = ['temp-' num2str(f) '-image'];
        %saveas(gcf, filename, 'tiffn');
        print(filename, '-dtiffn', '-r500');
        drawnow
    end

    % create the video writer with 1 fps
    saveVideoPath = [saveVideoPathName '.avi'];
    
    writerObj = VideoWriter(saveVideoPath);
    writerObj.FrameRate = 10;% set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    %for i=1:length(F)
    for i = 1:numFigs    
        filename = ['temp-' num2str(i) '-image.tif'];
        img = imread(filename);
        delete(filename);
        writeVideo(writerObj, img);
    end
    % close the writer object
    close(writerObj);
    close all;

>>>>>>> d54b3621b842140e74bdf8231a9fcd8d9e2bc858
    pass = 1;
end