function [] = makevideobarcode(image, points, numRounds, numChannels, threshold, displayRange, savePath)
% function makes a video of the points from a barcode experiment
%
% Date: 8/22/2019

    size = 4;
    for r = 1:numRounds
        for c = 1:numChannels
            hybIndex = (r - 1) * numChannels + c; 
            p = points{r}(c).channels;
            imshow(max(image{hybIndex},[],3), 'DisplayRange', displayRange); % get the original image: % 8000 for Control1 %1800 for Control2 %40000 for aggression and mating
            hold on
            scatter(p(:,1), p(:,2), size, 'mo', 'LineWidth', 0.05);
            legend(['Round-' num2str(r) '-Ch' num2str(c) '-Threshold:' num2str(threshold(r,c))], 'Location', 'southeast');
            hold off
            
            
            filename = [savePath '-R' num2str(r) '-C' num2str(c) '-imageSlice.tif'];
            print(filename, '-dtiffn', '-r300');
            %F(hybIndex) = getframe(gcf);
            drawnow
        end
                


    end

    % create the video writer with 1 fps
    saveVideoPath = [savePath '.avi'];
    writerObj = VideoWriter(saveVideoPath);
    writerObj.FrameRate = 10;
      % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for r = 1:numRounds
        for c = 1:numChannels 
            filename = [savePath '-R' num2str(r) '-C' num2str(c) '-imageSlice.tif'];
            img = imread(filename);
            delete(filename);
            writeVideo(writerObj, img);%frame);
        end
    end
    % close the writer object
    close(writerObj);

    close all;
end