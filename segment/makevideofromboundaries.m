function [] = makevideofromboundaries(boundaries, nissleImage, numZSlices, saveVideoFilePath, upperImageThreshold)
% function that makes video from the boundaries organized as zslices and
% cell with vertices 
% 
% Date: June 3, 2019
    

    for z = 1:numZSlices
        fprintf('Geting Zslice: %.0f\n', z);
        rawImage = nissleImage(:,:,z);

        % draw figure
        figure;
        imshow(rawImage, [0 upperImageThreshold]); % get the original image: % 8000 for Control1 %1800 for Control2 %40000 for aggression and mating
        hold on
        for c = 1:length(boundaries{z}) % fix the structure later
            if ~isempty(boundaries{z}{c})
                pgon = polyshape(boundaries{z}{c});
                [xlabel,ylabel] = centroid(pgon);
                labelCell = num2str(c);
                labelBoldCell = ['\bf' num2str(c)];
                % plot the borders and the associated cell number 
                plot(boundaries{z}{c}(:,1),boundaries{z}{c}(:,2), 'm.', 'MarkerSize',2);
                text((xlabel-10),ylabel,labelBoldCell,'Color','g','FontSize',9); 
                text((xlabel-10),ylabel,labelCell,'Color','k','FontSize',8);  
            end
        end
        legend('nissle', 'Location', 'southeast');
        hold off
        F(z) = getframe(gcf);
        drawnow

    end

    % create the video writer with 1 fps
    saveVideoPath = [saveVideoFilePath '.avi'];
    writerObj = VideoWriter(saveVideoPath);
    writerObj.FrameRate = 10;
      % set the seconds per image
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);

    close all;
end