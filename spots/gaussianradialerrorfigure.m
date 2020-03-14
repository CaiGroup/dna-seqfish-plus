function fig = gaussianradialerrorfigure(image, refPoints, fixedPoints, gaussPoints, radialPoints, maxDisplayRange, hyb, ch, tempSavePath)

    % colors for figures
    %colorR = {'m', [0.4660 0.6740 0.1880], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], 'c', [0.6350 0.0780 0.1840], 'b', 'r','y', [0 0.4470 0.7410], [0.3 0.4 0.1], [0.9 0.7 0.9], [0.9 1.0 0.7], [0.7 0.5 0.1], [0.5 0.7 0.8]};
    %% Figure for all three points
    close all
    fig = figure;
    imshow(max(image, [], 3), 'DisplayRange', [0 maxDisplayRange]);
    hold on
    s1 = scatter(fixedPoints(:,1), fixedPoints(:,2), 0.1, 'r+', 'LineWidth', 0.01);
    s1.MarkerFaceAlpha = 0.01;
    s2 = scatter(gaussPoints(:,1), gaussPoints(:,2), 0.1, 'g+', 'LineWidth', 0.01);
    s2.MarkerFaceAlpha = 0.01;
    s3 = scatter(radialPoints(:,1), radialPoints(:,2), 0.1, 'm+', 'LineWidth', 0.01);
    s3.MarkerFaceAlpha = 0.01;
    
    scatter(fixedPoints(:,1), fixedPoints(:,2), 5, 'ro', 'LineWidth', 0.01);
    scatter(gaussPoints(:,1), gaussPoints(:,2), 5, 'gs', 'LineWidth', 0.01);
    scatter(radialPoints(:,1), radialPoints(:,2), 5, 'm^', 'LineWidth', 0.01);
    legend('none','gauss','radial', 'none','gauss','radial', 'Location','southeast');
    hold off

    
    filename = fullfile(tempSavePath, ['temp-hyb' num2str(hyb) 'ch' num2str(ch) '-figure.tif']);
    %saveas(gcf, filename, 'tiffn');
    print(filename, '-dtiffn', '-r800');
    drawnow
    close all;
end