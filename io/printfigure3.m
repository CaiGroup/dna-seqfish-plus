function [] = printfigure3(rawImage, dotsLogical, dotsFiltered, saveFigPath, dotsThird)

    f = figure();
    % View image with dots: can use logFish or fish image
    imshow(max(rawImage, [], 3), 'DisplayRange', [0 2000], 'InitialMagnification', 'fit');
    hold on;
    [v2,v1] = find(max(dotsLogical, [], 3) == 1);
    scatter(v1(:), v2(:), 75);
    %[p2,p1] = find(max(dotsLogicalFiltered, [], 3) == 1);
    scatter(dotsFiltered(:,1), dotsFiltered(:,2), 50, 'mv');
    scatter(dotsThird(:,1), dotsThird(:,2), 50, 'g+');
    hold off;
    savefig(f, saveFigPath);
    close all;
end