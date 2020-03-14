function F = gauss2dfrobenium(x,sup,ydata)
% 2D Gaussian Function that calculates the Frobenium norm to calculate the
% difference between the current matrix with the gaussian matrix.
%
% Inputs: (x vector [background, amplitude, x0, y0, std_x, std_y, theta],
% dimensions of the matrix, and the original matrix)
%
% Outputs: Frobenium norm
%
% To do: Calculate the area under the gaussian fit function
%
% Author: Mike Lawson
% Date: May 2018
% Modified By: Nico Pierson
% Date: 1/30/2019

% A numberic function for calculating what the minimum
% using this as the starting point : x0  = [min max-min(range) 4 4 1 1 45] 
% sup is the size
% ydate are the pixel values
back = x(1); % Set to background as minimum pixel value
amp = x(2); % Height of the peak after subtracting background
mx = x(3); % Initial x0 and y0 values. Why is it set to 4?
my = x(4); % So (4,4) in the 7x7 matrix is the initial start value in the center
sx = x(5); % the standard deviation of x and y is 1
sy = x(6);
th = x(7); % theta is 45 degrees

a = ((cosd(th)^2) / (2*sx^2)) + ((sind(th)^2) / (2*sy^2));
b = -((sind(2*th)) / (4*sx^2)) + ((sind(2*th)) / (4*sy^2));
c = ((sind(th)^2) / (2*sx^2)) + ((cosd(th)^2) / (2*sy^2));
for i = 1:sup(1) % sup is the size of the matrix: for each i and j
    for j = 1:sup(2)
        % Calculate a pixel value: minimum pixel value + the 2D gaussian
        % function of the 
        
        % Returns coordinates of gaussian 2D function
        % Add background of pixel value to all pixels
        q(i,j) = (back + ...
            amp*exp(-(a*(i-mx)^2 + 2*b*(i-mx)*(j-my) + c*(j-my)^2))  );
    end
end
% q is a 7x7 matrix of the gaussian normalized curve

% Subtract the pixels of the 7 x 7 matrix with the gaussian values, and
% finally calculate the frobenium norm for the matrix and return it

% Subtract the current matrix with the gaussian curve:
% values over will be positive and values less will be negative...
F = norm((ydata-q),'fro'); 
% Find the curve that fits the best according to the region of
% pixels

if mx > sup(1) || mx < 0 || my > sup(2) || my < 0 || back < 0 || amp < 0
    F = 1000000000000; % Set Forbenium to high if variables exceeds boundaries
end