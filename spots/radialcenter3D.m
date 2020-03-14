% radialcenter3D.m
%
% Calculates the center of a 3D intensity distribution -- 3D generalization
% of the 2D "radialcenter.m"
%
% Method: Considers lines passing through each half-pixel point with slope
% parallel to the gradient of the intensity at that point.  Considers the
% distance of closest approach between these lines and the coordinate
% origin, and determines (analytically) the origin that minimizes the
% weighted sum of these distances-squared.
% Considers 2D gradients in slices, and then gradient along z
% Uses 3x3 smoothing of gradients in xy, not in z.
% 
% See radialcenter.m, and notes of March 23-30, 2012 on the 3D
% generalization
%
% Inputs
%   I  : 3D intensity distribution (i.e. a grayscale image)
%   zscaleratio : ratio of scale (um/px e.g.) in z-direction relative to xy
%            plane (default 1)
%
% Outputs
%   rc :     the center of radial symmetry (3-element vector: [xc, yc, zc]
%            px, from px #1 = left/topmost/slice1 pixel
%            So a shape centered in the middle of a 2*N+1 x 2*N+1 x 2*Nz+1
%            rectangle will return a center value at x0=y0=N+1, z = Nz+1.
%            Note that y increases with increasing row number (i.e. "downward")
%   sigma  : Rough measure of the width of the distribution (sqrt. of the 
%            second moment of I - min(I)).  (Not determined by the fit.)
%
% Raghuveer Parthasarathy
% The University of Oregon
% March 30, 2012 
% last modified March 30, 2012 
% Copyright 2012, Raghuveer Parthasarathy, The University of Oregon

function [rc, sigma] = radialcenter3D(I, zscaleratio)
%#codegen
% Default parameters
if ~exist('zscaleratio', 'var') || isempty(zscaleratio)
    zscaleratio = 1.0;  % ratio of z-direction scale (distance/px) to xy-drection scale
end

% Number of grid points
[Ny, Nx, Nz] = size(I);

% grid coordinates are -n:n, where Nx (or Ny) = 2*n+1
% grid midpoint coordinates are -n+0.5:n-0.5;
% note that in radialcenter.m I avoid repmat for speed; here I'll use it
% for the z-repetition, for simplicity.
xm_onerow = -(Nx-1)/2.0+0.5:(Nx-1)/2.0-0.5;
xm2D = xm_onerow(ones(Ny-1, 1), :);
xm = repmat(xm2D,[1 1 Nz-1]);
ym_onecol = (-(Ny-1)/2.0+0.5:(Ny-1)/2.0-0.5)';  % Note that y increases "downward"
ym2D = ym_onecol(:,ones(Nx-1,1));
ym = repmat(ym2D,[1 1 Nz-1]);
zm_onedepth = reshape((-(Nz-1)/2.0+0.5:(Nz-1)/2.0-0.5), 1 , 1, Nz-1);  
zm = repmat(zm_onedepth, [Nx-1 Ny-1 1]);

%% 
% For each slice, calculate the gradient in the x-y plane
% Nx-1 x Ny-1 arrays: dImag2D2, dI2Dx, dI2Dy
dI2Dx = zeros(Ny-1, Nx-1, Nz);
dI2Dy = zeros(Ny-1, Nx-1, Nz);
for j=1:Nz
    % Calculate derivatives along 45-degree shifted coordinates (u and v)
    % Note that y increases "downward" (increasing row number) -- we'll deal
    % with this when calculating "m" below.
    dIdu = I(1:Ny-1,2:Nx,j)-I(2:Ny,1:Nx-1,j);
    dIdv = I(1:Ny-1,1:Nx-1,j)-I(2:Ny,2:Nx,j);
    
    % Smoothing --
    h = ones(3)/9;  % simple 3x3 averaging filter
    fdu = conv2(dIdu, h, 'same');
    fdv = conv2(dIdv, h, 'same');
    % Note that we need a 45 degree rotation of
    % the u,v components to express the slope in the x-y coordinate system.
    % The negative sign "flips" the array to account for y increasing
    % "downward"
    dI2Dx(:,:,j) = fdu-fdv;
    dI2Dy(:,:,j) = -(fdv+fdu);
    
    % Delete various slope=undefined checks found in radialcenter.m
end

%% 
% For each pair of slices, calculate the gradient in z, and the average of
% the xy-plane gradients
dIx = zeros(Ny-1, Nx-1, Nz-1);
dIy = zeros(Ny-1, Nx-1, Nz-1);
dIz = zeros(Ny-1, Nx-1, Nz-1);
for j=1:Nz-1
    dIx(:,:,j) = 0.5*(dI2Dx(:,:,j) + dI2Dx(:,:,j+1)); 
    dIy(:,:,j) = 0.5*(dI2Dy(:,:,j) + dI2Dy(:,:,j+1)); 
    dIz(:,:,j) = (1/zscaleratio)*0.25*(I(1:Ny-1,1:Nx-1,j+1) - I(1:Ny-1,1:Nx-1,j) + ...
                       I(1:Ny-1,2:Nx,j+1) - I(1:Ny-1,2:Nx,j) + ...
                       I(2:Ny,1:Nx-1,j+1) - I(2:Ny,1:Nx-1,j) + ...
                       I(2:Ny,2:Nx,j+1) - I(2:Ny,2:Nx,j));
end
dImag2 = dIx.*dIx + dIy.*dIy + dIz.*dIz;
sqdImag2 = sqrt(dImag2);
mx = dIx ./ sqdImag2; % x-component of unit vector in the gradient direction
my = dIy ./ sqdImag2; % x-component of unit vector in the gradient direction
mz = dIz ./ sqdImag2; % x-component of unit vector in the gradient direction

% Weighting: weight by square of gradient magnitude and inverse 
% distance to gradient intensity centroid.
sdI2 = sum(dImag2(:));
prodx = dImag2.*xm;
prody = dImag2.*ym;
prodz = dImag2.*zm;
xcentroid = sum(prodx(:))/sdI2;
ycentroid = sum(prody(:))/sdI2;
zcentroid = sum(prodz(:))/sdI2;
w  = dImag2./sqrt((xm-xcentroid).*(xm-xcentroid)+(ym-ycentroid).*(ym-ycentroid)+(zm-zcentroid).*(zm-zcentroid));  

% least-squares minimization to determine the translated coordinate
% system origin (xc, yc, zc) such that gradient lines have
% the minimal total distance^2 to the origin:
% MOVE to function lsradialcenterfit3D (below)
% rc = [xc yc zc]
[rc] = lsradialcenterfit3D(xm, ym, zm, mx, my, mz, w);



%%
% Return output relative to upper left coordinate, slice 1
rc(1) = rc(1) + (Nx+1)/2.0;
rc(2) = rc(2) + (Ny+1)/2.0;
rc(3) = rc(3) + (Nz+1)/2.0;

% A rough measure of the particle width.
% Not at all connected to center determination, but may be useful for tracking applications; 
% could eliminate for (very slightly) greater speed
Isub = I - min(I(:));
[px,py,pz] = meshgrid(1:Nx,1:Ny,1:Nz);
xoffset = px - rc(1);
yoffset = py - rc(2);
zoffset = pz - rc(3);
r2 = xoffset.*xoffset + yoffset.*yoffset + zoffset.*zoffset;
sigma = sqrt(sum(sum(sum(Isub.*r2)))/sum(Isub(:)))/2;  % second moment is 2*Gaussian width
    

%%

    function [rc] = lsradialcenterfit3D(xm, ym, zm, mx, my, mz, wk)
        % least squares solution to determine the radial symmetry center
        
        % inputs mx, my, mz, w are defined on a grid
        % w are the weights for each point
        mz2pmy2w = (mz.*mz + my.*my).*wk;
        mx2pmz2w = (mx.*mx + mz.*mz).*wk;
        my2pmx2w = (my.*my + mx.*mx).*wk;
        mxyw = mx.*my.*wk;
        mxzw = mx.*mz.*wk;
        myzw = my.*mz.*wk;
        A = [sum(mz2pmy2w(:))  -sum(mxyw(:))      -sum(mxzw(:)); ...
             -sum(mxyw(:))      sum(mx2pmz2w(:))  -sum(myzw(:)); ...
             -sum(mxzw(:))     -sum(myzw(:))       sum(my2pmx2w(:))];
        xmz2pmy2w = xm.*mz2pmy2w;
        ymx2pmz2w = ym.*mx2pmz2w;
        zmy2pmx2w = zm.*my2pmx2w;
        ymxyw = ym.*mxyw;
        zmxzw = zm.*mxzw;
        xmxyw = xm.*mxyw;
        zmyzw = zm.*myzw;
        xmxzw = xm.*mxzw;
        ymyzw = ym.*myzw;
        rhs = [sum(xmz2pmy2w(:) - ymxyw(:)     - zmxzw(:)); ...
               sum(-xmxyw(:)    + ymx2pmz2w(:) - zmyzw(:)); ...
               sum(-xmxzw(:)    - ymyzw(:)     + zmy2pmx2w(:))];
        % best fit minimal distance, relative to image center
        rc = A\rhs;
    end

end
