function [pointsSR, intensitySR, sigmaSR] = SuperResPoints(points,image,xyPixSize,zPixSize, varargin)
%#codegen


 %% Set up optional Parameters
    numvarargs = length(varargin);
    if numvarargs > 1
        error('myfuns:SuperResPoints:TooManyInputs', ...
            'requires at most 1 optional inputs');
    end
    optargs = {true}; 
    optargs(1:numvarargs) = varargin;
    [filteroption] = optargs{:};

pointsSR = points;
intensitySR = zeros(size(points, 1), 1);
sigmaSR = zeros(size(points, 1), 1);

[y,x,z] = size(image);

for k = size(points,1):-1:1
    if points(k,2)-1 == 0 || points(k,2)+1 > y ||...
       points(k,1)-1 == 0 || points(k,1)+1 > x ||...
       points(k,3)-1 == 0 || points(k,3)+1 > z
        if filteroption
            pointsSR(k,:) = [];
            intensitySR(k) = [];
            sigmaSR(k,:) = [];
        else
            pointsSR(k,:) = [1 1 1];
            intensitySR(k) = 0;
            sigmaSR(k,:) = 0;
        end
    else   
        I = image(points(k,2)-1:points(k,2)+1,...
                               points(k,1)-1:points(k,1)+1,...
                               points(k,3)-1:points(k,3)+1);
        [rc, sigmaSR(k)] = radialcenter3D(double(I), zPixSize/xyPixSize);
        rc([1 2]) = rc([2 1]);
        pointsSR(k,:) =((rc-[2;2;3]+points(k,:)').*[xyPixSize;xyPixSize;zPixSize])';
        xv = round(pointsSR(k,1));
        yv = round(pointsSR(k,2));
        zv = round(pointsSR(k,3));
        if zv > z
            zv = z;
        elseif zv < 1
            zv = 1;
        end
        if yv > y
            yv = y;
        elseif yv < 1
            yv = 1;
        end
        if xv > x
            xv = x;
        elseif xv < 1
            xv = 1;
        end
        if (isnan(xv) || isnan(yv) || isnan(zv))||(isempty(xv) || isempty(yv) || isempty(zv))
            intensitySR(k) = 0;
        else
            intensitySR(k) = image(yv,xv,zv); % replace with sigma
        end
    end
end


