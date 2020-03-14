function tform = maketform2(x,y,z)
% get the offsets from the csv

    tform = affine3d(eye(4));
    tform.T(4,1) = x;
    tform.T(4,2) = y;
    tform.T(4,3) = z;

end