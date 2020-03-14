function tform = addtform(tformMoving, tformFixed)

    tform = tformFixed;
    tform.T(4,1) = tform.T(4,1) + tformMoving.T(4,1);
    tform.T(4,2) = tform.T(4,2) + tformMoving.T(4,2);
    tform.T(4,3) = tform.T(4,3) + tformMoving.T(4,3);

end