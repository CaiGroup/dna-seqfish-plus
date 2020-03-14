function chaTform = loadchabberationtform(chArray, refPoints)
    
    chaTform = cell(length(chArray),1);
    for i = chArray
        if i == 1
            chaTform{i} = affine3d(eye(4));
        else
            chaTform{i} = getglobaltform(refPoints(1).channels,refPoints(i).channels);
            
        end
    end
    
end