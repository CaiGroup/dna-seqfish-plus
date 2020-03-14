function minThreshold = findminthreshold(x1, x2, numIters, target, I, typedots)
% takes interval and divides the interval into certain number of bins,
% iterates over the function and finds the minimum value from the function.
%
% Date: 8/12/2019
    
    highThreshold = 9999999999;
    numIters1 = 5;
    stepInterval = (x2 - x1) / (numIters1-1);
    binArray = x1:stepInterval:x2;
    d = ones(length(binArray),1) * highThreshold;
    prevMin = highThreshold;
    x2Index = x2;
    x1Index = x1;
    logImage = [];
    regMax = [];
    
    % first interval
    for i = 1:length(binArray)
        threshold = binArray(i);
        % put same processed image in if using multiple times
        if i == 1
            [d(i), logImage, regMax] = getthreshold4numpoints(threshold, target, I, typedots,logImage, regMax);
        else
            [d(i), ~] = getthreshold4numpoints(threshold, target, I, typedots, logImage, regMax);
        end
        if prevMin <= d(i)
            if i <= 2
                x1Index = binArray(1);
            else
                x1Index = binArray(i-2);
            end
            x2Index = binArray(i);
            break;
        end
        prevMin = d(i);
    end
    
    
    % check the second interval
    stepInterval2 = (x2Index - x1Index) / (numIters-1);
    binArray2 = x1Index:stepInterval2:x2Index;
    d2 = ones(length(binArray2),1) * highThreshold;
    prevMin2 = highThreshold;
    for i = 1:length(binArray2)
        threshold = binArray2(i);
        [d2(i), ~] = getthreshold4numpoints(threshold, target, I, typedots, logImage, regMax);
        % prevMin9999 and d = 1000; prev=1000 d = 2000
        if prevMin2 <= d2(i)
            minIndex = i - 1;
            break;
        end
        prevMin2 = d2(i);
    end
    
    
    minThreshold = binArray2(minIndex);
end