function I = subtractmulticut(multicutMask, simpleMask)
% function subractmulticut subtracts the multicut mask from a simple
% segmentation mask.

    revisedSimpleMask = (double(simpleMask) - 2) .* -1;
    I = multicutMask .* uint32(revisedSimpleMask);


end