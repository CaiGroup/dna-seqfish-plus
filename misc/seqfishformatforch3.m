function pointsch3final = seqfishformatforch3(pointsch, folderArray, numRounds, numCh3pointshyb)
numHybs = length(folderArray);
numChannels = round(numHybs / numRounds);% barcoding pseudo-channels.
corrpoints = cell(numHybs,1);

for ro = 1:numRounds
    for bch = 1:numChannels
        corrpoints{(ro-1)*numChannels+bch,1}.channels = pointsch{ro,1}(bch).channels;
        corrpoints{(ro-1)*numChannels+bch,1}.intensity = pointsch{ro,1}(bch).intensity;
    end
end

pointsch3final = corrpoints(1:numCh3pointshyb,1);

end