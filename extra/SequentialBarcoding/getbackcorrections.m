function backcorrections = getbackcorrections(pathName, numOfChannels)
% getbackcorrections grabs all blank images and makes a normalized
% background that corrects for irregular illumination by returning a matrix
% of values between 0 and 1 to correct this.
%
% Author: Nico Pierson
% Email: nicogpt@caltech.edu
% Date: 1/22/2019

    %% Get the background corrections
    % Use Blanks to calculate corrections for each channel (2048 x 2048)
    backCorrPath = [pathName filesep 'Blanks_Hyper_Swapped'];
    backcorrections = backcorrintron(backCorrPath, numOfChannels);
end