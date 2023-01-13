function [xOut,offset] = meaninterval(dateIn,xIn,intStart,intEnd)
%Subract away the mean between intStart and intEnd, inclusive
% INPUTS:
%           dateIn: Nx1 datetime vector for xIn
%           xIn: Nx1 vector of time series values passed in
%           intStart: starting year
%           intEnd: ending year
% OUTPUTS:
%           xOut: xIn, demeaned between intStart and intEng

ind = dateIn.Year >= intStart & dateIn.Year <= intEnd;
offset = nanmean(xIn(ind));
xOut = xIn - offset;
end

