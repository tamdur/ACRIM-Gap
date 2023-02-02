function [mthVal,mthDate,nObs] = dailytomonthly(jd,vals,dateMS,dateME,idvals,id)
% DAILYTOMONTHLY convert vector of daily data with attached julian dates to
% a vector of averaged monthly values
if nargin < 5
    idvals = [];
    id = [];
end

jdS = juliandate(dateMS);
jdE = juliandate(dateME);
mthVal = NaN(length(jdS),1);
for iB = 1:length(jdS)
    if isempty(id)
        iMth = jd >= jdS(iB) & jd < jdE(iB);
    else
        iMth = jd >= jdS(iB) & jd < jdE(iB) & idvals == id;
    end
    if any(iMth)
        mthVal(iB) = nanmean(vals(iMth));
    end
    nObs(iB)=sum(iMth);
end
mthDate = mean([dateME dateMS],2);
%Remove entries with NaN values
inan=~isnan(mthVal);
mthVal=mthVal(inan);mthDate=mthDate(inan);nObs=nObs(inan);
end