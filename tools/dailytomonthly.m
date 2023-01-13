function [mthVal,mthDate,nObs] = dailytomonthly(jd,vals,idvals,id)
% DAILYTOMONTHLY convert vector of daily data with attached julian dates to
% a vector of averaged monthly values
if nargin < 3
    idvals = [];
    id = [];
end

%make a vector for every month in the dataset
stJD = min(jd); endJD = max(jd);
stDate = datejd(stJD); endDate = datejd(endJD);
dateM = stDate;
while juliandate(dateM(end)) < juliandate(endDate)
    dateM = [dateM; dateM(end) + calmonths(1)];
end
dateMS = dateshift(dateM,'start','month');
dateME = dateshift(dateM,'end','month');
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