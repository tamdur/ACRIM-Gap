function [yrVal,yr] = monthlytoyearly(datetimes,vals)
%make a vector for every month in the dataset
yrInit=datetimes.Year;
stYr = min(yrInit); endYr = max(yrInit);
yr=(stYr:endYr)';
yrVal=NaN(length(yr),size(vals,2));
for iB = 1:length(yr)
    iYr = yrInit == yr(iB);
    if sum(iYr) > 0
        yrVal(iB,:) = nanmean(vals(iYr,:),1);
    end
end
end
