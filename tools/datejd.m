function date = datejd(JD)
%DATEJD give datetime from Julian Date
% Ted Amdur 12/13/19
date = datetime(JD,'ConvertFrom','JulianDate');
end

