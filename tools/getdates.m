function dateS = getdates
%Load the range of dates defining solar cycles and the range of dates used
%in the BTSI study. dateS.all is the range used in analysis
dateS.cycles = [juliandate(datetime(1976,03,01)) juliandate(datetime(1986,8,31));
juliandate(datetime(1986,9,01)) juliandate(datetime(1996,7,31));
juliandate(datetime(1996,8,01)) juliandate(datetime(2008,11,31));
juliandate(datetime(2008,12,01)) juliandate(datetime(2019,11,31));
juliandate(datetime(2019,12,01)) juliandate(datetime(2030,11,31))];

dateS.allold = [juliandate(datetime(1976,03,01)) juliandate(datetime(2020,05,31))];
dateS.all = [juliandate(datetime(1978,11,1)) juliandate(datetime(2021,11,30))];
end