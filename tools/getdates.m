function dateS = getdates
%Load the range of dates defining solar cycles and the range of dates used
%in the BTSI study. dateS.all is the range used in analysis
dateS.cycles = [juliandate(datetime(1976,03,01)) juliandate(datetime(1986,8,31));
juliandate(datetime(1986,9,01)) juliandate(datetime(1996,7,31));
juliandate(datetime(1996,8,01)) juliandate(datetime(2008,11,31));
juliandate(datetime(2008,12,01)) juliandate(datetime(2019,11,31));
juliandate(datetime(2019,12,01)) juliandate(datetime(2030,11,31))];

dateS.acrim = [juliandate(datetime(1989,06,02)) juliandate(datetime(1991,10,03))];
dateS.acrimplusfive = [juliandate(datetime(1989-5,06,02)) juliandate(datetime(1991+5,10,03))];
end