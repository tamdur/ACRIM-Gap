%Convert all satellite records to JD, save in structure array
%
% Ted Amdur
% 08/17/2018

%Import each respective text file

%--------------------------------------------------
% HF/NIMBUS
%--------------------------------------------------
nimbusimport
sat(1).ID = 'NIMBUS-7/HF';
sat(1).JD = juliandate(datetime(NIMBUS.year,NIMBUS.month,NIMBUS.day));
sat(1).TSI = NIMBUS.TSI;
%--------------------------------------------------
% ACRIM1/SMM
%--------------------------------------------------
smmimport
sat(2).ID = 'SMM/ACRIM1';
sat(2).JD = juliandate(datetime(SMM.year,SMM.month,SMM.day));
sat(2).TSI = SMM.TSI;
sat(2).version = '1';
%--------------------------------------------------
% ERBE/ERBS
%--------------------------------------------------
erbsimport
sat(3).ID = 'ERBS/ERBE';
sat(3).JD = ERBS.JD;
sat(3).TSI = ERBS.TSI;
%--------------------------------------------------
% ACRIM2
%--------------------------------------------------
acrim2import
sat(4).ID = 'UARS/ACRIM2';
dates = ACRIM2.YEAR;
[years, days] = ConvertSerialYearToDate(dates);
days(days<1)=1; %Eliminate errors from rounding down to 0
dateStr = strcat(string(num2str(years)), repmat("-",length(years),1),...
    string(num2str(floor(days))));
sat(4).JD = juliandate(datetime(dateStr,'InputFormat', 'uuuu-D'));
sat(4).TSI = ACRIM2.TSI;

%--------------------------------------------------
% VIRGO/SOHO
%--------------------------------------------------
virgoimport2
sat(5).ID = 'SOHO/VIRGO';
sat(5).JD = VIRGO.JulianDate;
sat(5).TSI = VIRGO.VIRGOAB_correc; %Note there exist other TSI products from this sat
sat(5).version = 'V8_20220105';
%--------------------------------------------------
% ACRIM3/ACRIMSAT
%--------------------------------------------------
acrim3import
sat(6).ID = 'ACRIMSAT/ACRIM3';
sat(6).JD = juliandate(datetime(ACRIM3.YYMMDD,'InputFormat', 'yyMMdd'));
sat(6).TSI = ACRIM3.TSI;
sat(6).version= 'v1L2';
%--------------------------------------------------
% SORCE
%--------------------------------------------------
sorceimport
sat(7).ID = 'SORCE/TIM';
sat(7).JD = juliandate(datetime(int32(SORCE.YYYYMMDD), 'ConvertFrom', 'YYYYMMDD'));
sat(7).TSI = SORCE.TSI;
nanInd = find(sat(7).TSI < 1000);
sat(7).TSI(nanInd) = NaN; %set zeros to NaN
sat(7).version = 'v19L3';

%--------------------------------------------------
% PREMOS
%--------------------------------------------------
load('PREMOSdaily.mat')
sat(8).ID = 'PICARD/PREMOS';
sat(8).JD = PREMOSDAY.JD';
sat(8).TSI = PREMOSDAY.TSI'; %No NaNs recorded
sat(8).version = 'v2';

%--------------------------------------------------
% TCTE
%--------------------------------------------------
tcteimport
sat(9).ID = 'TCTE/TIM';
sat(9).JD = juliandate(datetime(int32(TCTE.YYYYMMDD), 'ConvertFrom', 'YYYYMMDD'));
sat(9).TSI = TCTE.TSI;
nanInd = find(sat(9).TSI < 1000);
sat(9).TSI(nanInd) = NaN; %set zeros to NaN
sat(9).version = 'v3L3';
%--------------------------------------------------
% TSIS
%--------------------------------------------------
tsisimport
sat(10).ID = 'TSIS-1';
sat(10).JD = juliandate(datetime(int32(TSIS.YYYYMMDD), 'ConvertFrom', 'YYYYMMDD'));
sat(10).TSI = TSIS.TSI;
nanInd = find(sat(10).TSI < 1000);
sat(10).TSI(nanInd) = NaN; %set zeros to NaN
sat(10).version = 'v3L3';
%--------------------------------------------------
% HF/NIMBUS PMOD CORRECTED
%--------------------------------------------------
nimbuscorrectedimport
sat(11).ID = 'NIMBUS-7/HF PMOD CORRECTED';
sat(11).JD = NIMBUSCORR.JD;
sat(11).TSI = NIMBUSCORR.TSI;
sat(11).TSI(sat(11).TSI < 0)=NaN;
sat(11).version = 'ftp.pmodwrc.ch/pub/data/irradiance/acrim/AdditionalFiles accessed 22/07/07';
%--------------------------------------------------
% ACRIM1 PMOD CORRECTED
%--------------------------------------------------
acrim1correctedimport
sat(12).ID = 'SMM/ACRIM1 PMOD CORRECTED';
sat(12).JD = ACRIM1CORR.JD;
sat(12).TSI = ACRIM1CORR.TSI;
sat(12).TSI(sat(12).TSI < 0)=NaN;
sat(12).version = 'ftp.pmodwrc.ch/pub/data/irradiance/acrim/AdditionalFiles accessed 22/07/07';
%--------------------------------------------------
% ACRIM2 PMOD CORRECTED
%--------------------------------------------------
acrim2correctedimport
sat(13).ID = 'UARS/ACRIM2 PMOD CORRECTED';
sat(13).JD = ACRIM2CORR.JD;
sat(13).TSI = ACRIM2CORR.TSI;
sat(13).TSI(sat(13).TSI < 0)=NaN;
sat(13).version = 'ftp.pmodwrc.ch/pub/data/irradiance/acrim/AdditionalFiles accessed 22/07/07';
%--------------------------------------------------
% ERBS PMOD CORRECTED
%--------------------------------------------------
erbscorrectedimport
sat(14).ID = 'ERBS/ERBE PMOD CORRECTED';
sat(14).JD = ERBSCORR.JD;
sat(14).TSI = ERBSCORR.TSI;
sat(14).TSI(sat(14).TSI < 0)=NaN;
sat(14).version = 'ftp.pmodwrc.ch/pub/data/irradiance/acrim/AdditionalFiles accessed 22/07/07';

%save('../data/codegenerated/satUpdated.mat', 'sat')

function [year, day] = ConvertSerialYearToDate( y )
  %Create the makings of a date object using a decimal year
  year = floor(y);
  partialYear = mod(y,1);
  date0 = datenum(num2str(year),'yyyy');
  date1 = datenum(num2str(year+1),'yyyy');
  daysInYear = date1 - date0;
  day =  partialYear .* daysInYear;
end