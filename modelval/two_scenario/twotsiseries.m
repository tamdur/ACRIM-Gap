function [tsi] = twotsiseries
%TWOTSISERIES Return monthly time series for the two scenarios being
%examined: the ACRIM composite and the PMOD composite. These will represent
%the underlying 'true' TSI of the system.
% Ted Amdur
% 2023/01/24
%OUTPUT:
%       tsi: a structure array containing vectors ACRIM: a Nx1 vector of
%       the ACRIM composite TSI observations, PMOD: a Nx1 vector of the
%       PMOD CPMDF observations, and date: an Nx1 datetime vector
%       containing the dates corresponding to the ACRIM and PMOD values


%Load external TSI values
load oTSI_22_10_27.mat %From readothertsi.m in the code_22_06 directory

%Load ACRIM values and datetime vector
ACRIM=oTSI(6).mthtsi;
dateACRIM=oTSI(6).mthdatetime;

%Load PMOD values and datetime vector
PMOD=oTSI(7).mthtsi;
datePMOD=oTSI(7).mthdatetime;


dateS=getdates;
dates=datejd(dateS.acrimplusfive);


%Find the indices of the start and end date in the datetime vectors
indACRIM = find(dateACRIM >= dates(1) & dateACRIM <= dates(2));
indPMOD = find(datePMOD >= dates(1) & datePMOD <= dates(2));

%Extract the relevant data
tsi.ACRIM = ACRIM(indACRIM);
tsi.PMOD = PMOD(indPMOD);
tsi.date = dateACRIM(indACRIM);

end

