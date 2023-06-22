function tsi = maketsiseries
%MAKETSISERIES Return monthly time series for the scenarios scenarios being
%examined: the ACRIM composite, the PMOD composite, and the Community-Consensus 
%model. These will represent the underlying 'true' TSI of the system.
% Ted Amdur
% 2023/06/21
%OUTPUT:
%       tsi: a structure array containing vectors ACRIM: a Nx1 vector of
%       the ACRIM composite TSI observations, PMOD: a Nx1 vector of the
%       PMOD CPMDF observations, SOLID: a Nx1 vector of the
%       SOLID observations and date: an Nx1 datetime vector
%       containing the dates corresponding to the ACRIM PMOD, and SOLID
%       values. Note, here we use SOLID corrected TSI values


%Load external TSI values
load oTSI_23_02_01.mat %From readothertsi.m in the code_22_06 directory

%Load ACRIM values and datetime vector
ACRIM=oTSI(6).mthtsi;
dateACRIM=oTSI(6).mthdatetime;

%Load PMOD values and datetime vector
PMOD=oTSI(7).mthtsi;
datePMOD=oTSI(7).mthdatetime;

%Load SOLID values and datetime vector
SOLID=oTSI(9).mthtsi;
dateSOLID=oTSI(9).mthdatetime;

dateS=getdates;
dates=datejd(dateS.acrimplusfive);


%Find the indices of the start and end date in the datetime vectors
indACRIM = find(dateACRIM >= dates(1) & dateACRIM <= dates(2));
indPMOD = find(datePMOD >= dates(1) & datePMOD <= dates(2));
indSOLID = find(dateSOLID >= dates(1) & dateSOLID <= dates(2));

%Extract the relevant data
tsi.ACRIM =ACRIM(indACRIM);
tsi.PMOD =PMOD(indPMOD);
tsi.SOLID=SOLID(indSOLID);
tsi.date = dateACRIM(indACRIM);
end

