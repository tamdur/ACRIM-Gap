function [valS,dateM] = gensynthobs_23_03_10(scenario,Ainit,epsilon,rho,t,oM)
%GENSYNTHOBS Generate a set of synthetic observations for ACRIM-Gap BTSI
%testing.
% INPUTS:
%           scenario: char vector, either 'ACRIM' or 'PMOD', that informs
%           which baseline TSI scenario to use. Error otherwise
%           Ainit: assumed true observation matrix, from
%           initobservationmodelparams.m
%           epsilon: assumed instrument noise, from initobservationmodelparams.m
%           rho: assumed autocorrelation parameter for satellite noise,
%           from initobservationmodelparams.m
%           t: array of instrumental ages, from initobservationmodelparams.m
%           oM: array of present and missing obs, from initobservationmodelparams.m
% OUTPUTS:
%           valS: valM observation array with synthetic data generated
%           using input assumptions. Proxies include offset and synthetic Gaussian
%           noise, and satellites include offset, linear drift, and AR(1)
%           autocorrelated residual error


%Load TSI series used for the assumed truth
tsi = twotsiseries;

if strcmp(scenario,'ACRIM')
    x=tsi.ACRIM;
elseif strcmp(scenario,'ACRIM/PMOD proxy')
    x=tsi.ACRIM;
    xprox=tsi.PMOD;
else 
    x=tsi.PMOD;
end
dateM=tsi.date;

%Create array to hold synthetic outputs
T=length(dateM);
N=size(Ainit,1);
valS=zeros(T,N);

%Determine which columns correspond to synthetic proxies
prox=round(Ainit(:,2))~=1;
proxInd=find(prox);satInd=find(~prox);

%First, create synthetic data for proxies (offset and Gaussian noise)
for ii=1:length(proxInd)
    if strcmp(scenario,'ACRIM/PMOD proxy')
        valS(:,proxInd(ii))=(xprox-mean(xprox)).*Ainit(proxInd(ii),2)+Ainit(proxInd(ii),1);
    else
        valS(:,proxInd(ii))=(x-mean(x)).*Ainit(proxInd(ii),2)+Ainit(proxInd(ii),1);
    end
    valS(:,proxInd(ii))=valS(:,proxInd(ii))+epsilon(proxInd(ii)).*randn(T,1);
end

%Next, create synthetic data for satellites (offset and linear drift and
%AR(1) error
for ii=1:length(satInd)
    valS(:,satInd(ii))=x.*Ainit(satInd(ii),2)+Ainit(satInd(ii),1) + Ainit(satInd(ii),3).*t(:,satInd(ii));
    error=zeros(T,1);
    error(1)=epsilon(satInd(ii)).*randn;
    sigma=sqrt(epsilon(satInd(ii)).^2.*(1-rho.^2));
    for iT=2:T
        %Get effective innovation from https://online.stat.psu.edu/stat501/book/export/html/995
        error(iT)=rho.*error(iT-1)+sigma.*randn;
    end
    valS(:,satInd(ii))=valS(:,satInd(ii))+error;
    valS(:,satInd(ii))=valS(:,satInd(ii))-nanmean(x); %Center
end

%Last, make observation matrix contain the same missing obs structure as
%original
valS(~oM)=NaN;
    

end

