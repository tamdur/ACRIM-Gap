function [Ainit,epsilon,rho,t,oM] = initobservationmodelparams
%INITOBSERVATIONMODELPARMS 
%OUTPUT:
%       Ainit: observation model truth (randomly generated)
%       epsilon: true observational noise from sigY BTSI fit
%       rho: True AR(1) parameter for satellites
%       t: array of observer ages in decades
%       

%Load the actual observers
obsmatrix='obs_23_02_01'; %Load data array, with colLabels corresponding to observer source for each column
load(obsmatrix);
T=size(valM,1);
N=size(valM,2);
oindex=[1 1 1 1 1 1 0]; %oindex=1 for observers with varying offset, 0 for fixed
tindex=[1 1 0 1 1 0 1]; %tindex=1 for observers with time dependent drift, 0 otherwise
sindex=[0 0 1 0 0 1 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
satindex=find(tindex); %satindex=1 for observers that are satellites, 0 otherwise

%Load the BTSI output
load ar2_23_01_14.mat

A0=mean(A,3);
Ainit=zeros(size(A0));
Ainit(:,2)=A0(:,2); Ainit(3,2)=0.01; Ainit(6,2)=150; %Create scaling for Mg-II, spots
Ainit(:,1)=0.5.*randn(7,1); %Create random offsets
Ainit(:,1)=Ainit(:,1).*Ainit(:,2); %Use scaling on offsets
Ainit(satindex,3)=A0(satindex(randperm(length(satindex))),3); %Randomly permute satellite drifts

%Calculate the average noise estimates for each observer, use that to
%create observational errors for those observers
epsilon=sqrt(mean(sigY,2));

%rho=0.9; %Assumed AR(1) parameter for satellites, can infer from obs in future editions

t =repmat(linspace(0,T./120,T)',[1 N]); %Make time rows for t
for ii=1:N
    TM=mean(t(oM(:,ii),ii));
    if isnan(TM)
        TM=0;
    end
    t(:,ii)=t(:,ii)-TM;
end

thresh=48; %Threshold for using overlapping observers be it years or months
%Generate prior hyperparameter estimates 
nObs = length(colLabels);
%Make labels for which records fall under satellite vs proxy
isProx=~satindex;
proxInd=find(isProx);
satInd=satindex;
T=size(valM,1);
t1 = (0:T)';
ind=1;
for ii=1:nObs
    if any(ii==satInd) %Only make comparisons for primary observers
        vSat = NaN(nObs,1);
        iS = ii+1;
        while iS <= nObs
            if any(iS==satInd) && sum(logical(oM(:,iS).*oM(:,ii)))>=thresh
                overlap = logical(oM(:,iS).*oM(:,ii));
                %Make a set of priors for the proxy observations using all satellites
                satdiff=valM(overlap,ii)-valM(overlap,iS);
                %Revised 9/8/21 to be in native units
                pred=[ones(sum(overlap),1) t1(overlap)];
                b = regress(satdiff,pred);
                residual(ind).vals=satdiff-pred*b;
                residual(ind).sat1=colLabels(iS);
                residual(ind).sat2=colLabels(ii);
                ind=ind+1;
            end
            iS = iS + 1;
        end
    end
end
%Estimate autocorrelation for each
for ii=1:length(residual)
    x=residual(ii).vals;
    rho(ii)=x(2:end)'*x(1:end-1)/(x(1:end-1)'*x(1:end-1));
end
rho=mean(rho); %use the average autocorrelation parameter seen amongst models

end