function [Ainit,epsilon,rho] = initobservationmodelparams
%INITOBSERVATIONMODELPARMS 
%OUTPUT:
%       tsi: a structure array containing vectors ACRIM: a Nx1 vector of
%       the ACRIM composite TSI observations, PMOD: a Nx1 vector of the
%       PMOD CPMDF observations, and date: an Nx1 datetime vector
%       containing the dates corresponding to the ACRIM and PMOD values

rng(1); %Use a consistent random seed for reproducibility

%Load the actual observers
obsmatrix='obs_23_01_13'; %Load data array, with colLabels corresponding to observer source for each column
load(obsmatrix);
oindex=[1 1 1 1 1 1 0]; %oindex=1 for observers with varying offset, 0 for fixed
tindex=[1 1 0 1 1 0 1]; %tindex=1 for observers with time dependent drift, 0 otherwise
sindex=[0 0 1 0 0 1 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
satindex=find(tindex); %satindex=1 for observers that are satellites, 0 otherwise

%Load the BTSI output
load ar2_23_01_14.mat

A0=mean(A,3);
A1=zeros(size(A0));
A1(:,2)=A0(:,2); A1(3,2)=0.01; A1(6,2)=150; %Create scaling for Mg-II, spots
A1(:,1)=0.5.*randn(7,1); %Create random offsets
A1(:,1)=A1(:,1).*A1(:,2); %Use scaling on offsets
A1(satindex,3)=A0(satindex(randperm(length(satindex))),3); %Randomly permute satellite drifts

%Calculate the average noise estimates for each observer, use that to
%create observational errors for those observers
e0=mean(sigY,2);

rho=0.9; %Assumed AR(1) parameter for satellites, can infer from obs in future editions



% thresh=48; %Threshold for using overlapping observers be it years or months
% %Generate prior hyperparameter estimates 
% nObs = length(colLabels);
% %Make labels for which records fall under satellite vs proxy
% isProx=~satindex;
% proxInd=find(isProx);
% satInd=find(satindex);
% T=size(valM,1);
% t = (0:T)';

% for ii=1:nObs
%     if any(ii==satInd) %Only make comparisons for primary observers
%         vSat = NaN(nObs,1);
%         iS = 1;
%         while iS <= nObs
%             if iS ~= ii && satindex(iS)
%                 overlap = logical(oM(:,iS).*oM(:,ii));
%                 %Make a set of priors for the proxy observations using all satellites
%                 satdiff=valM(overlap,ii)-valM(overlap,iS);
%                 %Revised 9/8/21 to be in native units
%                 pred=[ones(sum(overlap),1) t(overlap)];
%                 b = regress(satdiff,pred);
%                 residual=satdiff-b*pred;
%             end
%             iS = iS + 1;
%         end
%         %Include the data collected into a structure object
%         if isfield(opts,'type')
%             obsPrior(ii).type = opts.type;
%         else
%             obsPrior(ii).type = "sat";
%         end
%         obsPrior(ii).name = colLabels(ii);
%         obsPrior(ii).std = sqrt(nanmean(vSat))./sqrt(2); %Correct for satellite noise coming from two observers
%     end
% end

%Load the 


end