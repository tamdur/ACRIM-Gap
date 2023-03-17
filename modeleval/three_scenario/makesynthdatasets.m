function [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(rngN,rhoInput,savePath)
%Make synthetic datasets for ACRIM and PMOD scenarios
%
% Ted Amdur
% 2023/01/25
if ~exist('rngN','var') || isempty(rngN)
    rng(1)
else
    rng(rngN)
end
nSets=1000;
%Create starting parameters
[Ainit,epsilon,rho,t,oM] = initobservationmodelparams;

if exist('rhoInput','var') && ~isempty(rhoInput)
    %Change rho if desired
    rho=rhoInput;
end
%Create ACRIM and PMOD datasets
AP=struct;
ACRIM=struct;
PMOD=struct;
for ii=1:nSets
    ACRIM(ii).valM=gensynthobs('ACRIM',Ainit,epsilon,rho,t,oM);
    AP(ii).valM=gensynthobs('ACRIM/PMOD proxy',Ainit,epsilon,rho,t,oM);
    [PMOD(ii).valM,dateM]=gensynthobs('PMOD',Ainit,epsilon,rho,t,oM);
end
%Get other metadata
%Load data array, with colLabels corresponding to observer source for each column
%This is the obsmatrix used in initobservationmodelparams.m
obsmatrix='obs_23_02_01'; 
load(obsmatrix);
setInfo.dateM=dateM;
setInfo.Ainit=Ainit;
setInfo.epsilon=epsilon;
setInfo.rho=rho;
setInfo.t=t;
setInfo.oM=oM;
setInfo.colLabels=colLabels;
setInfo.obsmatrix=obsmatrix;
setInfo.notes='Fixes prior error where satellite observations were not centered to zero';
if exist('savePath','var') && ~isempty(savePath)
    save(savePath,'ACRIM','PMOD','AP','setInfo');
end
end