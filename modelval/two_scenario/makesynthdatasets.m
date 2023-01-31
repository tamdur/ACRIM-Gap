%Make synthetic datasets for ACRIM and PMOD scenarios
%
% Ted Amdur
% 2023/01/25
clearvars
rng(1); %Create common random seed for debugging/validation purposes
nSets=1000;
%Create starting parameters
[Ainit,epsilon,rho,t,oM] = initobservationmodelparams;
%Create ACRIM and PMOD datasets
ACRIM=struct;
PMOD=struct;
for ii=1:nSets
    ACRIM(ii).valM=gensynthobs('ACRIM/PMOD proxy',Ainit,epsilon,rho,t,oM);
    [PMOD(ii).valM,dateM]=gensynthobs('PMOD',Ainit,epsilon,rho,t,oM);
end
%Get other metadata
%Load data array, with colLabels corresponding to observer source for each column
%This is the obsmatrix used in initobservationmodelparams.m
obsmatrix='obs_23_01_13'; 
load(obsmatrix);
save('2scenario_23_01_31_PMODproxy.mat','ACRIM','PMOD','dateM','Ainit','epsilon','rho','t','oM','colLabels');