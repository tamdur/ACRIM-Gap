%Make synthetic datasets for ACRIM and PMOD scenarios
%
% Ted Amdur
% 2023/01/25

nSets=1000;
%Create starting parameters
[Ainit,epsilon,rho,t,oM] = initobservationmodelparams;
%Create ACRIM and PMOD datasets
ACRIM=struct;
PMOD=struct;
for ii=1:nSets
    ACRIM(ii).valM=gensynthobs('ACRIM',Ainit,epsilon,rho,t,oM);
    [PMOD(ii).valM,dateM]=gensynthobs('PMOD',Ainit,epsilon,rho,t,oM);
end
save('2scenario_23_01_25.mat','ACRIM','PMOD','dateM','Ainit','epsilon','rho','t');