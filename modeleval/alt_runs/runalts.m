%RUNALTS run alterantive formulations of BTSI to test robustness of model
% Ted Amdur
% 2023/02/15

lag1=0;
lag3=0;
synthAltH=1;
noERBE=0;
altACRIM1=0;

if lag1
    load obs_23_02_01.mat
    opts.burnin = 500; %Number of burn-in reps assumed for chain length analysis
    opts.chain=10500; %Total length of chain, including burn-in
    opts.randomizeChain=true; %Disperse chains for independent comparison
    opts.dispProgress=true;
    opts.lags=1;
    opts.saveFile='ar1_23_02_15.mat';
    runchain_23_01_13(valM,oM,colLabels,opts);
end
if lag3
    load obs_23_02_01.mat
    opts.burnin = 500; %Number of burn-in reps assumed for chain length analysis
    opts.chain=10500; %Total length of chain, including burn-in
    opts.randomizeChain=true; %Disperse chains for independent comparison
    opts.dispProgress=true;
    opts.lags=3;
    opts.saveFile='ar3_23_02_15.mat';
    runchain_23_01_13(valM,oM,colLabels,opts);
end
if synthAltH
    [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(2);
    runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],'threetestcluster_23_02_15.mat');
end
if noERBE
end
if altACRIM1
end