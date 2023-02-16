%RUNALTS run alterantive formulations of BTSI to test robustness of model
% Ted Amdur
% 2023/02/15
clearvars
lag1=0;
lag3=0;
synthGeneric=0;
synthAltH=0;
synthNoRho=0;
largeSynth=1;
noERBE=0;
altACRIM1=0;

%Make sure everything is visible
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/tools')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/modeleval')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/modeleval/alt_runs')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/modeleval/three_scenario')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/tools/BoE')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/mat_files')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/chain_output')

if lag1
    load obs_23_02_01.mat
    opts.burnin = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=1;
    opts.saveFile='ar1_23_02_15.mat';
    runchain_23_01_13(valM,oM,colLabels,opts);
end
if lag3
    load obs_23_02_01.mat
    opts.burnin = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=3;
    opts.saveFile='ar3_23_02_15.mat';
    runchain_23_01_13(valM,oM,colLabels,opts);
end
if synthGeneric
    [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(1,[],[]);
    runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],'threetestcluster_generic_23_02_15b.mat');
end
if synthAltH
    [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(2,[],'3scenario_rng2__23_02_15.mat');
    runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],'threetestcluster_23_02_15.mat');
end
if synthNoRho
    [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(1,0,'3scenario_23_02_15.mat');
    runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],'threetestcluster_norho_23_02_15.mat');
end
if largeSynth
    ii=8;
    savePth=['threetestcluster_rng' num2str(ii) '_23_02_16.mat'];
    [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(ii,0,[]);
    runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],savePth);
end
        
if noERBE
    load obs_23_02_01.mat
    valM(:,4)=NaN;oM(:,4)=false; %Turn off ERBE
    opts.burnin = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.logContributions=false;
    opts.saveFile='ar2_noERBE_23_02_15.mat';
    runchain_23_01_13(valM,oM,colLabels,opts);
end
if altACRIM1
    load obs_23_02_01.mat
    valM(1:20,1)=NaN;oM(1:20,1)=false; %Turn off first 20 months of ACRIM1
    opts.burnin = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.logContributions=false;
    opts.saveFile='ar2_shortACRIM1_23_02_15.mat';
    runchain_23_01_13(valM,oM,colLabels,opts);
end
