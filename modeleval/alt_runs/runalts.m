%RUNALTS run alterantive formulations of BTSI to test robustness of model
% Ted Amdur
% 2023/02/15
clearvars
lag1=0;
lag3=0;
synthGeneric=1;
synthAltH=0;
synthNoRho=0;
largeSynth=0;
noERBE=0;
altACRIM1=0;

%Make sure everything is visible
parpool('local',str2num(getenv('SLURM_CPUS_PER_TASK')))
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
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=1;
    opts.NRLTSIprior=true;
    opts.saveFile='ar1_23_02_21.mat';
    runchain_23_02_21(valM,oM,colLabels,opts);
end
if lag3
    load obs_23_02_01.mat
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=3;
    opts.NRLTSIprior=true;
    opts.saveFile='ar3_23_02_15.mat';
    runchain_23_02_21(valM,oM,colLabels,opts);
end
if synthGeneric
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor
    opts.NRLTSIprior=true;
    opts.normalize=true;
    opts.magDependent=true;
    [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(1,[],[]);
    runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],'threetestcluster_generic_23_03_17.mat',opts);
end
if synthAltH
    [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(2,[],'3scenario_rng2__23_02_15.mat');
    runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],'threetestcluster_23_02_15.mat');
end
if synthNoRho
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor
    opts.NRLTSIprior=true;
    opts.normalize=true;
    opts.magDependent=true;
    [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(1,0,[]);
    runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],'threetestcluster_norho_23_03_17.mat',opts);
end
if largeSynth
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=false;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor
    opts.NRLTSIprior=true;
    for ii=1:10
        savePth=['threetestcluster_rng' num2str(ii) '_23_03_17.mat'];
        [AP,ACRIM,PMOD,setInfo]=makesynthdatasets(ii,[],[]);
        runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,[],savePth,opts);
    end
end
        
if noERBE
    load obs_23_02_01.mat
    valM(:,4)=NaN;oM(:,4)=false; %Turn off ERBE
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.logContributions=false;
    opts.NRLTSIprior=true;
    opts.saveFile='ar2_noERBE_23_02_21.mat';
    runchain_23_02_21(valM,oM,colLabels,opts);
end
if altACRIM1
    load obs_23_02_01.mat
    valM(1:20,1)=NaN;oM(1:20,1)=false; %Turn off first 20 months of ACRIM1
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.logContributions=false;
    opts.NRLTSIprior=true;
    opts.saveFile='ar2_shortACRIM1_23_02_21.mat';
    runchain_23_02_21(valM,oM,colLabels,opts);
end
