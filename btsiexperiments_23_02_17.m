%BTSIEXPERIMENTS_23_02_17 Set-up and run the experiments that are cited in
%the ACRIM-Gap paper
%
% Ted Amdur
% 2023/02/17

clearvars
testBed=0; %Test different configurations of BTSI
mainExperiment=1; %Run full experiment that results in published results and figures
noEarlyACRIM=0; %Remove the first 20 months of observations from ACRIM1
noSpots=0; %Run BTSI experiment without sunspot record
noERBE=0; %Run without ERBE
dataPMOD=0; %Run with PMOD corrected data

if testBed
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=true;
    opts.normalize=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.magDependent=false; %Start by seeing if I can replicate standard BTSI
    opts.saveFile='ar2_23_03_17.mat';
    runchain_23_02_21([],[],[],[],opts);
end
if mainExperiment
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=true;
    opts.normalize=true;
    opts.magDependent=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='ar2main_23_03_17.mat';
    runchain_23_02_21([],[],[],[],opts);
    
end
if noEarlyACRIM
    load obs_23_02_01.mat
    valM(1:20,1)=NaN;oM(1:20,1)=false; %Turn off first 20 months of ACRIM1
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=false;
    opts.normalize=true;
    opts.magDependent=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='ar2_noearlyacrim_23_03_17.mat';
    runchain_23_02_21(valM,oM,dateM,colLabels,opts);
end
if noSpots
    load obs_23_02_01.mat
    valM(:,6)=NaN;oM(:,6)=false;%Turn off first 20 months of ACRIM1
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=true;
    opts.normalize=true;
    opts.magDependent=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='ar2_nospots_23_03_17.mat';
    runchain_23_02_21(valM,oM,dateM,colLabels,opts);
end
if noERBE
    load obs_23_02_01.mat
    valM(:,4)=NaN;oM(:,4)=false; %Turn off ERBE 
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=true;
    opts.normalize=true;
    opts.magDependent=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='ar2_noserbe_23_03_17.mat';
    runchain_23_02_21(valM,oM,dateM,colLabels,opts);
end
if dataPMOD
    load obsPMOD_23_02_01.mat
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.randomizeChain=false;
    opts.logContributions=true;
    opts.normalize=true;
    opts.magDependent=true;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor, 1 default
    opts.saveFile='ar2_pmod_23_03_17.mat';
    runchain_23_02_21(valM,oM,dateM,colLabels,opts);
end
