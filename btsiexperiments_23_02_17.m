%BTSIEXPERIMENTS_23_02_17 Set-up and run the experiments that are cited in
%the ACRIM-Gap paper
%
% Ted Amdur
% 2023/02/17

mainExperiment=1; %Run full experiment that results in published results and figures


if mainExperiment
    opts.burn = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=10500; %Total length of chain, including burn-in
    opts.dispProgress=true;
    opts.lags=2;
    opts.NRLTSIprior=true;
    opts.saveFile='ar2_23_02_21b.mat';
    opts.randomizeChain=false;
    opts.excludeFliers=false;
    opts.logContributions=true;
    opts.normalize=false;
    opts.HsigScale=1; %Change the variance parameters of Hsig by scaling factor
                      %1 to maintain default
    runchain_23_02_21([],[],[],opts);
end