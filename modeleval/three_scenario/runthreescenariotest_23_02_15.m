function runthreescenariotest_23_02_15(ACRIM,PMOD,AP,setInfo,rngN,savePath,opts)

%Run three scenario synthetic data test in parallel computing environment
% Ted Amdur
% 2023/02/15

parpool('local',str2num(getenv('SLURM_CPUS_PER_TASK')))
%Make sure everything is visible
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/tools')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/tools/gibbs_functions')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/modeleval')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/modeleval/alt_runs')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/modeleval/three_scenario')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/tools/BoE')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/mat_files')

if ~exist('ACRIM','var') || ~exist('PMOD','var') || ~exist('AP','var') ||...
        isempty(ACRIM) || isempty(PMOD) || isempty(AP)
    %Load the synthetic datasets to be examined
    load 2scenario_23_01_31_PMODProxyb.mat
    AP=ACRIM;
    load 2scenario_23_01_31.mat
end

if ~exist('rngN','var') || isempty(rngN)
    rng(1)
else
    rng(rngN)
end

if ~exist('opts','var') || isempty(opts)
    opts.burnin = 500; %Number of burn-in reps assumed for chain length analysis
    opts.reps=1500; %Total length of chain, including burn-in
    opts.dispProgress=false;
    opts.normalize=true;
    opts.magDependent=false;
end


oM=setInfo.oM;
colLabels=setInfo.colLabels;
dateM=setInfo.dateM;
tN=length(ACRIM); %Number of synthetic datasets to be inferred
parfor ii=1:tN
    %Run BTSI on ACRIM scenario
    tic;
    [xAll,sigY,~,~,~,A,~,~] = runchain_23_02_21(ACRIM(ii).valM,oM,dateM,colLabels,opts);
    [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap]=returnscenarioinfo(xAll,sigY,A,dateM);
    threeTest(ii).ACRIM.Aout=Aout;
    threeTest(ii).ACRIM.sigYOut=sigYOut;
    threeTest(ii).ACRIM.AUnc=AUnc;
    threeTest(ii).ACRIM.sigYUnc=sigYUnc;
    threeTest(ii).ACRIM.muGap=muGap;
    threeTest(ii).ACRIM.uncGap=uncGap;
    threeTest(ii).ACRIM.xm=mean(xAll,2);
    threeTest(ii).ACRIM.tRun=toc;
    
    %Run BTSI on PMOD scenario
    tic;
    [xAll,sigY,~,~,~,A,~,~] = runchain_23_02_21(PMOD(ii).valM,oM,dateM,colLabels,opts);
    [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap]=returnscenarioinfo(xAll,sigY,A,dateM);
    threeTest(ii).PMOD.Aout=Aout;
    threeTest(ii).PMOD.sigYOut=sigYOut;
    threeTest(ii).PMOD.AUnc=AUnc;
    threeTest(ii).PMOD.sigYUnc=sigYUnc;
    threeTest(ii).PMOD.muGap=muGap;
    threeTest(ii).PMOD.uncGap=uncGap;
    threeTest(ii).PMOD.xm=mean(xAll,2);
    threeTest(ii).PMOD.tRun=toc;
    
    
    %Run BTSI on A/P scenario
    tic;
    [xAll,sigY,~,~,~,A,~,~] = runchain_23_02_21(AP(ii).valM,oM,dateM,colLabels,opts);
    [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap]=returnscenarioinfo(xAll,sigY,A,dateM);
    threeTest(ii).AP.Aout=Aout;
    threeTest(ii).AP.sigYOut=sigYOut;
    threeTest(ii).AP.AUnc=AUnc;
    threeTest(ii).AP.sigYUnc=sigYUnc;
    threeTest(ii).AP.muGap=muGap;
    threeTest(ii).AP.uncGap=uncGap;
    threeTest(ii).AP.xm=mean(xAll,2);
    threeTest(ii).AP.tRun=toc;
end
scriptName=mfilename;
save(savePath,'threeTest','scriptName')
end

function [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap]=returnscenarioinfo(xAll,sigY,A,dateM)
%Return observation model values
Aout=mean(A,3);
sigYOut=mean(sigY,2);
AUnc=std(A,0,3);
sigYUnc=std(sigY,0,2);

%Return muGap and uncGap
%Get the variables of interest from this single scenario
dateS=getdates;
smoothWindow=1; %smoothing of data (if value different than 1)
tsix = prctile(xAll',50)';
tsiAll=xAll;
tsix(:) = smoothPH(tsix(:),smoothWindow);
xm=tsix;

%Create time windows for year before and year after acrim gap
stInt=datejd([dateS.acrim(1)-365 dateS.acrim(1)]);
endInt=datejd([dateS.acrim(end) dateS.acrim(end)+365]);

%Get TSI at the beginning and end interval for each
stBTSI=dateM >= stInt(1)& dateM <=stInt(end);
endBTSI=dateM >=  endInt(1) & dateM<=endInt(end);
muGap=[mean(xm(stBTSI));
    mean(xm(endBTSI))];
muGap=diff(muGap);

%Get max and min change
BTSIAll=zeros(size(tsiAll,2),3);
for ii=1:size(tsiAll,2)
    BTSIAll(ii,1:2)=[mean(tsiAll(stBTSI,ii)),mean(tsiAll(endBTSI,ii))];
    BTSIAll(ii,3)=diff(BTSIAll(ii,1:2));
end
%Use z-score table for normal distribution to estimate uncertainty
uncGap=(quantile(BTSIAll(:,3),0.975)-quantile(BTSIAll(:,3),0.025))./(2.*1.96);
end


