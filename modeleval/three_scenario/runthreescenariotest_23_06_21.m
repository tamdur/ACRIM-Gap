function runthreescenariotest_23_06_21(sc1,sc2,sc3,scenarios,setInfo,rngN,savePath,opts)

%Run three scenario synthetic data test in parallel computing environment
% Ted Amdur
% 2023/02/15
%
% sc1 is the first scenario, where TSI follows ACRIM, sc2 is the second
% scenario, where TSI follows either SOLID or PMOD, and sc3 is the hybrid
% scenario where proxies reflect PMOD or SOLID and satellites reflect
% ACRIM.

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

if ~exist('sc1','var') || ~exist('sc2','var') || ~exist('sc3','var') ||...
        isempty(sc1) || isempty(sc2) || isempty(sc3)
    %Load the synthetic datasets to be examined
    load 2scenario_23_01_31_PMODProxyb.mat
    sc3=ACRIM;
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
    opts.magDependent=true;
end


oM=setInfo.oM;
colLabels=setInfo.colLabels;
dateM=setInfo.dateM;
tN=length(sc1); %Number of synthetic datasets to be inferred
parfor ii=1:tN
    %Run BTSI on ACRIM scenario
    tic;
    [xAll,sigY,~,~,~,A,~,outDat] = runchain_23_02_21(sc1(ii).valM,oM,dateM,colLabels,opts);
    [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap,H0]=returnscenarioinfo(xAll,sigY,A,dateM,outDat);
    threeTest(ii).sc1.Aout=Aout;
    threeTest(ii).sc1.sigYOut=sigYOut;
    threeTest(ii).sc1.AUnc=AUnc;
    threeTest(ii).sc1.sigYUnc=sigYUnc;
    threeTest(ii).sc1.muGap=muGap;
    threeTest(ii).sc1.uncGap=uncGap;
    threeTest(ii).sc1.H0=H0;
    threeTest(ii).sc1.xm=mean(xAll,2);
    threeTest(ii).sc1.tRun=toc;
    
    %Run BTSI on PMOD or SOLID scenario
    tic;
    [xAll,sigY,~,~,~,A,~,outDat] = runchain_23_02_21(sc2(ii).valM,oM,dateM,colLabels,opts);
    [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap,H0]=returnscenarioinfo(xAll,sigY,A,dateM,outDat);
    threeTest(ii).sc2.Aout=Aout;
    threeTest(ii).sc2.sigYOut=sigYOut;
    threeTest(ii).sc2.AUnc=AUnc;
    threeTest(ii).sc2.sigYUnc=sigYUnc;
    threeTest(ii).sc2.muGap=muGap;
    threeTest(ii).sc2.uncGap=uncGap;
    threeTest(ii).sc2.H0=H0;
    threeTest(ii).sc2.xm=mean(xAll,2);
    threeTest(ii).sc2.tRun=toc;
    
    
    %Run BTSI on A/P scenario
    tic;
    [xAll,sigY,~,~,~,A,~,outDat] = runchain_23_02_21(sc3(ii).valM,oM,dateM,colLabels,opts);
    [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap,H0]=returnscenarioinfo(xAll,sigY,A,dateM,outDat);
    threeTest(ii).sc3.Aout=Aout;
    threeTest(ii).sc3.sigYOut=sigYOut;
    threeTest(ii).sc3.AUnc=AUnc;
    threeTest(ii).sc3.sigYUnc=sigYUnc;
    threeTest(ii).sc3.muGap=muGap;
    threeTest(ii).sc3.uncGap=uncGap;
    threeTest(ii).sc3.H0=H0;
    threeTest(ii).sc3.xm=mean(xAll,2);
    threeTest(ii).sc3.tRun=toc;
end
scriptName=mfilename;
save(savePath,'threeTest','scriptName','scenarios')
end

function [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap,H0]=returnscenarioinfo(xAll,sigY,A,dateM,outDat)
%Return observation model values
Aout=mean(A,3);
sigYOut=mean(sigY,2);
AUnc=std(A,0,3);
sigYUnc=std(sigY,0,2);
H0=outDat.H0;

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



