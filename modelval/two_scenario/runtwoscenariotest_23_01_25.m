%Run two scenario synthetic data test in parallel computing environment
% Ted Amdur
% 2023/01/25

clearvars
parpool('local',str2num(getenv('SLURM_CPUS_PER_TASK')))
%Make sure everything is visible
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/tools')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/tools/BoE')
addpath('/net/rcstorenfs02/ifs/rc_labs/huybers_lab/tamdur/ACRIM-Gap/mat_files')
rng(1)

%Load the synthetic datasets to be examined
load 2scenario_23_01_25b.mat
tN=length(ACRIM); %Number of synthetic datasets to be inferred
parfor ii=1:tN
    tic;
    %Run BTSI on ACRIM scenario
    [xAll,sigY,~,~,~,A,~,~] = runchain_23_01_13(ACRIM(ii).valM,oM,colLabels,[],false);
    [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap]=returnscenarioinfo(xAll,sigY,A,dateM);
    twoTest(ii).ACRIM.Aout=Aout;
    twoTest(ii).ACRIM.sigYOut=sigYOut;
    twoTest(ii).ACRIM.AUnc=AUnc;
    twoTest(ii).ACRIM.sigYUnc=sigYUnc;
    twoTest(ii).ACRIM.muGap=muGap;
    twoTest(ii).ACRIM.uncGap=uncGap;
    twoTest(ii).ACRIM.tRun=toc;
    
    %Run BTSI on PMOD scenario
    [xAll,sigY,~,~,~,A,~,~] = runchain_23_01_13(PMOD(ii).valM,oM,colLabels,[],false);
    [Aout,sigYOut,AUnc,sigYUnc,muGap,uncGap]=returnscenarioinfo(xAll,sigY,A,dateM);
    twoTest(ii).PMOD.Aout=Aout;
    twoTest(ii).PMOD.sigYOut=sigYOut;
    twoTest(ii).PMOD.AUnc=AUnc;
    twoTest(ii).PMOD.sigYUnc=sigYUnc;
    twoTest(ii).PMOD.muGap=muGap;
    twoTest(ii).PMOD.uncGap=uncGap;
    twoTest(ii).PMOD.tRun=toc;
end
scriptName=mfilename;
save('twotestcluster_23_01_25b.mat','twoTest','scriptName')

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
uncGap=(quantile(BTSIAll(:,3),0.975)-quantile(BTSIAll(:,3),0.025))./1.96;
end


