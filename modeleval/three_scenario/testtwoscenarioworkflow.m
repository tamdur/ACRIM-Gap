%Testbed script to make sure all scripts work together for two scenario
%synthetic data test

% Ted Amdur
% 2023/01/25

clearvars

%Load the synthetic datasets to be examined
load 2scenario_23_01_25.mat

%Run the chain for a single scenario
[xAll,sigY,sigX,theta,a,A,t,outDat] = runchain_23_01_13(PMOD(5).valM,oM,colLabels,[],true);

%Get the variables of interest from this single scenario
dateS=getdates;
smoothWindow=1;
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
BTSI=[mean(xm(stBTSI));
    mean(xm(endBTSI))];
BTSI=diff(BTSI);

%Get max and min change
BTSIAll=zeros(size(tsiAll,2),3);
for ii=1:size(tsiAll,2)
    BTSIAll(ii,1:2)=[mean(tsiAll(stBTSI,ii)),mean(tsiAll(endBTSI,ii))];
    BTSIAll(ii,3)=diff(BTSIAll(ii,1:2));
end
BTSIUNC=(quantile(BTSIAll(:,3),0.975)-quantile(BTSIAll(:,3),0.025))./1.96;

    
    