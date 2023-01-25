

obsmatrix='obs_23_01_13.mat';
load(obsmatrix); %From makeobsmatrix.m
%load excludeMask22_6_29.mat %from exclude_fliers_22_04_26.m
% valM(excludeMask) = NaN;
% oM(excludeMask) = false;

inCtAll=0;
allCtAll=0;
figure
for obs=1:12
%obs=2;%observer
nSim=size(oTest(obs).x,1);
inCt=zeros(size(valM,1),1);
allCt=zeros(size(valM,1),1);
for ii=1:nSim
    iV=~isnan(oTest(obs).p(ii,:,1));
    p=squeeze(oTest(obs).p(ii,iV,:));
    iX=oTest(obs).iX(ii,iV);
    ci95=prctile(p',[2.5 97.5]);
    i95=ci95(1,:)<valM(iX,obs)' & valM(iX,obs)'<ci95(2,:);
    inCt(iX(i95))=inCt(iX(i95))+1;
    allCt(iX)=allCt(iX)+1;
end
inCtAll=inCtAll+sum(inCt);
allCtAll=allCtAll+sum(allCt);
subplot(3,4,obs)
plot(dateM,inCt./allCt,'.')
titleStr=strcat(colLabels(obs)', ': ', num2str(100.*sum(inCt)./sum(allCt)), '%');
title(titleStr)
set(gca,'FontSize',16)
end
    
    