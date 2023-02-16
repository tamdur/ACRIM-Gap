%Create a set of plots for evaluating BTSI long output
%
% Ted Amdur
% 11/30/22

clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each conditional, set to 1 to plot/calculate that figure/table, 0
% otherwise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
tsiComparison =0; %plot Figure 1 of manuscript
priorposterior=0; %Plot the priors and posteriors for each observer
priorposterior2=0; %Plot the priors and posteriors for each observer
obsContributions=0; %Plot the relative contribution of each observer to BTSI over time
twoScenario=0; %Plot results of synthetic data experiment for ACRIM and PMOD gaps
threeScenario=1; %twoScenario, but with ACRIM-Sat/CPMDF-Proxy scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER CALCULATIONS
gapChange=0; %Calculate change in TSI between two periods
trendUnc=0;%Calculate uncertainty in linear drift from BTSI
posteriorParams=0; %Calculate posterior parameter values and confidence interval
uncBTSI=0;%Calculate and plot the uncertainty in BTSI
table1=0; %Calculate parameter values for Table 1 of manuscript
table2=0; %Calculate values for Table 2, the posterior model parameters
tableSynthH=0; %Show observer errors used in synthetic experiment
autocorr=0; %Calculate autocorrelation of BTSI vs other TSI reconstructions
PMODCorrections=0; %Calculate and plot the corrections made by Frohlich
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fSize = 20;
load ar2_23_01_14.mat; %Select the output chain to plot/analyze
obsmatrix='obs_23_02_01.mat';

%--------------------------------------------------------------------------
% UNCOMMENT IN ORDER TO USE CARRINGTON ROTATION RATHER THAN MONTHLY SUBSETS
% load ar2carrington_23_02_02.mat; %Select the output chain to plot/analyze
% obsmatrix='obscarrington_23_02_02.mat';
%--------------------------------------------------------------------------


load(obsmatrix); %From makeobsarray.m
valAll=valM-offsets; %Remove mean offset
dateS=getdates;
dates=dateS.acrimplusfive;
%dateI=dateM>=datejd(dates(1))&dateM<datejd(dates(2)+4);
%oM=oM(dateI,:);valAll=valAll(dateI,:);dateM=dateM(dateI);

%Color scheme
c1 = [51; 136; 68; 17; 153; 221; 204; 136; 170];
c2 = [34; 204; 170; 119; 153; 204; 102; 34; 68];
c3 = [136; 238; 153; 51; 51; 119; 119; 85; 153];
c = [c1 c2 c3]; c = c./255;
c(10,:) = c(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if tsiComparison 
    %Order as proxies followed by satellites in chronological order
    lI=[6;3;5;1;4;2;7];
    A=A(lI,:,:);
    colLabels=colLabels(lI);
    offsets=offsets(lI);
    oM=oM(:,lI);
    sigY=sigY(lI,:);
    t=t(:,lI);
    valM=valM(:,lI); valMAll=valM;
    dateMAll=dateM;
    
    %Create x-axis points for cycle demarcation
    ca=datejd([dateS.acrim, fliplr(dateS.acrim)]);

    showTrend = 1;
    smoothWindow = 6; %set smoothing (months)
    
    figure2('Position',[10 10 1000 1500])
    subplot('position',[.09 .85 .85 .13]) %Plot of proxy observations
    yyaxis left
    fill(ca,[0 0 350 350],[.85 .85 .85],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    h(1)=plot(dateM,valM(:,1)+offsets(1),'.','MarkerSize',10);
    ylabel('Sunspot number')
    yyaxis right
    h(2)=plot(dateM,valM(:,2)+offsets(2),'.','MarkerSize',10);
    legend(h,'Silso sunspot number','Mg-II')
    legend boxoff
    ylabel('Mg-II index')
    xlim([datejd(dates(1)) datejd(dates(2))])
    yyaxis left
    ylim([0 300])
    set(gca,'FontSize',fSize)
    
    subplot('position',[.09 .575 .85 .25]) %Plot of satellite observations
    fill(ca,[1360 1360 1375 1375],[.85 .85 .85],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    ind = 1;
    for ii = 3:numel(lI) %Iterate over satellite observations
        if showTrend
            hold on
            [trends,offsets2]= returntrend(A,t,ii);
            tM = mean(trends,2) + mean(offsets2)+offsets(ii);
            t995 = quantile(trends,.995,2)+quantile(offsets2,.995)+offsets(ii);
            t005 = quantile(trends,.005,2)+quantile(offsets2,.005)+offsets(ii);
            x2 = [dateM(oM(:,ii))', flip(dateM(oM(:,ii)))'];
            fill(x2,[t995(oM(:,ii))',flip(t005(oM(:,ii)))'], ...
                [1 .85 .85],'FaceAlpha',0.5,'LineStyle','none');
            hold on
            plot(dateM(oM(:,ii)),tM(oM(:,ii)),'Color','r')
            hold on
        end
        plot(dateMAll,valMAll(:,ii)+offsets(ii),'o','MarkerSize',4,...
            'Color',c(ind,:));
        hold on
        hh(ind) = plot(dateM,valM(:,ii)+offsets(ii),'.','MarkerSize',10,...
            'Color',c(ind,:));
        ind = ind + 1;
    end
    xlim([datejd(dates(1)) datejd(dates(2))])
    ylim([1360 1374])
    legend(hh,colLabels(3:end))
    legend boxoff
    ylabel('TSI (W/m^{2})')
    set(gca,'FontSize',fSize)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot of reconstructions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot('position',[.09 .06 .85 .49]) 
    fill(ca,[-2 -2 2 2],[.85 .85 .85],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    cColor = get(gca,'colororder');
    %Get CI of our estimate,plot
    tsix = prctile(xAll',[.5 5 50 95 99.5])';
    for iS = 1:size(tsix,2)
        tsix(:,iS) = smoothPH(tsix(:,iS),smoothWindow);
    end
    [~,tsiO] = meaninterval(dateM,tsix(:,3),1985,1995);
    tsix = tsix-tsiO;
    x2 = [dateM', fliplr(dateM')];
    fill(x2,[tsix(:,1)',fliplr(tsix(:,end)')], [.85 .85 .85],'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[tsix(:,2)',fliplr(tsix(:,4)')], [.75 .75 .75],'FaceAlpha',0.5,...
        'LineStyle','none');
    
    %Plot other TSI reconstructions
    ind = 1;
    load oTSI_23_02_01.mat %From readothertsi.m in the code_22_06 directory
    %Plot PMOD
    hold on
    tsiA = smoothPH(oTSI(7).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(7).mthdatetime,tsiA,1985,1995);
    h(ind) = plot(oTSI(7).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(2,:));
    legendtxt(ind) = string(oTSI(7).product);
    ind = ind + 1;
    
    %Plot SOLID
    hold on
    tsiA = smoothPH(oTSI(9).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(9).mthdatetime,tsiA,1985,1995);
    h(ind) = plot(oTSI(9).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(7,:));
    legendtxt(ind) = "Comm.-Consensus Corr.";
    ind = ind + 1;
    
    %Plot SOLID Uncorrected
    hold on
    tsiA = smoothPH(oTSI(10).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(10).mthdatetime,tsiA,1985,1995);
    h(ind) = plot(oTSI(10).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(8,:));
    legendtxt(ind) = "Comm.-Consensus Uncorr.";
    ind = ind + 1;
    
    %Plot ACRIM
    hold on
    tsiA = smoothPH(oTSI(6).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(6).mthdatetime,tsiA,1985,1995);
    h(ind) = plot(oTSI(6).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(4,:));
    legendtxt(ind) = string(oTSI(6).product);
    ind=ind+1;
    

    hold on
    h(ind) = plot(dateM,tsix(:,3),'LineWidth',3,...
        'Color','k');
    legendtxt(ind) = "BTSI";

    legend(h,legendtxt)
    legend boxoff
    set(gca,'FontSize',fSize)
    xlabel('Year')
    ylabel('TSI anomaly from 1985-1995 mean (W/m^{2})')
    xlim([datejd(dates(1)) datejd(dates(2))])
    ylim([-0.9 1.25])
    saveas(gcf,'plots/tsicompare_23_02_01.png')
end
if priorposterior
    datesave='23_02_01'; %Date for figure name
    satindex=outDat.satindex;
    obsUsed=satindex;
    obsUsed(3)=true; %Turn on Bremen Mg-II
    obsUsed(6)=true; %Turn on SILSO spots
    satI=find(satindex);
    %Pull output from simulation
    reps = outDat.reps;
    sC = squeeze(A(:,2,:));
    bC = squeeze(A(:,1,:));
    mC = squeeze(A(:,3,:));
    aP = a;
    oP = squeeze(bC)';
    lP = squeeze(mC(satI,:))';
    rP = sigY;
    H0=outDat.H0;
    Hsig=outDat.Hsig;
    T0=outDat.T0;
    th0=outDat.th0;
    %Make the data that goes in the varpdf fields for each of 3 plots
    
    %First, do the offset variables
    varNames1 = ["o1";"o2";"o4";"o5"];
    th1M1 = H0(satI(1:end-1),1);
    th2M1 = sqrt(Hsig(satI(1:end-1),1));
    vals1 = oP(:,satI(1:end-1));
    %Then, do the drift variables
    varNames2 = ["l1";"l2";"l4";"l5";"l7"];
    th1M2 = H0(satI,3);
    th2M2 = sqrt(Hsig(satI,3));
    vals2 = lP;
    %Last, do the noise variables
    nN=length(colLabels(satI)); %Number of noisy observers
    varNames3=[];
    for ii = 1:nN
        varNames3 = [varNames3;strcat("\epsilon_{",colLabels(satI(ii)),"}")];
    end
    th1M3 = T0(satI);
    th2M3 = th0(satI);
    vals3 = rP(satI,:)'; 
    vals3(vals3 <= 0) = NaN;
    vals3 = vals3.^0.5;%Note 0.5 is choice to show std v var
    distType3 = repmat("invgamma",[nN 1]);

    %------------------------------------------------------------------
    % First, plot estimated offsets 
    %------------------------------------------------------------------
    figure2('Position',[10 10 1600 600])
    %subplot('position',[.09 .7 .85 .27])
    offsetsI = satI';
    numPlots = length(varNames1);
    pDim = ceil(sqrt(numPlots));
    for ii = 1:numPlots
        subplot(pDim,pDim,ii)
        hold on
        th1 = th1M1(ii) + offsets(offsetsI(ii));
        th2 = th2M1(ii);
        %Plot the prior distribution
        x = vals1;
        xL = min(x);xH = max(x);
        xPrior = randraw('norm', [th1 th2], [1e5 1] );
        xL = min([xPrior; vals1(:,ii)+ offsets(offsetsI(ii))]);
        xH = max([xPrior;vals1(:,ii)+ offsets(offsetsI(ii))]);
        pEdges = linspace(xL,xH,1000);
        [yP,xP] = histcounts(xPrior,pEdges); 
        [xP,yP,~] = histtoplot(yP,xP,50);
        %Plot the prior
        plot(xP,yP);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals1(:,ii)+ offsets(offsetsI(ii)),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,10);
        hold on
        plot(xPost,yPost)
        xlim([quantile(xPrior,0.001) quantile(xPrior,0.999)])
        ylim([0 2.*max(yP)])
        title(colLabels(offsetsI(ii)))
        xlabel("W/m^{2}")
        set(gca,'ytick',[])
        set(gca,'FontSize',fSize)
    end
    saveas(gcf,['plots/priorposterior1_' datesave '.png'])
    %------------------------------------------------------------------
    % Next, plot estimated linear drifts 
    %------------------------------------------------------------------
    figure2('Position',[10 10 1600 600])
    offsetsI = satI';
    numPlots = length(varNames2);
    for ii = 1:numPlots
        subplot(2,3,ii)
        hold on
        th1 = th1M2(ii);
        th2 = th2M2(ii);
        %Plot the prior distribution
        x = vals2;
        xL = min(x);xH = max(x);
        xPrior = randraw('norm', [th1 th2], [1e5 1] );
        xL = min([xPrior; vals2(:,ii)]);
        xH = max([xPrior; vals2(:,ii)]);
        %Plot the prior
        pEdges = linspace(xL,xH,1000);
        [yP,xP] = histcounts(xPrior,pEdges);
        [xP,yP,~] = histtoplot(yP,xP,50);
        plot(xP,yP);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals2(:,ii),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,20);
        hold on 
        plot(xPost,yPost)
        title(colLabels(offsetsI(ii)))
        xlim([quantile(vals2(:,ii),0.001) quantile(vals2(:,ii),0.999)])
        xlabel("W/m^{2}/decade")
        set(gca,'ytick',[])
        set(gca,'FontSize',fSize)
    end
    saveas(gcf,['plots/priorposterior2_' datesave '.png'])
    %------------------------------------------------------------------
    % Last, plot noise estimates
    %------------------------------------------------------------------
    figure2('Position',[10 10 1600 600])
    numPlots = length(varNames3);
    for ii = 1:numPlots
        subplot(2,3,ii)
        hold on
        th1 = th1M3(ii);
        th2 = th2M3(ii);
        %Plot the prior distribution
        x = vals3;
        xL = min(x);xH = max(x);
        %From randraw.m, form of drawing is y = randraw('gamma', [1 3 2], [1e5 1] );
        %where inputs are location, shape, and scale parameters,
        %respectively. Second input is number of draws
        if strcmp(distType3(ii),"normal")
            xPrior = randraw('norm', [th1 th2], [1e5 1] );
            xPrior(xPrior < 0) = NaN;
        else
            xPrior = drawgamma(0, th1, th2, 1e5);
        end
        xPrior = xPrior.^0.5;
        xL = max([min([xL(ii); min(xPrior)]); 1e-5]);
        xH = max([xH(ii); xPrior]);
        %Plot the prior
        pEdges = logspace(log10(xL),log10(xH),1000);
        [yP,xP] = histcounts(xPrior,pEdges);
        [xP,yP,~] = histtoplot(yP,xP,50);
        plot(xP,yP);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals3(:,ii),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,50);
        hold on
        plot(xPost,yPost)
        title(varNames3(ii))
        xmin = max([min([quantile(vals3(:,ii),0.001); quantile(xPrior,0.01)]); 1e-5]);
        xmax = max([quantile(vals3(:,ii),0.999); quantile(xPrior,0.99)]);
        xlim([xmin xmax])
        set(gca,'ytick',[])
        set(gca,'Xtick',linspace(xmin,xmax,4))
        xtickformat('%.2f')
        xlabel("W/m^{2}")
        set(gca,'FontSize',fSize)
    end
    saveas(gcf,['plots/priorposterior3_' datesave '.png'])
end
if priorposterior2
    datesave='23_02_14'; %Date for figure name
    satindex=outDat.satindex;
    obsUsed=satindex;
    obsUsed(3)=true; %Turn on Bremen Mg-II
    obsUsed(6)=true; %Turn on SILSO spots
    satI=find(satindex);
    %Pull output from simulation
    reps = outDat.reps;
    sC = squeeze(A(:,2,:));
    bC = squeeze(A(:,1,:));
    mC = squeeze(A(:,3,:));
    aP = a;
    oP = squeeze(bC)';
    lP = squeeze(mC(satI,:))';
    rP = sigY;
    H0=outDat.H0;
    Hsig=outDat.Hsig;
    T0=outDat.T0;
    th0=outDat.th0;
    %Make the data that goes in the varpdf fields for each of 3 plots
    
    %First, do the offset variables
    th1M1 = H0(satI(1:end-1),1);
    th2M1 = sqrt(Hsig(satI(1:end-1),1));
    vals1 = oP(:,satI(1:end-1));
    %Then, do the drift variables
    varNames2 = ["l1";"l2";"l4";"l5";"l7"];
    th1M2 = H0(satI,3);
    th2M2 = sqrt(Hsig(satI,3));
    vals2 = lP;
    %Last, do the noise variables
    nN=length(colLabels(satI)); %Number of noisy observers
    varNames3=[];
    for ii = 1:nN
        varNames3 = [varNames3;strcat("\epsilon_{",colLabels(satI(ii)),"}")];
    end
    th1M3 = T0(satI);
    th2M3 = th0(satI);
    vals3 = rP(satI,:)'; 
    vals3(vals3 <= 0) = NaN;
    vals3 = vals3.^0.5;%Note 0.5 is choice to show std v var
    distType3 = repmat("invgamma",[nN 1]);

    %------------------------------------------------------------------
    % First, plot estimated offsets 
    %------------------------------------------------------------------
    figure2('Position',[10 10 1000 1000])
    subplot('position',[.09 .73 .85 .27])
    offsetsI = satI';
    for ii = 1:4
        hold on
        th1 = th1M1(ii) + offsets(offsetsI(ii));
        th2 = th2M1(ii);
        %Plot the prior distribution
        x = vals1;
        xL = min(x);xH = max(x);
        xPrior = randraw('norm', [th1 th2], [1e5 1] );
        xL = min([xPrior; vals1(:,ii)+ offsets(offsetsI(ii))]);
        xH = max([xPrior;vals1(:,ii)+ offsets(offsetsI(ii))]);
        pEdges = linspace(xL,xH,1000);
        [yP,xP] = histcounts(xPrior,pEdges); 
        [xP,yP,~] = histtoplot(yP,xP,50);
        %Plot the prior
        h(ii)=plot(xP,yP,'Color',c(2*ii,:),'LineWidth',2);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals1(:,ii)+ offsets(offsetsI(ii)),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,10);
        hold on
        plot(xPost,yPost,'Color',c(2*ii,:),'LineWidth',2)
    end
    h(5)=line([offsets(7) offsets(7)],[0 6],'Color',c(10,:),'LineWidth',2);
    legend(h,colLabels(offsetsI),'Location','NorthWest')
    legend boxoff
    xlabel("W/m^{2}")
    set(gca,'ytick',[])
    set(gca,'FontSize',fSize)
    xlim([1357 1372])
    %------------------------------------------------------------------
    % Next, plot estimated linear drifts 
    %------------------------------------------------------------------
    subplot('position',[.09 .39 .85 .27])
    offsetsI = satI';
    numPlots = length(varNames2);
    for ii = 1:numPlots
        hold on
        th1 = th1M2(ii);
        th2 = th2M2(ii);
        %Plot the prior distribution
        x = vals2;
        xL = min(x);xH = max(x);
        xPrior = randraw('norm', [th1 th2], [1e5 1] );
        xL = min([xPrior; vals2(:,ii)]);
        xH = max([xPrior; vals2(:,ii)]);
        %Plot the prior
        pEdges = linspace(xL,xH,200);
        [yP,xP] = histcounts(xPrior,pEdges);
        [xP,yP,~] = histtoplot(yP,xP,50);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals2(:,ii),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,20);
        hold on 
        h(ii)=plot(xPost,yPost,'Color',c(2*ii,:),'LineWidth',3);
    end
    h(ii+1)=plot(xP,yP,'Color','k','LineWidth',5);
    legendtxt=[colLabels(offsetsI);"Prior Distribution"];
    legend(h,legendtxt,'Location','NorthWest')
    legend boxoff
    %xlim([quantile(vals2(:,ii),0.001) quantile(vals2(:,ii),0.999)])
    xlabel("W/m^{2}/decade")
    set(gca,'ytick',[])
    set(gca,'FontSize',fSize)
    xlim([-1 1])
    %------------------------------------------------------------------
    % Last, plot noise estimates
    %------------------------------------------------------------------
    subplot('position',[.09 .08 .85 .24])
    numPlots = length(varNames3);
    for ii = 1:numPlots
        hold on
        th1 = th1M3(ii);
        th2 = th2M3(ii);
        %Plot the prior distribution
        x = vals3;
        xL = min(x);xH = max(x);
        %From randraw.m, form of drawing is y = randraw('gamma', [1 3 2], [1e5 1] );
        %where inputs are location, shape, and scale parameters,
        %respectively. Second input is number of draws
        if strcmp(distType3(ii),"normal")
            xPrior = randraw('norm', [th1 th2], [1e5 1] );
            xPrior(xPrior < 0) = NaN;
        else
            xPrior = drawgamma(0, th1, th2, 1e5);
        end
        xPrior = xPrior.^0.5;
        xL = max([min([xL(ii); min(xPrior)]); 1e-5]);
        xH = max([xH(ii); xPrior]);
        pEdges = logspace(log10(xL),log10(xH),200);
        %Plot the posterior distribution
        [yPost,xPost] = histcounts(vals3(:,ii),pEdges); 
        [xPost,yPost,~] = histtoplot(yPost,xPost,50);
        hold on
        plot(xPost,yPost,'Color',c(2*ii,:),'LineWidth',3);
    end
    legend(varNames3,'Location','NorthWest')
    legend boxoff
    xlim([0 0.4])
    set(gca,'ytick',[])
    set(gca,'Xtick',linspace(0,0.4,5))
    xtickformat('%.2f')
    xlabel("W/m^{2}")
    set(gca,'FontSize',fSize)
    saveas(gcf,['plots/priorposterior_' datesave '.png'])
end
if obsContributions
     %Plot the relative contribution of each observer to the estimate of TSI
    smoothWindow = 36;
    conChain=outDat.contributionChain;
    %Reorient
    lI=[6;3;5;1;4;2;7];
    A=A(lI,:,:);
    colLabels=colLabels(lI);
    offsets=offsets(lI);
    oM=oM(:,lI);
    sigY=sigY(lI,:);
    t=t(:,lI);
    valM=valM(:,lI);
    conChain=conChain(:,:,lI);%Reorient to agree with formatting of this plotting function
    cn=squeeze(mean(abs(conChain),1));
    for ii=1:size(valM,2)
        cn(:,ii)=movmean(cn(:,ii),smoothWindow,'omitnan');
        cn(~oM(:,ii),ii)=NaN;
    end
    cn=cn./nansum(cn,2);
    
    %Plot
    figure2('Position',[10 10 1600 700])
    %Create x-axis points for cycle demarcation
    ca=datejd([dateS.acrim, fliplr(dateS.acrim)]);
    fill(ca,[0 0 0.5 0.5],[.85 .85 .85],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    ind=1;
    for ii=1:2
        h(ind)=plot(dateM,cn(:,ind),'--','LineWidth',4);
        hold on
        ind=ind+1;
    end
    for ii=1:5
        h(ind)=plot(dateM,cn(:,ind),'LineWidth',3);
        hold on
        ind=ind+1;
    end
    legend(h,colLabels,'Location','NorthWest','NumColumns',2)
    legend boxoff
    xlabel('Year')
    ylabel('Fractional Contribution')
    set(gca,'FontSize',fSize)
    saveas(gcf,'plots/obscontribution_23_01_18.png')
    
end
if twoScenario
    fSize=16;
    %Plot a panel with the error structure of the satellites
    %load 2scenario_23_01_25b.mat %Get hyperparameters, main structure
    %load 2scenario_23_01_31_PMODproxy.mat
    load 2scenario_23_01_31_PMODProxyb.mat
    tsi = twotsiseries;
    Ainit=setInfo.Ainit;
    figure2('Position',[10 10 1150 1000])
    subplot('position',[.09 .635 .85 .35]) %Plot of proxy observations
    %First, plot the error structure sans AR(1) 
    satIndex=[1;2;4;5;7];
    ind=1;
    for ii=1:length(satIndex)
        pred=Ainit(satIndex(ii),1)+Ainit(satIndex(ii),3).*t(:,satIndex(ii));
        pred=pred(oM(:,satIndex(ii)));
        h(ind)=plot(dateM(oM(:,satIndex(ii))),pred,'Color',c(2*ii,:),'LineWidth',2);
        hold on
        plot(dateM,ACRIM(1).valM(:,satIndex(ii))-(tsi.ACRIM-nanmean(tsi.ACRIM)),'--','Color',c(2*ii,:))
        hold on
        ind=ind+1;
    end
    legend(h,colLabels(satIndex),'Location','NorthWest')
    legend boxoff
    xlabel('Year')
    ylabel('Satellite observer error (W/m^{2})')
    set(gca,'FontSize',fSize)
    
    %Plot a panel with the ACRIM Gap outputs
    clear h
    
    %load twotestcluster_23_01_25b.mat
    %load twotestcluster_ACRIMsatPMODprox_23_01_31.mat
    load twotestcluster_23_01_31.mat
    PMODGAP=0.0159; %FROM THE gapChange calculation
    ACRIMGAP=0.7057; %From the gapChange calculation
    %Plot a panel with the correct ACRIM-Gap for ACRIM, what the model finds
    ACRIM=struct;
    PMOD=struct;
    for ii=1:size(twoTest,2)
        ACRIM.gap(ii)=twoTest(ii).ACRIM.muGap;
        ACRIM.gapUnc(ii)=twoTest(ii).ACRIM.uncGap;
        ACRIM.sigY1(ii,:)=twoTest(ii).ACRIM.sigYOut;
        PMOD.gap(ii)=twoTest(ii).PMOD.muGap;
        PMOD.gapUnc(ii)=twoTest(ii).PMOD.uncGap;
        PMOD.sigY1(ii,:)=twoTest(ii).PMOD.sigYOut;
        %See if gap is within expected uncertainty
        ACRIM.within(ii)=abs(ACRIM.gap(ii)-ACRIMGAP)<(2.*ACRIM.gapUnc(ii));
        PMOD.within(ii)=abs(PMOD.gap(ii)-PMODGAP)<(2.*PMOD.gapUnc(ii));
    end
    %Run again for ACRIM satellite/PMOD composite
    load twotestcluster_ACRIMsatPMODprox_23_01_31.mat
    for ii=1:size(twoTest,2)
        AP.gap(ii)=twoTest(ii).ACRIM.muGap;
        AP.gapUnc(ii)=twoTest(ii).ACRIM.uncGap;
        AP.sigY1(ii,:)=twoTest(ii).ACRIM.sigYOut;
        %See if gap is within expected uncertainty
        AP.within(ii)=abs(ACRIM.gap(ii)-ACRIMGAP)<(2.*ACRIM.gapUnc(ii));
    end
    subplot('position',[.09 .07 .85 .5]) %Plot of proxy observations
    h(1)=histogram(ACRIM.gap,'BinWidth',0.025);
    hold on
    line([ACRIMGAP ACRIMGAP],[0 180],'LineWidth',2,'Color',[ 0    0.4470    0.7410]);
    disp(['ACRIM 95\% CI:' num2str(sum(ACRIM.within)./size(twoTest,2),'%.2f')])
    hold on
    h(2)=histogram(PMOD.gap,'BinWidth',0.025);
    hold on
    line([PMODGAP PMODGAP],[0 180],'LineWidth',2,'Color',[0.8500    0.3250    0.0980]);
    hold on
    h(3)=histogram(AP.gap,'BinWidth',0.025,'FaceColor',[1 0.01 0.01],...
        'FaceAlpha',0.5,'LineStyle','--','LineWidth',1);
    
    legend(h,'ACRIM-All','CPMDF-All','ACRIM-Satellite/CPMDF-Proxy','Location','NorthWest')
    legend boxoff
    disp(['PMOD 95\% CI:' num2str(sum(PMOD.within)./size(twoTest,2),'%.2f')])
    xlabel('ACRIM-Gap magnitude (W/m^{2})')
    ylabel('Number of simulations')
    set(gca,'FontSize',fSize)
    xlim([-0.36 1])
    
    %Calculate the statistical power of BTSI in inferring ARIM Gap
    PMOD_sig=quantile(PMOD.gap,.975);
    sig_ct=sum(AP.gap>PMOD_sig);
    stat_power=sig_ct./length(ACRIM.gap);
    saveas(gcf,'plots/twoscenario_23_02_11.png')
    
end
if threeScenario
    fSize=16;
    %Plot a panel with the error structure of the satellites
    %Load synthetic datasets and results
    load 3scenario_23_02_15.mat
    tsi = twotsiseries;
    Ainit=setInfo.Ainit;
    figure2('Position',[10 10 1150 1000])
    subplot('position',[.09 .635 .85 .35]) %Plot of proxy observations
    %First, plot the error structure sans AR(1) 
    satIndex=[1;2;4;5;7];
    ind=1;
    for ii=1:length(satIndex)
        pred=Ainit(satIndex(ii),1)+Ainit(satIndex(ii),3).*t(:,satIndex(ii));
        pred=pred(oM(:,satIndex(ii)));
        h(ind)=plot(dateM(oM(:,satIndex(ii))),pred,'Color',c(2*ii,:),'LineWidth',2);
        hold on
        plot(dateM,ACRIM(1).valM(:,satIndex(ii))-(tsi.ACRIM-nanmean(tsi.ACRIM)),'--','Color',c(2*ii,:))
        hold on
        ind=ind+1;
    end
    legend(h,colLabels(satIndex),'Location','NorthWest')
    legend boxoff
    xlabel('Year')
    ylabel('Satellite observer error (W/m^{2})')
    set(gca,'FontSize',fSize)
    
    %Plot a panel with the ACRIM Gap outputs
    clear h
    load threetestcluster_generic_23_02_15b.mat
    PMODGAP=0.0159; %FROM THE gapChange calculation
    ACRIMGAP=0.7057; %From the gapChange calculation
    %Plot a panel with the correct ACRIM-Gap for ACRIM, what the model finds
    ACRIM=struct;
    PMOD=struct;
    AP=struct;
    for ii=1:size(threeTest,2)
        ACRIM.gap(ii)=threeTest(ii).ACRIM.muGap;
        ACRIM.gapUnc(ii)=threeTest(ii).ACRIM.uncGap;
        ACRIM.sigY1(ii,:)=threeTest(ii).ACRIM.sigYOut;
        PMOD.gap(ii)=threeTest(ii).PMOD.muGap;
        PMOD.gapUnc(ii)=threeTest(ii).PMOD.uncGap;
        PMOD.sigY1(ii,:)=threeTest(ii).PMOD.sigYOut;
        %See if gap is within expected uncertainty
        ACRIM.within(ii)=abs(ACRIM.gap(ii)-ACRIMGAP)<(2.*ACRIM.gapUnc(ii));
        PMOD.within(ii)=abs(PMOD.gap(ii)-PMODGAP)<(2.*PMOD.gapUnc(ii));
        AP.gap(ii)=threeTest(ii).AP.muGap;
        AP.gapUnc(ii)=threeTest(ii).AP.uncGap;
        AP.sigY1(ii,:)=threeTest(ii).AP.sigYOut;
        %See if gap is within expected uncertainty
        AP.within(ii)=abs(AP.gap(ii)-ACRIMGAP)<(2.*AP.gapUnc(ii));
    end
    subplot('position',[.09 .07 .85 .5]) %Plot of proxy observations
    h(1)=histogram(ACRIM.gap,'BinWidth',0.025);
    hold on
    line([ACRIMGAP ACRIMGAP],[0 180],'LineWidth',2,'Color',[ 0    0.4470    0.7410]);
    disp(['ACRIM 95\% CI:' num2str(sum(ACRIM.within)./size(threeTest,2),'%.2f')])
    hold on
    h(2)=histogram(PMOD.gap,'BinWidth',0.025);
    hold on
    line([PMODGAP PMODGAP],[0 180],'LineWidth',2,'Color',[0.8500    0.3250    0.0980]);
    hold on
    h(3)=histogram(AP.gap,'BinWidth',0.025,'FaceColor',[1 0.01 0.01],...
        'FaceAlpha',0.5,'LineStyle','--','LineWidth',1);
    
    legend(h,'ACRIM-All','CPMDF-All','ACRIM-Satellite/CPMDF-Proxy','Location','NorthWest')
    legend boxoff
    disp(['PMOD 95\% CI:' num2str(sum(PMOD.within)./size(threeTest,2),'%.2f')])
    xlabel('ACRIM-Gap magnitude (W/m^{2})')
    ylabel('Number of simulations')
    set(gca,'FontSize',fSize)
    xlim([-0.36 1])
    
    %Calculate the statistical power of BTSI in inferring ARIM Gap
    PMOD_sig=quantile(PMOD.gap,.975);
    sig_ct=sum(AP.gap>PMOD_sig);
    stat_power=sig_ct./length(ACRIM.gap);
    saveas(gcf,'plots/twoscenario_23_02_11.png')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if gapChange
    smoothWindow = 1;
    %Load BTSI reconstruction
    tsix = prctile((xAll+offsets(7))',[.5 5 50 95 99.5])';
    tsiAll=xAll+offsets(7);
    for iS = 1:size(tsix,2)
        tsix(:,iS) = smoothPH(tsix(:,iS),smoothWindow);
    end
    xm=tsix(:,3);
    %Get other TSI reconstructions
    load oTSI_23_02_01.mat
    
    %Create time windows for year before and year after acrim gap
    stInt=datejd([dateS.acrim(1)-365 dateS.acrim(1)]);
    endInt=datejd([dateS.acrim(end) dateS.acrim(end)+365]);
    
    %Get TSI at the beginning and end interval for each
    stBTSI=dateM >= stInt(1)& dateM <=stInt(end);
    endBTSI=dateM >=  endInt(1) & dateM<=endInt(end);
    BTSI=[mean(xm(stBTSI));
        mean(xm(endBTSI))];
    BTSI=[BTSI;diff(BTSI)];
    
    %Get max and min change
    BTSIAll=zeros(size(tsiAll,2),3);
    for ii=1:size(tsiAll,2)
        BTSIAll(ii,1:2)=[mean(tsiAll(stBTSI,ii)),mean(tsiAll(endBTSI,ii))];
        BTSIAll(ii,3)=diff(BTSIAll(ii,1:2));
    end
    BTSILB=quantile(BTSIAll,0.025,1);
    BTSIUB=quantile(BTSIAll,0.975,1);
        
    
    
    %Get PMOD (7)
    xPMOD=oTSI(7).tsi;datePMOD=oTSI(7).datetime;
    PMOD=[mean(xPMOD(datePMOD >= stInt(1)& datePMOD <=stInt(end)));
        mean(xPMOD(datePMOD >=  endInt(1) & datePMOD<=endInt(end)))];
    PMOD=[PMOD;diff(PMOD)];
    
    %Get SOLID (9)
    xSOLID=oTSI(9).tsi;dateSOLID=oTSI(9).datetime;
    SOLID=[mean(xSOLID(dateSOLID >= stInt(1)& dateSOLID <=stInt(end)));
        mean(xSOLID(dateSOLID >=  endInt(1) & dateSOLID<=endInt(end)))];
    SOLID=[SOLID;diff(SOLID)];
    %Get ACRIM (6)
    xACRIM=oTSI(6).tsi;dateACRIM=oTSI(6).datetime;
    ACRIM=[mean(xACRIM(dateACRIM >= stInt(1)& dateACRIM <=stInt(end)));
        mean(xACRIM(dateACRIM >=  endInt(1) & dateACRIM<=endInt(end)))];
    ACRIM=[ACRIM;diff(ACRIM)];
end
if trendUnc
    smoothWindow = 1;
    %Load BTSI reconstruction
    tsix = prctile((xAll+offsets(7))',[.5 5 50 95 99.5])';
    tsiAll=xAll+offsets(7);
    for iS = 1:size(tsix,2)
        tsix(:,iS) = smoothPH(tsix(:,iS),smoothWindow);
    end
    xm=tsix(:,3);
    %Get other TSI reconstructions
    load oTSI_23_02_01.mat
    
    %Create time windows first and last year of analysis
    stInt=[datetime(1986,1,1) datetime(1986,12,31)];
    endInt=[datetime(1996,1,1) datetime(1996,12,31)];
    tElapsed=years(mean(endInt)-mean(stInt))./10;
    
    %Get TSI at the beginning and end interval for each
    stBTSI=dateM >= stInt(1)& dateM <=stInt(end);
    endBTSI=dateM >=  endInt(1) & dateM<=endInt(end);
    BTSI=[mean(xm(stBTSI));
        mean(xm(endBTSI))];
    BTSI=[BTSI;diff(BTSI)];
    BTSI(3)=BTSI(3)./tElapsed;
    
    %Get max and min change
    BTSIAll=zeros(size(tsiAll,2),3);
    for ii=1:size(tsiAll,2)
        BTSIAll(ii,1:2)=[mean(tsiAll(stBTSI,ii)),mean(tsiAll(endBTSI,ii))];
        BTSIAll(ii,3)=diff(BTSIAll(ii,1:2));
    end
    BTSIAll(:,3)=BTSIAll(:,3)./tElapsed;
    BTSILB=quantile(BTSIAll,0.025,1);
    BTSIUB=quantile(BTSIAll,0.975,1);
        
    
    
    %Get PMOD (7)
    xPMOD=oTSI(7).tsi;datePMOD=oTSI(7).datetime;
    PMOD=[mean(xPMOD(datePMOD >= stInt(1)& datePMOD <=stInt(end)));
        mean(xPMOD(datePMOD >=  endInt(1) & datePMOD<=endInt(end)))];
    PMOD=[PMOD;diff(PMOD)];
    
    %Get SOLID (9)
    xSOLID=oTSI(9).tsi;dateSOLID=oTSI(9).datetime;
    SOLID=[mean(xSOLID(dateSOLID >= stInt(1)& dateSOLID <=stInt(end)));
        mean(xSOLID(dateSOLID >=  endInt(1) & dateSOLID<=endInt(end)))];
    SOLID=[SOLID;diff(SOLID)];
    %Get ACRIM (6)
    xACRIM=oTSI(6).tsi;dateACRIM=oTSI(6).datetime;
    ACRIM=[mean(xACRIM(dateACRIM >= stInt(1)& dateACRIM <=stInt(end)));
        mean(xACRIM(dateACRIM >=  endInt(1) & dateACRIM<=endInt(end)))];
    ACRIM=[ACRIM;diff(ACRIM)];
end
if posteriorParams
    alpha=0.95; %Specify width of confidence interval
    expA=squeeze(mean(A,3));
    intA=quantile(A,1-(1-alpha)/2,3)-quantile(A,(1-alpha)/2,3);
end
if uncBTSI
    alpha=0.95; %CI width for analysis
    tsix = quantile(xAll,[(1-alpha)/2,1-(1-alpha)/2],2);
    figure
    plot(dateM,tsix(:,2)-tsix(:,1))
end
if table1
    %Order as proxies followed by satellites in chronological order
    lI=[6;3;5;1;4;2;7];
    A=A(lI,:,:);
    colLabels=colLabels(lI);
    offsets=offsets(lI);
    oM=oM(:,lI);
    sigY=sigY(lI,:);
    theta=theta(lI,:);
    t=t(:,lI);
    valM=valM(:,lI); valMAll=valM;
    dateMAll=dateM;
    %Pull output from simulation
    reps = outDat.reps;
    sC = squeeze(A(:,2,:));
    bC = squeeze(A(:,1,:));
    mC = squeeze(A(:,3,:));
    aP = a;
    oP = squeeze(bC)';
    rP = sigY;
    H0=outDat.H0(lI);
    Hsig=outDat.Hsig(lI);
    T0=outDat.T0(lI);
    th0=outDat.th0(lI)';
    
    obs=sum(oM,1);
    T=T0+obs;
    theta=th0+mean(theta,2);
    
    %To calculate prior expectation, use mean relationship for
    %inverse-gamma: mu=theta/(T-1)
    e0=th0./(T0'-1);
    epsilon=theta./(T'-1);
end
if table2
    clear tabout
    offsets([5,1,4,2,7])=offsets([5,1,4,2,7])-offsets(7);
    table2order=[5;1;4;2;7;6;3];
    for ii=1:length(colLabels)
        ind=table2order(ii);
        Asub=squeeze(A(ind,:,:));
        tabout(:,ii)=[prctile(Asub(1,:),[2.5 50 97.5])'+offsets(ind);prctile(Asub(3,:),[2.5 50 97.5])';...
            prctile(Asub(2,:),[2.5 50 97.5])';...
            prctile(sqrt(sigY(ind,:))./Asub(2,:),[2.5 50 97.5])'];
    end
    rows=["mean 2.5";"mean 50";"mean 97.5";"drift 2.5";"drift 50";"drift 97.5";"scaling 2.5";"scaling 50";"scaling 97.5";"error 2.5";"error 50";"error 97.5"];
    displayTable=array2table(tabout,'VariableNames',colLabels(table2order));
end
if tableSynthH
    load 3scenario_23_02_15.mat
    table2order=[5;1;4;2;7;6;3];
    offsets([5,1,4,2,7])=offsets([5,1,4,2,7])-offsets(7);
    Ainit=setInfo.Ainit;
    Ainit([5,1,4,2,7],1)=Ainit([5,1,4,2,7],1)+offsets([5,1,4,2,7])';
    Ainit=Ainit(table2order,:);
end
if autocorr
    smoothWindow=1;
    xm=mean(xAll,2);
    xm=smoothPH(xm-nanmean(xm),smoothWindow);
    rhoBTSI=rhoAR1(xm(2:end),xm(1:end-1));
    %Assess other reconstructions
    %indices: SOLID 9 PMOD CMDF7 ACRIM 6
    load oTSI_23_02_01.mat %From readothertsi.m in the code_22_06 directory
    stDate=dateM(1);endDate=dateM(end);
    SOLIDDate=find(oTSI(9).mthdatetime>=stDate & oTSI(9).mthdatetime<=endDate);
    PMODDate=find(oTSI(7).mthdatetime>=stDate & oTSI(7).mthdatetime<=endDate);
    ACRIMDate=find(oTSI(6).mthdatetime>=stDate & oTSI(6).mthdatetime<=endDate);
    SOLID=oTSI(9).mthtsi(SOLIDDate);SOLID=smoothPH(SOLID-nanmean(SOLID),smoothWindow);
    PMOD=oTSI(7).mthtsi(PMODDate);PMOD=smoothPH(PMOD-nanmean(PMOD),smoothWindow);
    ACRIM=oTSI(6).mthtsi(ACRIMDate);ACRIM=smoothPH(ACRIM-nanmean(ACRIM),smoothWindow);
    rhoSOLID=rhoAR1(SOLID(2:end),SOLID(1:end-1));
    rhoPMOD=rhoAR1(PMOD(2:end),PMOD(1:end-1));
    rhoACRIM=rhoAR1(ACRIM(2:end),ACRIM(1:end-1));
    rSOLID=xcorr(SOLID,1,'coeff');rSOLID=rSOLID(3);
    rPMOD=xcorr(PMOD,1,'coeff');rPMOD=rPMOD(3);
    rACRIM=xcorr(ACRIM,1,'coeff');rACRIM=rACRIM(3);
    rBTSI=xcorr(xm,1,'coeff');rBTSI=rBTSI(3);
end
if PMODCorrections
    %Outline:
    %-Load obsmatrix for PMOD corrected observations
    load obsPMOD_23_02_01.mat
    valP=valM;
    %-Load obsmatrix for observations
    load obs_23_02_01.mat
    valO=valM;
    %-Plot a difference plot of the two observational matrices
    figure2
    plot(dateM,valO-valP);
    legend(colLabels)
    %Calculate the drift by doing linear fit
    valD=valO-valP;
    valD(1:20,1)=NaN;oM(1:20,1)=false; %Turn off first 20 months of ACRIM1
    for ii=1:length(colLabels)
        pred=[ones(sum(oM(:,ii)),1) t(oM(:,ii),ii)];
        b=regress(valD(oM(:,ii),ii),pred);
        linD(ii)=b(2);
    end
        
    
    
end
    
    


function xout=sorceanomaly(x,sorceI)
for ii=1:size(x,2)
    xout(:,ii)=x(:,ii)-mean(x(sorceI,ii));
    xout(xout(:,ii)<-1.*mean(x(sorceI,ii)),ii)=-1.*mean(x(sorceI,ii));
end
end
function rho=rhoAR1(xt,xtminus1)
    rho=(xtminus1'*xt)/(xtminus1'*xtminus1);
end
% function rho=a
    
    