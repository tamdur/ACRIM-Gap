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
tsiComparison = 0; %plot Figure 1 of manuscript
priorposterior=0; %Plot the priors and posteriors for each observer
obsContributions=0; %Plot the relative contribution of each observer to BTSI over time
twoScenario=1; %Plot results of synthetic data experiment for ACRIM and PMOD gaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER CALCULATIONS
gapChange=0; %Calculate change in TSI between two periods
posteriorParams=0; %Calculate posterior parameter values and confidence interval
uncBTSI=0;%Calculate and plot the uncertainty in BTSI
table1=0; %Calculate parameter values for Table 1 of manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nonlinearrelation=0; %Plot the nonlinear relationship for each TSI contributor
pltChain=0; %Plot the chain of a chosen variable
pltTimeseries=0; %Plot the time series for TSI, other time series
obsContTSI=0; %Plot the contribution from each observer to TSI
obsContFac=0; %Plot the contribution from each observer to Facular brightening
obsContSpt=0; %Plot the contribution from each observer to Sunspots
obsContOSF=0; %Plot the contribution from each observer to OSF
residualFac=0; %Plot the residuals of facular indices vs the hidden state process in time
residualSpt=0; %Plot the residuals of sunspot indices vs the hidden state process in time
residualOsf=0; %Plot the residuals of osf indices vs the hidden state process in time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS FOR AGU 2022
modernTSI=0; %Plot of modern TSI challenges
pltAllTimeseries=0; %Plot other time series and BTSI on same plot
pltOtherTimeseries=0; %Plot other time series of solar forcing
pltOSFCompare=0; %Show overlap of TSI with uncertain OSF record
squareSatObs=0; %Plot just the Satellites of modern TSI plot
squareOsf=0; %Plot OSF with inferred process curve
squareNonlinear=0; %Plot of just one of the nonlinearrelation components
squareTSITimeSeries=0; %Plot just time series of BTSI with comp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fSize = 20;
load ar2_23_01_14.mat; %Select the output chain to plot/analyze
obsmatrix='obs_23_01_13.mat';
load(obsmatrix); %From makeobsarray.m
valAll=valM-offsets; %Remove mean offset
dateS=getdates;
dates=dateS.acrimplusfive;
dateI=dateM>=datejd(dates(1))&dateM<datejd(dates(2));
oM=oM(dateI,:);valAll=valAll(dateI,:);dateM=dateM(dateI);

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
    load oTSI_22_10_27.mat %From readothertsi.m in the code_22_06 directory
    %Plot PMOD
    hold on
    tsiA = smoothPH(oTSI(7).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(7).mthdatetime,tsiA,1985,1995);
    h(ind) = plot(oTSI(7).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(7,:));
    legendtxt(ind) = string(oTSI(7).product);
    ind = ind + 1;
    
    %Plot SOLID
    hold on
    tsiA = smoothPH(oTSI(9).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(9).mthdatetime,tsiA,1985,1995);
    h(ind) = plot(oTSI(9).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(8,:));
    legendtxt(ind) = "Comm.-Consensus";
    ind = ind + 1;
    
    %Plot ACRIM
    hold on
    tsiA = smoothPH(oTSI(6).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(6).mthdatetime,tsiA,1985,1995);
    h(ind) = plot(oTSI(6).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(5,:));
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
    saveas(gcf,'plots/tsicompare_23_01_14.png')
end
if priorposterior
    datesave='23_01_14'; %Date for figure name
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
    %Plot a panel with the error structure of the satellites
    load 2scenario_23_01_25.mat %Get hyperparameters, main structure
    figure2
    %First, plot the error structure sans AR(1) 
    satIndex=[1;2;4;5;7];
    for ii=1:length(satIndex)
        pred=Ainit(satIndex(ii),1)+Ainit(satIndex(ii),3).*t(:,satIndex(ii));
        pred=pred(oM(:,satIndex(ii)));
        plot(dateM(oM(:,satIndex(ii))),pred,'Color',c(2*ii,:),'LineWidth',2)
        hold on
        plot(dateM,ACRIM(1).valM(:,satIndex(ii))-nanmean(ACRIM(1).valM(:,satIndex(ii))),'--','Color',c(2*ii,:))
        hold on
    end
    
    load twotest_23_01_25.mat
    PMODGAP=0.0159; %FROM THE gapChange calculation
    ACRIMGAP=0.7057; %From the gapChange calculation
    %Plot a panel with the correct ACRIM-Gap for ACRIM, what the model finds
    ACRIM=struct;
    PMOD=struct;
    for ii=1:size(twoTest,2)
        ACRIM.gap(ii)=twoTest(ii).ACRIM.muGap;
        ACRIM.gapUnc(ii)=twoTest(ii).ACRIM.uncGap;
        PMOD.gap(ii)=twoTest(ii).PMOD.muGap;
        PMOD.gapUnc(ii)=twoTest(ii).PMOD.uncGap;
    end
    figure2
    histogram(ACRIM.gap)
    hold on
    line([ACRIMGAP ACRIMGAP],[0 size(twoTest,2)./10],'LineWidth',2)
    figure2
    histogram(PMOD.gap)
    hold on
    line([PMODGAP PMODGAP],[0 size(twoTest,2)./10],'LineWidth',2)
    
    
    
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
    load oTSI_22_10_27.mat
    
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other Plots 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nonlinearrelation
    c95=[.65 .65 .65];
    c80=[.4 .4 .4];
    xbar=mean(outDat.xbar);
    YEO21=1358.7; %Solar constant from Yeo et al. 2021
    sc=YEO21;
    %Get mean tsi output
    x=mean(xAll,2);
    cmpStYr=1978; %cmpStYr=outDat.cmpStYr
    xPlot=x(dateM.Year>=cmpStYr)+sc;
    regX=(sc:0.01:sc+5)'; %input for linear regression of indices
    %Get mean output for each process after the cutoff date
    %Get facular brightening contributions
    infI=dateM(outDat.dateFac).Year>=cmpStYr;
    yFac=outDat.facAll(infI,:);
    yFacMean=mean(yFac,2);
    facMu=offsets(4); %Bremen offset
    facStd=nanstd(valAll(:,4));
    ffFac=squeeze(outDat.ff(:,4,:));
    %Get sunspot darkening contributions
    infI=dateM(outDat.dateSpt).Year>=cmpStYr;
    ySpt=outDat.sptAll(infI,:);
    ySptMean=mean(ySpt,2);
    sptMu=offsets(17); %SILSO offset
    sptStd=nanstd(valAll(:,17));
    ffSpt=squeeze(outDat.ff(:,9,:));
    %Get OSF contributions
    infI=dateM(outDat.dateOsf).Year>=cmpStYr;
    osfMu=offsets(8); %Krivova OMF offset
    osfStd=nanstd(valAll(:,8));
    yOsf=outDat.osfAll(infI,:);
    yOsfMean=mean(yOsf,2);
    ffOsf=squeeze(outDat.ff(:,6,:));
    
    
    figure2('Position',[10 10 1000 1000])
    %Plot facular contribution
    subplot(1,3,1)
    plot(xPlot,yFacMean.*facStd+facMu,'.')
    hold on
    xx=xbar-5:0.01:xbar+3;
    glm= @(b,x) b(4)+(b(1)-b(4))./(1+exp(b(2)*(x-(xbar+b(3)))));
    for ii=1:size(ffFac,1)
        aFac(:,ii)=glm(ffFac(ii,:),xx);
    end
    aFacPrctile=prctile(aFac',[2.5 10 50 90 97.5])';
    facMu=offsets(4); %Bremen offset
    facStd=nanstd(valAll(:,4));
    aFacPrctile=aFacPrctile*facStd + facMu;
    x2 = [xx+sc, fliplr(xx)+sc];
    fill(x2,[aFacPrctile(:,1)',fliplr(aFacPrctile(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[aFacPrctile(:,2)',fliplr(aFacPrctile(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    %Plot linear regression
    pred=[ones(length(xPlot),1) xPlot];
    b=regress(yFacMean,pred);
    plotPred=[ones(length(regX),1) regX];
    plot(regX,(plotPred*b).*facStd+facMu,'--','Color','k','LineWidth',2)
    hold on
    plot(xx+sc,aFacPrctile(:,3),'Color','r','LineWidth',2)
    xlabel('TSI (W/m^{2})')
    ylabel('Facular Index')
    set(gca,'FontSize',fSize)
    xlim([min(xx+sc)+3 max(xx+sc)])
    xticks([1359:2:1363])
    ylim([-0.1 0.8])
    
    %Plot sunspot contribution
    subplot(1,3,2)
    plot(xPlot,ySptMean*sptStd+sptMu,'.')
    hold on
    glm= @(b,x) b(4)+(b(1)-b(4))./(1+exp(b(2)*(x-(xbar+b(3)))));
    for ii=1:size(ffSpt,1)
        aSpt(:,ii)=glm(ffSpt(ii,:),xx);
    end
    aSptPrctile=(prctile(aSpt',[2.5 10 50 90 97.5])')*sptStd+sptMu;
    x2 = [xx, fliplr(xx)]+sc;
    fill(x2,[aSptPrctile(:,1)',fliplr(aSptPrctile(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[aSptPrctile(:,2)',fliplr(aSptPrctile(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
     %Plot linear regression
    pred=[ones(length(xPlot),1) xPlot];
    b=regress(ySptMean,pred);
    plotPred=[ones(length(regX),1) regX];
    plot(regX,(plotPred*b).*sptStd+sptMu,'--','Color','k','LineWidth',2)
    hold on
    plot(xx+sc,aSptPrctile(:,3),'Color','r','LineWidth',2)
    xlabel('TSI (W/m^{2})')
    ylabel('International Sunspot Number')
    xlim([min(xx+sc)+3 max(xx+sc)])
    ylim([-50 300])
    xticks([1359:2:1363])
    set(gca,'FontSize',fSize)
    
    %Plot osf contribution
    subplot(1,3,3)
    plot(xPlot,yOsfMean*osfStd+osfMu,'.')
    hold on
    glm= @(b,x) b(4)+(b(1)-b(4))./(1+exp(b(2)*(x-(xbar+b(3)))));
    for ii=1:size(ffSpt,1)
        aOsf(:,ii)=glm(ffOsf(ii,:),xx);
    end
    aOsfPrctile=(prctile(aOsf',[2.5 10 50 90 97.5])')*osfStd+osfMu;
    x2 = [xx, fliplr(xx)]+sc;
    fill(x2,[aOsfPrctile(:,1)',fliplr(aOsfPrctile(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[aOsfPrctile(:,2)',fliplr(aOsfPrctile(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    %Plot linear regression
    pred=[ones(length(xPlot),1) xPlot];
    b=regress(yOsfMean,pred);
    plotPred=[ones(length(regX),1) regX];
    plot(regX,(plotPred*b).*osfStd+osfMu,'--','Color','k','LineWidth',2)
    hold on
    plot(xx+sc,aOsfPrctile(:,3),'Color','r','LineWidth',2)
    xlabel('TSI (W/m^{2})')
    ylabel('Open Solar Flux (10^{14} Wb)')
    xlim([min(xx+sc)+3 max(xx+sc)])
    ylim([0 10])
    xticks([1359:2:1363])
    set(gca,'FontSize',fSize)
    
    
end
if pltChain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Potential parameters of interest
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ffFac=squeeze(outDat.ff(:,4,:));
    ffSpt=squeeze(outDat.ff(:,9,:));
    ffOsf=squeeze(outDat.ff(:,6,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x=ffFac(:,4); %Select parameter draws to plot
    figure2('Position',[10 10 1000 1000])
    plot(x)
    xlabel('Iteration')
    set(gca,'FontSize',fSize)
    
end
if pltTimeseries
    fSize=18;
    c95=[.65 .65 .65];
    c80=[.4 .4 .4];
    %Convert the proxies back to native units
    facMu=offsets(4); %Bremen offset
    facStd=nanstd(valAll(:,4));
    sptMu=offsets(17); %SILSO offset
    sptStd=nanstd(valAll(:,17));
    osfMu=offsets(8); %Krivova OMF offset
    osfStd=nanstd(valAll(:,8));
    
    figure2('Position',[10 10 1600 900])
    subplot('position',[.05 .69 .26 .3]) 
    facAll=outDat.facAll;
    facX=prctile(facAll',[2.5 10 50 90 97.5])'.*facStd+facMu;
    x2 = [dateM(outDat.dateFac)', fliplr(dateM(outDat.dateFac)')];
    fill(x2,[facX(:,1)',fliplr(facX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[facX(:,2)',fliplr(facX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    plot(dateM(outDat.dateFac),facX(:,3),'r')
    ylim([-1.5 3.2])
    xlabel('Year')
    ylabel('Facular Index')
    set(gca,'FontSize',fSize)
    
    subplot('position',[.375 .69 .26 .3]) 
    sptAll=outDat.sptAll;
    sptX=prctile(sptAll',[2.5 10 50 90 97.5])'.*sptStd+sptMu;
    x2 = [dateM(outDat.dateSpt)', fliplr(dateM(outDat.dateSpt)')];
    fill(x2,[sptX(:,1)',fliplr(sptX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[sptX(:,2)',fliplr(sptX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    plot(dateM(outDat.dateSpt),sptX(:,3),'r')
    ylim([-1.5 3.2])
    xlabel('Year')
    ylabel('Sunspot Index')
    set(gca,'FontSize',fSize)
    
    subplot('position',[.70 .69 .26 .3]) 
    osfAll=outDat.osfAll;
    osfX=prctile(osfAll',[2.5 10 50 90 97.5])'.*osfStd+osfMu;
    x2 = [dateM(outDat.dateOsf)', fliplr(dateM(outDat.dateOsf)')];
    fill(x2,[osfX(:,1)',fliplr(osfX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[osfX(:,2)',fliplr(osfX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    plot(dateM(outDat.dateOsf),osfX(:,3),'r')
    xlabel('Year')
    ylabel('Open Solar Flux Index')
    set(gca,'FontSize',fSize)
    
    subplot('position',[.05 .07 .91 .56]) 
    tsiAll=xAll;
    sorceI=oM(:,18);
    tsiAll=sorceanomaly(tsiAll,sorceI);
    tsiX=prctile(tsiAll',[2.5 10 50 90 97.5])';
    x2 = [dateM', fliplr(dateM')];
    fill(x2,[tsiX(:,1)',fliplr(tsiX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[tsiX(:,2)',fliplr(tsiX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    plot(dateM,tsiX(:,3),'r')
    legend('95% CI', '80% CI', 'BTSI Median Estimate','Location','NorthWest')
    legend boxoff
    ylim([-10 1.8])
    xlabel('Year')
    ylabel('TSI Anomaly (W/m^{2})')
    set(gca,'FontSize',fSize)
end
if obsContTSI
    satindex=logical([1 1 1 0 0 1 0 0 0 0 0 1 1 0 0 0 0 1 0 1 1 0 0 1]);
    obsUsed=satindex;
    obsUsed(4)=true; %Turn on Bremen Mg-II
    obsUsed(17)=true; %Turn on SILSO spots
    obsUsed(8)=true; %Turn on Krivova OMF using ISN
    colLabels=colLabels(obsUsed);
    %Plot the relative contribution of each observer to the estimate of TSI
    smoothWindow = 24;
    %First, try just plotting an example
    conChain=outDat.contributionChain;
    cn=squeeze(mean(abs(conChain),1));
    cn=cn./nansum(cn,2);
    for ii=1:13
        cn(:,ii)=movmean(cn(:,ii),smoothWindow,'omitnan');
    end
    
    %Plot
    figure2('Position',[10 10 1600 700])
    plot(dateM,cn(:,[4 6 9]),'--','LineWidth',3)
    hold on
    plot(dateM,cn(:,[1 2 3 5 7 8 10 11 12 13]),'LineWidth',2)
    legend(colLabels([4 6 9 1 2 3 5 7 8 10 11 12 13]),'Location','NorthWest')
    legend boxoff
    xlabel('Year')
    ylabel('Fractional Contribution')
    set(gca,'FontSize',fSize)
    saveas(gcf,'plots/obscontribution_22_12_01.png')
end
if obsContFac
    satindex=logical([1 1 1 0 0 1 0 0 0 0 0 1 1 0 0 0 0 1 0 1 1 0 0 1]);
    obsUsed=satindex;
    obsUsed(4)=true; %Turn on Bremen Mg-II
    obsUsed(18)=true; %Turn on SILSO spots
    obsUsed(9)=true; %Turn on Krivova OMF using ISN
    colLabels=colLabels(obsUsed);
    %Plot the relative contribution of each observer to the estimate of TSI
    smoothWindow = 24;
    %First, try just plotting an example
    conChain=outDat.cChainFac;
    cn=squeeze(mean(abs(conChain),1));
    cn=cn./nansum(cn,2);
    for ii=1:size(cn,2)
        cn(:,ii)=movmean(cn(:,ii),smoothWindow,'omitnan');
    end
    
    %Plot
    figure2('Position',[10 10 1600 700])
    plot(dateM(outDat.dateFac),cn,'LineWidth',2)
    legend(outDat.colsFac,'Location','NorthWest')
    legend boxoff
    xlabel('Year')
    ylabel('Fractional Contribution')
    set(gca,'FontSize',fSize)
    saveas(gcf,'plots/obscontributionfac_22_12_05.png')
end
if obsContSpt
    satindex=logical([1 1 1 0 0 1 0 0 0 0 0 1 1 0 0 0 0 1 0 1 1 0 0 1]);
    obsUsed=satindex;
    obsUsed(4)=true; %Turn on Bremen Mg-II
    obsUsed(18)=true; %Turn on SILSO spots
    obsUsed(9)=true; %Turn on Krivova OMF using ISN
    colLabels=colLabels(obsUsed);
    %Plot the relative contribution of each observer to the estimate of TSI
    smoothWindow = 24;
    %First, try just plotting an example
    conChain=outDat.cChainSpt;
    cn=squeeze(mean(abs(conChain),1));
    cn=cn./nansum(cn,2);
    for ii=1:size(cn,2)
        cn(:,ii)=movmean(cn(:,ii),smoothWindow,'omitnan');
    end
    
    %Plot
    figure2('Position',[10 10 1600 700])
    plot(dateM(outDat.dateSpt),cn,'LineWidth',2)
    legend(outDat.colsSpt,'Location','NorthWest')
    legend boxoff
    xlabel('Year')
    ylabel('Fractional Contribution')
    set(gca,'FontSize',fSize)
    saveas(gcf,'plots/obscontributionspt_22_12_05.png')
end
if obsContOSF
    satindex=logical([1 1 1 0 0 1 0 0 0 0 0 1 1 0 0 0 0 1 0 1 1 0 0 1]);
    obsUsed=satindex;
    obsUsed(4)=true; %Turn on Bremen Mg-II
    obsUsed(18)=true; %Turn on SILSO spots
    obsUsed(9)=true; %Turn on Krivova OMF using ISN
    colLabels=colLabels(obsUsed);
    %Plot the relative contribution of each observer to the estimate of TSI
    smoothWindow = 24;
    %First, try just plotting an example
    conChain=outDat.cChainOsf;
    cn=squeeze(mean(abs(conChain),1));
    cn=cn./nansum(cn,2);
    for ii=1:size(cn,2)
        cn(:,ii)=movmean(cn(:,ii),smoothWindow,'omitnan');
    end
    
    %Plot
    figure2('Position',[10 10 1600 700])
    plot(dateM(outDat.dateOsf),cn,'LineWidth',2)
    legend(outDat.colsOsf,'Location','NorthWest')
    legend boxoff
    xlabel('Year')
    ylabel('Fractional Contribution')
    set(gca,'FontSize',fSize)
    saveas(gcf,'plots/obscontributionosf_22_12_05.png')
end
if residualFac
    facX=mean(outDat.facAll,2);
    colsFac=outDat.colsFac;
    dateFac=dateM(outDat.dateFac);
    
    %Plot the facular index and all contributors over time
    figure
    plot(dateFac,facX,'LineWidth',3)
    hold on
    plot(dateFac,fac.data,'.')
    legendtxt=["fac index";colsFac];
    legend(legendtxt)
    xlabel('Year')
    ylabel('Normalized anomaly')
    set(gca,'FontSize',fSize);

    figure
    for ii=1:length(colsFac)
        pred=[ones(length(dateFac),1) fac.data(:,ii)]*fac.hload(ii,:)';
        plot(dateFac,pred-facX,'.')
        legend(colsFac(ii))
        xlabel('Year')
        ylabel('Residual')
        set(gca,'FontSize',fSize);
        pause
    end
end
if residualSpt
    sptX=mean(outDat.sptAll,2);
    colsSpt=outDat.colsSpt;
    dateSpt=dateM(outDat.dateSpt);
    
    %Plot the facular index and all contributors over time
    figure
    plot(dateSpt,sptX,'LineWidth',3)
    hold on
    plot(dateSpt,spt.data,'.')
    legendtxt=["spot index";colsSpt];
    legend(legendtxt)
    xlabel('Year')
    ylabel('Normalized anomaly')
    set(gca,'FontSize',fSize);

    figure
    for ii=1:length(colsSpt)
        pred=[ones(length(dateSpt),1) spt.data(:,ii)]*spt.hload(ii,:)';
        plot(dateSpt,pred-sptX,'.')
        legend(colsSpt(ii))
        xlabel('Year')
        ylabel('Residual')
        set(gca,'FontSize',fSize);
        pause
    end
end
if residualOsf
    osfX=mean(outDat.osfAll,2);
    colsOsf=outDat.colsOsf;
    dateOsf=dateM(outDat.dateOsf);
    
    %Plot the facular index and all contributors over time
    figure
    plot(dateOsf,osfX,'LineWidth',3)
    hold on
    plot(dateOsf,osf.data,'.')
    legendtxt=["OSF index";colsOsf];
    legend(legendtxt)
    xlabel('Year')
    ylabel('Normalized anomaly')
    set(gca,'FontSize',fSize);

    figure
    for ii=1:length(colsOsf)
        pred=[ones(length(dateOsf),1) osf.data(:,ii)]*osf.hload(ii,:)';
        plot(dateOsf,pred-osfX,'.')
        legend(colsOsf(ii))
        xlabel('Year')
        ylabel('Residual')
        set(gca,'FontSize',fSize);
        pause
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS FOR AGU 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if modernTSI
fSize=20;
%Order for proxies is sunspots first then MgII
%Order for sats is chronological from first observation: HF, ACRIM1, ERBE, ACRIM2,
%VIRGO/SOHO, ACRIM3, TIM/SORCE, PREMOS/PICARD, TCTE, TSIS-1
%Order as proxies followed by satellites in chronological order
load ar2_22_11_03.mat 
load(outDat.obsmatrix); %From makeobsmatrix.m
% lI=[7;1;5;2;13;3;10;8;11;12];
lI=[8;4;6;1;5;2;12;3;9;7;10;11];
%Reorder arrays used in this plotting script

A=A(lI,:,:);
colLabels=colLabels(lI);
offsets=offsets(lI);
oM=oM(:,lI);
sigY=sigY(lI,:);
t=t(:,lI);
valM=valM(:,lI);

    showTrend = 1;
    smoothWindow = 6; %set smoothing (months)
    %Create x-axis points for cycle demarcation
    c21=datejd([dateS.cycles(1,:), fliplr(dateS.cycles(1,:))]);
    c23=datejd([dateS.cycles(3,:), fliplr(dateS.cycles(3,:))]);
    c25=datejd([dateS.cycles(5,:), fliplr(dateS.cycles(5,:))]);
    
    figure2('Position',[10 10 1000 1200])
    
    subplot('position',[.09 .6 .85 .37]) %Plot of satellite observations
    fill(c21,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c23,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c25,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
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
        hold on
        hh(ind) = plot(dateM,valM(:,ii)+offsets(ii),'.','MarkerSize',10,...
            'Color',c(ind,:));
        ind = ind + 1;
    end
    xlim([datetime(1980,2,1) datetime(2022,1,1)])
    ylim([1360 1374])
    legend(hh,colLabels(3:end),'NumColumns',2)
    legend boxoff
    ylabel('TSI (W/m^{2})')
    set(gca,'FontSize',fSize)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot of reconstructions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot('position',[.09 .09 .85 .475]) 
    fill(c21,[-1.2 -1.2 1.8 1.8],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c23,[-1.2 -1.2 1.8 1.8],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c25,[-1.2 -1.2 1.8 1.8],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    cColor = get(gca,'colororder');
    %Get CI of our estimate,plot
    tsix = prctile(xAll',[.5 5 50 95 99.5])';
    for iS = 1:size(tsix,2)
        tsix(:,iS) = smoothPH(tsix(:,iS),smoothWindow);
    end
    [~,tsiO] = meaninterval(dateM,tsix(:,3),1990,2010);
    tsix = tsix-tsiO;
    x2 = [dateM', fliplr(dateM')];
    fill(x2,[tsix(:,1)',fliplr(tsix(:,end)')], [.85 .85 .85],'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[tsix(:,2)',fliplr(tsix(:,4)')], [.75 .75 .75],'FaceAlpha',0.5,...
        'LineStyle','none');
    
    %     %Get satellite only reconstruction
    load ar2_22_11_04_satnodrift.mat
    xms = mean(xAll,2);
    xms = smoothPH(xms,smoothWindow);
    [~,xmsO] = meaninterval(dateM,xms,1990,2010);
    xms = xms-xmsO; 
    
%     %Get MLR model reconstruction
    xmmlr=mean(xMLR,2);
    xmmlr=smoothPH(xmmlr,smoothWindow);
    [~,xmmlrO] = meaninterval(dateM,xmmlr,1990,2010);
    xmmlr = xmmlr-xmmlrO;
    
    %Plot other TSI reconstructions
    ind = 1;
    load oTSI_22_10_27.mat
    %Plot PMOD
    hold on
    tsiA = smoothPH(oTSI(7).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(7).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(7).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(7,:));
    legendtxt(ind) = string(oTSI(7).product);
    ind = ind + 1;
    
    %Plot SOLID
    hold on
    tsiA = smoothPH(oTSI(9).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(9).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(9).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(8,:));
    legendtxt(ind) = "Comm.-Consensus";
    ind = ind + 1;
    
    %Plot NRLTSI2
    hold on
    tsiA = smoothPH(oTSI(4).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(4).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(4).mthdatetime,tsiA-tsiAO,'LineWidth',2.5,...
        'Color',c(5,:));
    legendtxt(ind) = string(oTSI(4).product);
    ind = ind + 1;
    

    %Plot SATIRE-S during whole interval
    hold on
    tsiA = smoothPH(oTSI(5).mthtsi,smoothWindow);
    [~,tsiAO] = meaninterval(oTSI(5).mthdatetime,tsiA,1990,2010);
    plot(oTSI(5).mthdatetime,tsiA-tsiAO,'--','LineWidth',2.5,...
        'Color',c(4,:));
    legendtxt(ind)=string(oTSI(5).product);
    hold on %Plot SATIRE-S again during period of high-quality MDI/HMI data
    dateMDIHMI=oTSI(5).mthdatetime(oTSI(5).mthdatetime>datetime(1999,02,02));
    h(ind) = plot(dateMDIHMI,tsiA(oTSI(5).mthdatetime>datetime(1999,02,02))-tsiAO,'LineWidth',2.5,...
        'Color',c(4,:));
    
    
    
    ind = ind+1; %Index our contribution
    hold on
    h(ind) = plot(dateM,tsix(:,3),'LineWidth',3.5,'Color','k');
    legendtxt(ind) = "BTSI";
    ind = ind+1; %Index our contribution
    hold on
    h(ind)=plot(dateM,xms,'LineWidth',2.5,'Color',c(1,:));
    legendtxt(ind)="BTSI sat. only";
    ind=ind+1;
    hold on
    h(ind)=plot(dateM,xmmlr,'LineWidth',2.5,'Color',c(2,:));
    legendtxt(ind)="BTSI proxy model";

    

    legend(h,legendtxt,'NumColumns',2)
    legend boxoff
    set(gca,'FontSize',fSize)
    xlabel('Year')
    ylabel('TSI anomaly from 1990-2010 mean (W/m^{2})')
    xlim([datetime(1980,2,1) datetime(2022,1,1)])
    ylim([-0.9 1.25])
    text(datejd(dateS.cycles(1,1))+years(4.5),-0.8,'Cycle 21','FontSize',14)
    text(datejd(dateS.cycles(2,1))+years(3.25),-0.8,'Cycle 22','FontSize',14)
    text(datejd(dateS.cycles(3,1))+years(4.75),-0.8,'Cycle 23','FontSize',14)
    text(datejd(dateS.cycles(4,1))+years(4.25),-0.8,'Cycle 24','FontSize',14)
end
if pltAllTimeseries
    fSize=17;
    load oTSI_long_22_12_07.mat
    c95=[.65 .65 .65];
    c80=[.4 .4 .4];
    %Convert the proxies back to native units
    facMu=offsets(4); %Bremen offset
    facStd=nanstd(valAll(:,4));
    sptMu=offsets(17); %SILSO offset
    sptStd=nanstd(valAll(:,17));
    osfMu=offsets(8); %Krivova OMF offset
    osfStd=nanstd(valAll(:,8));
    
    figure2('Position',[10 10 1600 900])
    subplot('position',[.05 .69 .26 .3]) 
    facAll=outDat.facAll;
    facX=prctile(facAll',[2.5 10 50 90 97.5])'.*facStd+facMu;
    x2 = [dateM(outDat.dateFac)', fliplr(dateM(outDat.dateFac)')];
    fill(x2,[facX(:,1)',fliplr(facX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[facX(:,2)',fliplr(facX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    plot(dateM(outDat.dateFac),facX(:,3),'r')
    ylim([0 0.8])
    xlabel('Year')
    ylabel('Facular Index')
    set(gca,'FontSize',fSize)
    
    subplot('position',[.375 .69 .26 .3]) 
    sptAll=outDat.sptAll;
    sptX=prctile(sptAll',[2.5 10 50 90 97.5])'.*sptStd + sptMu;
    x2 = [dateM(outDat.dateSpt)', fliplr(dateM(outDat.dateSpt)')];
    fill(x2,[sptX(:,1)',fliplr(sptX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[sptX(:,2)',fliplr(sptX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    plot(dateM(outDat.dateSpt),sptX(:,3),'r')
    ylim([0 300])
    xlabel('Year')
    ylabel('Sunspot Number')
    set(gca,'FontSize',fSize)
    
    subplot('position',[.70 .69 .26 .3]) 
    osfAll=outDat.osfAll;
    osfX=prctile(osfAll',[2.5 10 50 90 97.5])'.*osfStd+osfMu;
    x2 = [dateM(outDat.dateOsf)', fliplr(dateM(outDat.dateOsf)')];
    fill(x2,[osfX(:,1)',fliplr(osfX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[osfX(:,2)',fliplr(osfX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    plot(dateM(outDat.dateOsf),osfX(:,3),'r')
    xlabel('Year')
    ylabel('Open Solar Flux (10^{14} Wb)')
    ylim([0 10])
    set(gca,'FontSize',fSize)
    
    subplot('position',[.05 .07 .91 .56]) 
    tsiAll=xAll;
    sorceI=oM(:,18);
    tsiAll=sorceanomaly(tsiAll,sorceI)+nanmean(valM(:,18));
    tsiX=prctile(tsiAll',[2.5 10 50 90 97.5])';
    x2 = [dateM', fliplr(dateM')];
    fill(x2,[tsiX(:,1)',fliplr(tsiX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[tsiX(:,2)',fliplr(tsiX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    legendtxt=["95% CI"; "80% CI"];
    for ii=1:length(oTSI)
        plot(oTSI(ii).yrdatetime,oTSI(ii).yrtsi,'Color',c(ii,:),'LineWidth',2)
        hold on
        legendtxt=[legendtxt; string(oTSI(ii).product)];
    end
    hold on
    plot(dateM,tsiX(:,3),'r','LineWidth',3)
    legendtxt=[legendtxt;"BTSI Median Estimate"];
    legend(legendtxt,'Location','SouthEast')
    legend boxoff
    xlabel('Year')
    ylabel('TSI (W/m^{2})')
    ylim([1357.5 1362.75])
    yticks([1358:1362])
    set(gca,'FontSize',fSize)
end
if pltOtherTimeseries
    fSize=20;
    load oTSI_long_22_12_07.mat
    legendtxt=[];
    figure2('Position',[10 10 1600 900])
    for ii=1:length(oTSI)
        plot(oTSI(ii).yrdatetime,oTSI(ii).yrtsi,'Color',c(ii,:),'LineWidth',2)
        hold on
        legendtxt=[legendtxt; string(oTSI(ii).product)];
    end
    legend(legendtxt,'Location','NorthWest')
    legend boxoff
    xlabel('Year')
    ylabel('TSI (W/m^{2})')
    set(gca,'FontSize',fSize)
    
end
if pltOSFCompare
    %Save OSF values
    osfindex=logical([0 0 0 0 0 0 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1 0]);
    valsOsf=valAll(:,osfindex)+offsets(osfindex);
    oMOsf=oM(:,osfindex);
    for ii=1:sum(osfindex)
        valsOsf(:,ii)=interp1(dateM(oMOsf(:,ii)),valsOsf(oMOsf(:,ii),ii),dateM);
    end
    dateOsf=dateM;
    colsOsf=colLabels(osfindex);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fSize=20;
    load ar2_22_11_03.mat
    load(outDat.obsmatrix); %From makeobsmatrix.m
    smoothWindow = 6; %set smoothing (months)
    %Create x-axis points for cycle demarcation
    c21=datejd([dateS.cycles(1,:), fliplr(dateS.cycles(1,:))]);
    c23=datejd([dateS.cycles(3,:), fliplr(dateS.cycles(3,:))]);
    c25=datejd([dateS.cycles(5,:), fliplr(dateS.cycles(5,:))]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Plot of reconstructions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure2('Position',[10 10 1600 900])
    subplot('position',[.09 .6 .8 .375]) 
    set(gca,'ytick',[])
    set(gca,'xtick',[])
    yyaxis right
    cColor = get(gca,'colororder');
    %Get CI of our estimate,plot
    tsix = prctile(xAll',[.5 5 50 95 99.5])';
    for iS = 1:size(tsix,2)
        tsix(:,iS) = smoothPH(tsix(:,iS),smoothWindow);
    end
    [~,tsiO] = meaninterval(dateM,tsix(:,3),1990,2010);
    tsix = tsix-tsiO;
    x2 = [dateM', fliplr(dateM')];

    
    %Get satellite only reconstruction
    load ar2_22_11_04_satnodrift.mat
    xms = mean(xAll,2);
    xms = smoothPH(xms,smoothWindow);
    [~,xmsO] = meaninterval(dateM,xms,1990,2010);
    xms = xms-xmsO; 
    
    %Get MLR model reconstruction
    xmmlr=mean(xMLR,2);
    xmmlr=smoothPH(xmmlr,smoothWindow);
    [~,xmmlrO] = meaninterval(dateM,xmmlr,1990,2010);
    xmmlr = xmmlr-xmmlrO;
    
    %Plot other TSI reconstructions
    ind = 1;
    load oTSI_22_10_27.mat
    %Plot PMOD
    hold on
    tsiA = smoothPH(oTSI(7).mthtsi,smoothWindow);
    dateSat=oTSI(7).mthdatetime>datejd(dateS.all(1));
    [~,tsiAO] = meaninterval(oTSI(7).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(7).mthdatetime(dateSat),tsiA(dateSat)-tsiAO,'-','LineWidth',2.5,...
        'Color',c(7,:));
    legendtxt(ind) = string(oTSI(7).product);
    ind = ind + 1;
    
    %Plot SOLID
    hold on
    tsiA = smoothPH(oTSI(9).mthtsi,smoothWindow);
    dateSat=oTSI(9).mthdatetime>datejd(dateS.all(1));
    [~,tsiAO] = meaninterval(oTSI(9).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(9).mthdatetime(dateSat),tsiA(dateSat)-tsiAO,'-','LineWidth',2.5,...
        'Color',c(8,:));
    legendtxt(ind) = "Comm.-Consensus";
    ind = ind + 1;
    
    %Plot NRLTSI2
    hold on
    tsiA = smoothPH(oTSI(4).mthtsi,smoothWindow);
    dateSat=oTSI(4).mthdatetime>datejd(dateS.all(1));
    [~,tsiAO] = meaninterval(oTSI(4).mthdatetime,tsiA,1990,2010);
    h(ind) = plot(oTSI(4).mthdatetime(dateSat),tsiA(dateSat)-tsiAO,'-','LineWidth',2.5,...
        'Color',c(5,:));
    legendtxt(ind) = string(oTSI(4).product);
    ind = ind + 1;
    

    %Plot SATIRE-S during whole interval
    hold on
    tsiA = smoothPH(oTSI(5).mthtsi,smoothWindow);
    dateSat=oTSI(5).mthdatetime>datejd(dateS.all(1));
    [~,tsiAO] = meaninterval(oTSI(5).mthdatetime,tsiA,1990,2010);
    plot(oTSI(5).mthdatetime(dateSat),tsiA(dateSat)-tsiAO,'--','LineWidth',2.5,...
        'Color',c(4,:));
    legendtxt(ind)=string(oTSI(5).product);
    hold on %Plot SATIRE-S again during period of high-quality MDI/HMI data
    dateMDIHMI=oTSI(5).mthdatetime(oTSI(5).mthdatetime>datetime(1999,02,02));
    h(ind) = plot(dateMDIHMI,tsiA(oTSI(5).mthdatetime>datetime(1999,02,02))-tsiAO,'-','LineWidth',2.5,...
        'Color',c(4,:));
    
    
    
    ind = ind+1; %Index our contribution
    hold on
    h(ind) = plot(dateM,tsix(:,3),'-','LineWidth',3.5,'Color','k');
    legendtxt(ind) = "BTSI";
    ind = ind+1; %Index our contribution
    hold on
    h(ind)=plot(dateM,xms,'-','LineWidth',2.5,'Color',c(1,:));
    legendtxt(ind)="BTSI sat. only";
    ind=ind+1;
    hold on
    h(ind)=plot(dateM,xmmlr,'-','LineWidth',2.5,'Color',c(2,:));
    legendtxt(ind)="BTSI proxy model";

    

    legend(h,legendtxt,'Location','Northwest','NumColumns',2)
    legend boxoff
    set(gca,'FontSize',fSize)
    ylabel('TSI (W/m^{2})')
    
    xlim([datetime(1610,1,1) datetime(2022,1,1)])
    ylim([-.75 1.25])
    ax=gca;
    ax.YColor='k';
    subplot('position',[.09 .05 .8 .55]) 
    set(gca,'ytick',[])
    yyaxis right
    for ii=1:size(valsOsf,2)
        plot(dateOsf,valsOsf(:,ii),'-','LineWidth',2,'Color',c(ii,:))
        hold on
    end
    colsOsf=["Frost et al. 2022";"Krivova et al. 2021";"Lockwood et al. 2022";...
        "McCracken and Beer 2015"; "Usoskin et al. 2021"];
    legend(colsOsf,'Location','Northwest','NumColumns',2)
    legend boxoff
    xlabel('Year')
    ylabel('Open Solar Flux (10^{14} Wb)')
    xlim([datetime(1610,1,1) datetime(2022,1,1)])
    ylim([0 12])
    set(gca,'ytick',[0:2:10])
    ax=gca;
    ax.YColor='k';
    set(gca,'FontSize',fSize)
end
if squareSatObs
    fSize=20;
    load ar2_22_11_03.mat 
load(outDat.obsmatrix); %From makeobsmatrix.m
% lI=[7;1;5;2;13;3;10;8;11;12];
lI=[8;4;6;1;5;2;12;3;9;7;10;11];
%Reorder arrays used in this plotting script

A=A(lI,:,:);
colLabels=colLabels(lI);
offsets=offsets(lI);
oM=oM(:,lI);
sigY=sigY(lI,:);
t=t(:,lI);
valM=valM(:,lI);

    showTrend = 1;
    smoothWindow = 6; %set smoothing (months)
    %Create x-axis points for cycle demarcation
    c21=datejd([dateS.cycles(1,:), fliplr(dateS.cycles(1,:))]);
    c23=datejd([dateS.cycles(3,:), fliplr(dateS.cycles(3,:))]);
    c25=datejd([dateS.cycles(5,:), fliplr(dateS.cycles(5,:))]);
    
    figure2('Position',[10 10 1000 1200])
    fill(c21,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c23,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
        0.4,'LineStyle','none');
    hold on
    fill(c25,[1360 1360 1374 1374],[.96 .96 .863],'FaceAlpha',...
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
        hold on
        hh(ind) = plot(dateM,valM(:,ii)+offsets(ii),'.','MarkerSize',10,...
            'Color',c(ind,:));
        ind = ind + 1;
    end
    xlim([datetime(1980,2,1) datetime(2022,1,1)])
    ylim([1360 1374])
    legend(hh,colLabels(3:end),'NumColumns',2)
    legend boxoff
    ylabel('TSI (W/m^{2})')
    set(gca,'FontSize',fSize)
end
if squareOsf
    fSize=20;
    c95=[.65 .65 .65];
    c80=[.4 .4 .4];
    osfMu=offsets(8); %Krivova OMF offset
    osfStd=nanstd(valAll(:,8));
    osfAll=outDat.osfAll.*osfStd+osfMu;
    osfX=prctile(osfAll',[2.5 10 50 90 97.5])';
    x2 = [dateM(outDat.dateOsf)', fliplr(dateM(outDat.dateOsf)')];
    figure
    fill(x2,[osfX(:,1)',fliplr(osfX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[osfX(:,2)',fliplr(osfX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    
    hold on
    osfindex=logical([0 0 0 0 0 0 1 1 0 1 1 0 0 0 0 0 0 0 0 0 0 0 1 0]);
    valsOsf=valAll(:,osfindex)+osfMu;
    oMOsf=oM(:,osfindex);
    for ii=1:sum(osfindex)
        valsOsf(:,ii)=interp1(dateM(oMOsf(:,ii)),valsOsf(oMOsf(:,ii),ii),dateM);
        h(ii)=plot(dateM,valsOsf(:,ii),'LineWidth',2);
        hold on
    end
    dateOsf=dateM;
    colsOsf=["Frost et al. 2022";"Krivova et al. 2021";"Lockwood et al. 2022";...
        "McCracken and Beer 2015"; "Usoskin et al. 2021"];
    
    
    hold on
    h(ii+1)=plot(dateM(outDat.dateOsf),osfX(:,3),'r','LineWidth',2);
    colsOsf=[colsOsf;"OSF BTSI estimate"];
    xlabel('Year')
    ylabel('Open Solar Flux (10^{14} Wb)')
    legend(h,colsOsf,'Location','NorthWest')
    legend boxoff
    ylim([0 14])
    set(gca,'FontSize',fSize)
end
if squareNonlinear
    c95=[.65 .65 .65];
    c80=[.4 .4 .4];
    xbar=2.2;
    YEO21=1358.7; %Solar constant from Yeo et al. 2021
    sc=YEO21;
    regX=(sc:0.01:sc+5)'; %input for linear regression of indices
    %Get mean tsi output
    x=mean(xAll,2);
    cmpStYr=1978; %cmpStYr=outDat.cmpStYr
    xPlot=x(dateM.Year>=cmpStYr)+sc;
    %Get mean output for each process after the cutoff date
    %Get OSF contributions
    infI=dateM(outDat.dateOsf).Year>=cmpStYr;
    osfMu=offsets(8); %Krivova OMF offset
    osfStd=nanstd(valAll(:,8));
    yOsf=outDat.osfAll(infI,:);
    yOsfMean=mean(yOsf,2)*osfStd+osfMu;
    ffOsf=squeeze(outDat.ff(:,6,:));
    
    
    
    figure2('Position',[10 10 1000 1000])
    plot(xPlot,yOsfMean,'.')
    hold on
    xx=0:0.01:5;
    glm= @(b,x) b(4)+(b(1)-b(4))./(1+exp(b(2)*(x-(xbar+b(3)))));
    for ii=1:size(ffOsf,1)
        aOsf(:,ii)=glm(ffOsf(ii,:),xx);
    end
    aOsfPrctile=prctile(aOsf',[2.5 10 50 90 97.5])'.*osfStd+osfMu;
    x2 = [xx, fliplr(xx)]+sc;
    fill(x2,[aOsfPrctile(:,1)',fliplr(aOsfPrctile(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[aOsfPrctile(:,2)',fliplr(aOsfPrctile(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    %Plot linear regression
    pred=[ones(length(xPlot),1) xPlot];
    b=regress(yOsfMean,pred);
    plotPred=[ones(length(regX),1) regX];
    plot(regX,(plotPred*b),'--','Color','k','LineWidth',2)
    hold on
    plot(xx+sc,aOsfPrctile(:,3),'Color','r','LineWidth',2)
    xlabel('TSI (W/m^{2})')
    ylabel('Open Solar Flux (10^{14} Wb)')
    xlim([sc sc+5])
    ylim([min(aOsfPrctile(:)) max(aOsfPrctile(:))])
    xticks([1359:2:1363])
    set(gca,'FontSize',fSize)
    
end
if squareTSITimeSeries
    fSize=19;
    load oTSI_long_22_12_07.mat
    c95=[.65 .65 .65];
    c80=[.4 .4 .4];
    figure2('Position',[10 10 1600 900])
    tsiAll=xAll;
    sorceI=oM(:,18);
    tsiAll=sorceanomaly(tsiAll,sorceI)+nanmean(valM(:,18));
    tsiX=prctile(tsiAll',[2.5 10 50 90 97.5])';
    x2 = [dateM', fliplr(dateM')];
    fill(x2,[tsiX(:,1)',fliplr(tsiX(:,end)')], c95,'FaceAlpha',...
        0.5,'LineStyle','none');
    hold on
    fill(x2,[tsiX(:,2)',fliplr(tsiX(:,4)')], c80,'FaceAlpha',0.5,...
        'LineStyle','none');
    hold on
    legendtxt=["95% CI"; "80% CI"];
    for ii=1:length(oTSI)
        plot(oTSI(ii).yrdatetime,oTSI(ii).yrtsi,'Color',c(ii,:),'LineWidth',2)
        hold on
        legendtxt=[legendtxt; string(oTSI(ii).product)];
    end
    hold on
    plot(dateM,tsiX(:,3),'r','LineWidth',3)
    legendtxt=[legendtxt;"BTSI Median Estimate"];
    legend(legendtxt,'Location','NorthWest')
    legend boxoff
    xlabel('Year')
    ylabel('TSI (W/m^{2})')
    ylim([1357.5 1362.75])
    yticks([1358:1362])
    set(gca,'FontSize',fSize)
end

function xout=sorceanomaly(x,sorceI)
for ii=1:size(x,2)
    xout(:,ii)=x(:,ii)-mean(x(sorceI,ii));
    xout(xout(:,ii)<-1.*mean(x(sorceI,ii)),ii)=-1.*mean(x(sorceI,ii));
end
end

    
    
    