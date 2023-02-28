function[H0, Hsig, T0, th0,X0,Xsig] = getpriors(NN,oindex,tindex,sindex,valM,oM,...
    satindex,colLabels,opts,priorpath)
thresh=24; %Threshold for using overlapping observers be it years or months
savePath='mat_files/obspriors_23_01_13.mat'; %Save path if priors are to be saved
%Generate prior hyperparameter estimates
if nargin<10   
    nObs = length(colLabels);
    %Make labels for which records fall under satellite vs proxy
    isProx=~satindex;
    proxInd=find(isProx);
    satInd=find(satindex);
    
    %Use NRLTSI relationship between sunspots, mg-II and TSI to set priors
    %and also the innovation uncertainty
    if isfield(opts,'NRLTSIprior') && opts.NRLTSIprior==true
        load oTSI_23_02_01.mat
        valT=valM;oMT=oM; %Temporarily save input values 
        load obs_23_02_01.mat
        dateS=getdates;
        dateR=dateS.acrimplusfive;
        tsiNRL=oTSI(4).mthtsi;
        dateNRL=oTSI(4).mthdatetime;
        x=tsiNRL(dateNRL>=datejd(dateR(1)) & dateNRL <= datejd(dateR(end)));
        x=x-nanmean(x);
        for ii=1:sum(isProx)
            overlap = oM(:,proxInd(ii));
            y=valM(overlap,proxInd(ii));
            pred=[ones(length(y),1) x(overlap)];
            b=regress(y,pred);
            res = y - pred*b;
            obsPrior(proxInd(ii)).type = "process";
            obsPrior(proxInd(ii)).name = colLabels(proxInd(ii));
            obsPrior(proxInd(ii)).std = nanstd(res);
            obsPrior(proxInd(ii)).b = b(1);
            obsPrior(proxInd(ii)).m = b(2);
            %One sigma uncertainty in offset (note this is large due to
            %lack of centering of proxy records)
            obsPrior(proxInd(ii)).bsig = b(2).*.5; 
            %One sigma uncertainty in scaling
            obsPrior(proxInd(ii)).msig = b(2).*.25; 
        end
        
        %Get residual errors of autoregressive AR(2) model
        YCycleAll=x-movmean(x,12.*11,'omitnan');
        YCycleAll=YCycleAll-min(YCycleAll); %Set Y-intercept to minimum TSI value
        Y=x;
        X=[];
        for iL=1:opts.lags %Create columns for lagged estimates of x
            X=[X lag0(x,iL)];
        end
        X=[ones(size(Y,1),1) X];%To ensure a stable linear regression, estimate mean of TSI
        Y=Y(2:end,:);
        X=X(2:end,:);
        
        M=inv(X'*X)*(X'*Y);
        errorsx=Y-X*M;
        
        %Estimate regression parameters of noise as a fn. of solar
        %cycle magnitude
        [X0,bint]=regress(abs(errorsx),[ones(length(errorsx),1) YCycleAll(2:end)]);
        Xsig=(bint(:,2)-bint(:,1))./4; %Convert from 95\% CI to one sigma
        
        
        valM=valT;oM=oMT; %Restore original values
        
    %Use raw observations to set priors, RISK OF REGRESSION DILUTION
    else
        %Make a set of priors for the proxy observations using all satellites
        for ii = 1:sum(isProx)
            ind=proxInd(ii);
            iS = 1;
            vProx = [];
            bProx = [];
            bIntProx = [];
            threshUsed=thresh;
            while iS <= nObs
                if ~any(iS==proxInd) %only operate on satellites
                    overlap = and(oM(:,iS),oM(:,ind));
                    if sum(overlap) >= threshUsed %5 year cutoff for comparison to be made
                        %Revised 9/8/21 to be in native units
                        pred=[ones(sum(overlap),1) valM(overlap,iS)-nanmean(valM(overlap,iS))];
                        [b,bint] = regress(valM(overlap,ind),pred);
                        res = valM(overlap,ind) - pred*b;
                        vProx = [vProx; nanvar(res,1)];
                        bProx = [bProx; b'];
                        bIntProx = [bIntProx; bint(2,:)];
                    end
                end
                iS = iS + 1;
            end
            obsPrior(proxInd(ii)).type = "process";
            obsPrior(proxInd(ii)).name = colLabels(ind);
            obsPrior(proxInd(ii)).std = sqrt(nanmean(vProx));
            obsPrior(proxInd(ii)).b = mean(bProx(:,1));
            obsPrior(proxInd(ii)).m = mean(bProx(:,2));
            obsPrior(proxInd(ii)).bsig = std(bProx(:,1)); %One sigma uncertainty in offset
            obsPrior(proxInd(ii)).msig = std(bProx(:,2)); %One sigma uncertainty in scaling
        end
    end
    
    %Then, make a set of comparisons between all the primary observers
    for ii=1:nObs
        if any(ii==satInd) %Only make comparisons for primary observers
            vSat = NaN(nObs,1);
            iS = 1;
            while iS <= nObs
                if iS ~= ii && satindex(iS)
                    overlap = logical(oM(:,iS).*oM(:,ii));
                    if opts.normalize
                        s1 = normPH(valM(overlap,ii));
                        s2 = normPH(valM(overlap,iS)); 
                    else
                        s1 = valM(overlap,ii) - mean(valM(overlap,ii));
                        s2 = valM(overlap,iS) - mean(valM(overlap,iS));
                    end
                    dSat = s2-s1;
                    vSat(iS) = nanvar(dSat,1);
                end
                iS = iS + 1;
            end
            %Include the data collected into a structure object
            if isfield(opts,'type')
                obsPrior(ii).type = opts.type;
            else
                obsPrior(ii).type = "sat";
            end
            obsPrior(ii).name = colLabels(ii);
            obsPrior(ii).std = sqrt(nanmean(vSat))./sqrt(2); %Correct for satellite noise coming from two observers
        end
    end
    colsPrior=colLabels;
    if isfield(opts,'save')
        save(savePath,'obsPrior','colsPrior')
    end
else
    load(priorpath);
end

priFrac = 0.25; %Fraction of the 'observations' for sigX, sigY to be from prior
offSig = 5; %prior satellite variance for a (offset)
mSig = 0.25; %prior variance for c (satellite drift)

H0=zeros(NN,3);H0(:,2) = 1; %prior for scaling H of form [a_i b_i c_i]
Hsig=1E-12.*ones(NN,3);
Hsig=Hsig.*H0(:,2);
Hsig(:,1)=Hsig(:,1)+(offSig.*oindex');%prior sigma for offset
Hsig(:,2)=Hsig(:,2)+(1000.*sindex');%prior sigma for scaling
Hsig(:,3)=Hsig(:,3)+(mSig.*tindex'); %prior sigma for drift
%Make changes for proxies if linearity assumed (first row sunspots, second row mg)
proxI=find(~satindex);%indices for proxies
for iP=1:length(proxI)
    Hsig(proxI(iP),1)=obsPrior(proxI(iP)).bsig.^2; %Var uncertainty for proxy offset
    Hsig(proxI(iP),2)=obsPrior(proxI(iP)).msig.^2; %var uncertainty for proxy scaling
    H0(proxI(iP),2)=obsPrior(proxI(iP)).m; %Proxy Scaling uncertainty prior
    Hsig(proxI(iP),3)=1E-12.*H0(proxI(iP),2).*ones(length(proxI(iP)),1); %Fix slope
end

if isfield(opts,'HsigScale')
    Hsig=Hsig.*opts.HsigScale;
end
%Next, draw a range of initial hyperparameters for the prior noise variance
iObs=sum(oM,1);
T0=ceil(iObs.*(priFrac./(1-priFrac))); %Assumes std from intercomparison
for ii = 1:NN %Prior for proxy noise
    th0(ii) = (T0(ii)-1).*obsPrior(ii).std.^2;
end

end

