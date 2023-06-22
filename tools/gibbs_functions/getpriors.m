function [H0, Hsig, T0, th0, eta0, Xsig] = getpriors(NN, oindex, tindex, pindex, valM, ...
    oM, satindex, colLabels, opts, priorpath, satDateM)
% GETPRIORS - generate prior hyperparameter estimates for a Bayesian linear regression model
% [H0, Hsig, T0, th0, X0, Xsig] = getpriors(NN, oindex, tindex, pindex, valM, oM, satindex, 
%  colLabels, opts, priorpath, satDateM)
%
% INPUTS:
% - NN: scalar integer representing the number of columns in the data matrix
% - oindex: vector of length NN representing whether each observer requires offsets (1) or not (0)
% - tindex: vector of length NN representing whether each observer requires a time variable (1) or not (0)
% - pindex: vector of length NN representing whether each observer is a proxy variable (1) or not (0)
% - valM: data matrix of size N x NN, where N is the number of observations
% - oM: valid observation matrix of size N x NN, where N is the number of
% observations
% - satindex: vector of length NN representing whether each column is a satellite (1) or not (0)
% - colLabels: string array of length NN containing column labels
% - opts: struct containing options for the function, including:
%   - priorthresh: threshold for using overlapping observers be it years or months (default: 24)
%   - NRLTSIprior: whether to use NRLTSI estimates of TSI to set priors for proxies (default: false)
%   - normalize: whether to normalize observations before comparing them (default: false)
%   - type: type of observation (default: 'sat')
%   - save: whether to save priors to a file (default: false)
%   - HsigScale: scale factor for prior sigma for H (default: 1)
% - priorpath: path to file containing priors (default: empty)
% - satDateM: vector containing satellite dates
%
% OUTPUTS:
% - H0: prior for scaling H of form [a_i b_i c_i]
% - Hsig: prior covariance hyperparameters for scaling H
% - T0: vector of initial hyperparameters for the prior noise variance
% - th0: prior for proxy noise variance
% - eta0: initial guess for eta
% - Xsig: prior sigma for eta

if ~isfield(opts,'priorthresh')
    opts.priorthresh=24; %Threshold for using overlapping observers be it years or months
end

%Generate prior hyperparameter estimates
if ~exist('priorpath', 'var') || isempty(priorpath)
    nObs = size(valM,2);
    %Make labels for which records fall under satellite vs proxy
    proxInd=find(pindex);
    satInd=find(satindex);
    
    %Use NRLTSI estimates of TSI to set priors for proxies
    %and also the innovation uncertainty
    if isfield(opts,'NRLTSIprior') && opts.NRLTSIprior==true
        % If the NRLTSIprior option is set to true in the opts struct, use NRLTSI data
        % to set priors for proxies
        
        % Load the NRLTSI data
        load oTSI_23_02_01.mat
        tsiNRL=oTSI(4).mthtsi;
        dateNRL=oTSI(4).mthdatetime;
        
        % Filter the NRLTSI data for the satellite period
        x=tsiNRL(dateNRL>=satDateM(1)-caldays(5) & dateNRL <= satDateM(end)+caldays(5));
        x=x-nanmean(x);
        % Loop over all the proxy observations
        for ii=1:sum(pindex)
            % Get the overlapping observations between the proxy and the
            % satellites
            overlap = oM(:,proxInd(ii));
            y=valM(overlap,proxInd(ii));
            
            % Perform linear regression between the overlapping NRLTSI data and the proxy
            pred=[ones(length(y),1) x(overlap)];
            b=regress(y,pred);
            
            % Calculate the residuals between the observed proxy data and the
            % regression estimate
            res = y - pred*b;
            
            % Store the prior information in a struct
            obsPrior(proxInd(ii)).type = "process";
            obsPrior(proxInd(ii)).name = colLabels(proxInd(ii));
            obsPrior(proxInd(ii)).std = nanstd(res);
            obsPrior(proxInd(ii)).b = b(1);
            obsPrior(proxInd(ii)).m = b(2);
            
            % One sigma uncertainty in offset (note this is large due to
            % lack of centering of proxy records)
            obsPrior(proxInd(ii)).bsig = b(2)*.5;
            % One sigma uncertainty in scaling
            obsPrior(proxInd(ii)).msig = b(2)*.25;
        end
        
        % Get residual errors of autoregressive AR model
        % Key onto the cyclic component of the NRLTSI data by subtracting a
        % moving average of window size 12*11 (11 years)
        YCycleAll=x-movmean(x,12.*11,'omitnan');
        % Set Y-intercept to minimum TSI value
        YCycleAll=YCycleAll-min(YCycleAll);
        Y=x;
        X=[];
        for iL=1:opts.lags % Create columns for lagged estimates of x
            X=[X lag0(x,iL)];
        end
        X=[ones(size(Y,1),1) X];
        
        Y=Y(2:end,:);
        X=X(2:end,:);
        
        % Calculate the regression coefficients of the residuals as a
        % function of solar cycle magnitude
        M=inv(X'*X)*(X'*Y);
        errorsx=Y-X*M;
        
        %Estimate regression parameters of noise as a fn. of solar
        %cycle magnitude
        if opts.magDependent
            [eta0,bint]=regress(errorsx.^2,[ones(length(errorsx),1) YCycleAll(2:end)]);
            if eta0(1)<(0.05.^2)
                eta0(1)=0.05.^2; %Set lower baseline of 0.05 W/m^2 innovation
            end
            Xsig=(bint(:,2)-bint(:,1))./4; %Convert from 95\% CI to one sigma
        else
            eta0=mean(errorsx.^2);
            Xsig=NaN;
        end
        
        
    %Use raw observations to set priors, RISK OF REGRESSION DILUTION
    else
        %If not using NRLTSI, set initial state vector to zero and estimate
        %uncertainty from intercomparisons
        for ii = 1:sum(~satindex)
            ind=proxInd(ii);
            iS = 1;
            vProx = [];
            bProx = [];
            bIntProx = [];
            threshUsed=priorthresh;
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
                    s1 = valM(overlap,ii) - mean(valM(overlap,ii));
                    s2 = valM(overlap,iS) - mean(valM(overlap,iS));
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
            %Correct for satellite noise coming from two observers
            obsPrior(ii).std = sqrt(nanmean(vSat))./sqrt(2); 
        end
    end
    colsPrior=colLabels;
    if isfield(opts,'save')
        %Create save path if we're saving the file
        currDate=datetime;
        savePath=['mat_files/obspriors_' num2str(mod(currDate.Year,100)) '_' ...
            num2str(currDate.Month) '_' num2str(currDate.Day) '.mat'];
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
if isfield(opts,'type')
    if strcmp(opts.type,'fac')
        Hsig(:,1)=Hsig(:,1)+(1.*oindex');%prior variance for offset
        Hsig(:,2)=Hsig(:,2)+(0.25.*pindex');%prior variance for scaling
    end
    if strcmp(opts.type,'spot')
        Hsig(:,1)=Hsig(:,1)+(1.*oindex');%prior variance for offset
        Hsig(:,2)=Hsig(:,2)+(0.25.*pindex');%prior variance for scaling
        for ii=1:length(obsPrior)
            if strcmp(obsPrior(ii).name,"SATIRE-S Sunspot Contribution")
                H0(ii,2)=-1;
            end
        end
    end
    if strcmp(opts.type,'osf')
        Hsig(:,1)=Hsig(:,1)+(1.*oindex');%prior variance for offset
        Hsig(:,2)=Hsig(:,2)+(0.25.*pindex');%prior variance for scaling
    end
else
    Hsig(:,1)=Hsig(:,1)+(offSig.*oindex');%prior variance for offset
    Hsig(:,2)=Hsig(:,2)+(1000.*pindex');%prior variance for scaling
    Hsig(:,3)=Hsig(:,3)+(mSig.*tindex'); %prior variance for drift
end
%Make changes for proxies if linearity assumed (first row sunspots, second row mg)
proxI=find(~satindex);%indices for proxies
for iP=1:length(proxI)
    Hsig(proxI(iP),1)=obsPrior(proxI(iP)).bsig.^2; %Var uncertainty for proxy offset
    Hsig(proxI(iP),2)=obsPrior(proxI(iP)).msig.^2; %var uncertainty for proxy scaling
    Hsig(proxI(iP),3)=1E-12.*H0(proxI(iP),2).*ones(length(proxI(iP)),1); %Fix slope
    H0(proxI(iP),2)=obsPrior(proxI(iP)).m; %Proxy Scaling uncertainty prior
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

