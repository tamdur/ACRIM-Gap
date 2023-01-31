function[H0, Hsig, T0, th0] = getpriors(NN,oindex,tindex,sindex,valM,oM,...
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
Hsig(:,1)=Hsig(:,1)+(offSig.*oindex');%prior sigma for offset
Hsig(:,2)=Hsig(:,2)+(1000.*sindex');%prior sigma for scaling
Hsig(:,3)=Hsig(:,3)+(mSig.*tindex'); %prior sigma for drift

%Make changes for proxies if linearity assumed (first row sunspots, second row mg)
proxI=find(~satindex);%indices for proxies
for iP=1:length(proxI)
    Hsig(proxI(iP),1)=obsPrior(proxI(iP)).bsig.^2; %Var uncertainty for proxy offset
    Hsig(proxI(iP),2)=obsPrior(proxI(iP)).msig.^2; %var uncertainty for proxy scaling
    Hsig(proxI(iP),3)=1E-12.*ones(length(proxI(iP)),1); %Fix slope
    H0(proxI(iP),2)=obsPrior(proxI(iP)).m; %Proxy Scaling uncertainty prior
end
%Next, draw a range of initial hyperparameters for the prior noise variance
iObs=sum(oM,1);
T0=ceil(iObs.*(priFrac./(1-priFrac))); %Assumes std from intercomparison
for ii = 1:NN %Prior for proxy noise
    th0(ii) = (T0(ii)-1).*obsPrior(ii).std.^2;
end

end

