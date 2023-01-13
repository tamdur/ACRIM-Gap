%Determine a most-likely estimate for the prior uncertainty of observers
%and the underlying TSI process, including all satellites and proxies
%
% Ted Amdur
% 8/23/21, updated 11/10/21 to homogenize process of determining priors for
% satellites and proxies. Updated 3/11/22 to use revised scripts and
% updated data

clearvars

%First, get the array of all observations to be compared
load obs_22_6_18.mat %From makeobsmatrix.m
nObs = length(colLabels);
%Make labels for which records fall under satellite vs proxy
proxInd=[4 8];
isSat = true(nObs,1);
isSat(proxInd) = false; 
isProx = false(nObs,1);
isProx(proxInd) = true;

%Make a set of priors for the proxy observations using all satellites
% labProx = colLabels(isProx);
% valProx = valM(:,isProx);
for ii = 1:2
    ind=proxInd(ii);
    pRec = nansum(valM(:,((ii-1)*4+1):(ii*4)),2);
    oRec = nansum(oM(:,((ii-1)*4+1):(ii*4)),2);
    vProx = NaN(nObs,1);
    iS = 1;
    vProx = [];
    bProx = [];
    bIntProx = [];
    while iS <= nObs
        if ~any(iS==proxInd) %only operate on satellites
            overlap = and(oM(:,iS),oM(:,ind));
            if sum(overlap) > 60 %5 year cutoff for comparison to be made
                %Revised 9/8/21 to be in native units
                [b,bint] = regress(valM(overlap,ind),...
                    [ones(sum(overlap),1) valM(overlap,iS)-nanmean(valM(overlap,iS))]);
                res = valM(overlap,ind) - b(1).*ones(sum(overlap),1) - b(2).*valM(overlap,iS);
                vProx = [vProx; nanvar(res,1)];
                bProx = [bProx; b'];
                bIntProx = [bIntProx; bint(2,:)];
            end
        end
        iS = iS + 1;
    end
    obsPrior(proxInd(ii)).type = "proxy";
    obsPrior(proxInd(ii)).name = colLabels(ind);
    obsPrior(proxInd(ii)).std = sqrt(nanmean(vProx));
    obsPrior(proxInd(ii)).b = mean(bProx(:,1));
    obsPrior(proxInd(ii)).m = mean(bProx(:,2));
    obsPrior(proxInd(ii)).bsig = std(bProx(:,1)); %One sigma uncertainty in offset
    obsPrior(proxInd(ii)).msig = std(bProx(:,2)); %One sigma uncertainty in scaling
end

%Then, make a set of comparisons between all the satellites
for ii=1:nObs
    if ~any(ii==proxInd) %Only make comparisons for satellites
        vSat = NaN(nObs,1);
        iS = 1;
        while iS <= nObs
            if iS ~= ii && isSat(iS)
                overlap = logical(oM(:,iS).*oM(:,ii));
                s1 = valM(overlap,ii) - mean(valM(overlap,ii));
                s2 = valM(overlap,iS) - mean(valM(overlap,iS));
                dSat = s2-s1;
                vSat(iS) = nanvar(dSat,1);
            end
            iS = iS + 1;
        end
        %Include the data collected into a structure object
        obsPrior(ii).type = "sat";
        obsPrior(ii).name = colLabels(ii);
        obsPrior(ii).std = sqrt(nanmean(vSat))./sqrt(2); %Correct for satellite noise coming from two observers
    end
end

save('mat_files/obspriors_22_06_23.mat','obsPrior')

