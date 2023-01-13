function[H0, Hsig, T0, th0] = getar2priorseval(NN,oindex,tindex,sindex,iObs,colLabels,obsPrior)
priFrac = 0.25; %Fraction of the 'observations' for sigX, sigY to be from prior
offSig = 5; %prior variance for a (offset)
mSig = 0.25; %prior variance for c (drift)
H0=zeros(NN,3);H0(:,2) = 1; %prior for scaling H of form [a_i b_i c_i]
Hsig=1E-12.*ones(NN,3);
Hsig(:,1)=Hsig(:,1)+(offSig.*oindex');%prior sigma for offset
Hsig(:,2)=Hsig(:,2)+(1000.*sindex');%prior sigma for scaling
Hsig(:,3)=Hsig(:,3)+(mSig.*tindex'); %prior sigma for drift
%Make changes for proxies (first row sunspots, second row mg)
proxI=find(sindex);%indices for proxies
spotI=find(strcmp(colLabels,"SILSO"));
mgI=find(strcmp(colLabels,"BremenMgII"));
Hsig(spotI,1)=obsPrior(spotI).bsig.^2; %Var Uncertainty for sunspots
Hsig(mgI,1)=obsPrior(mgI).bsig.^2; %Var Uncertainty for mg
H0(spotI,2)=obsPrior(spotI).m; %sunspot scaling
H0(mgI,2)=obsPrior(mgI).m; %mgII scaling
Hsig(spotI,2)=obsPrior(spotI).msig.^2; %var for spot scaling
Hsig(mgI,2)=obsPrior(mgI).msig.^2; %var for mgII scaling
Hsig(proxI,3)=1E-12.*ones(2,1); %Fix slope
%Next, draw a range of initial hyperparameters for the prior noise variance
T0=ceil(iObs.*(priFrac./(1-priFrac))); %Assumes std from intercomparison
for ii = 1:NN %Prior for proxy noise
    th0(ii) = (T0(ii)-1).*obsPrior(ii).std.^2;
end

end

