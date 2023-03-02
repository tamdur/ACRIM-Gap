function [xAll,sigY,sigX,theta,a,A,tau,outDat] = runchain_23_02_21(valM,oM,colLabels,opts)
%runchain_23_02_21.m PRODUCE GIBBS SAMPLING OUTPUT OF BAYESIAN KALMAN
%FILTER MODEL USING SATLELLITE AND PROXY DATA.
%
% Ted Amdur 2023/02/21
%Based upon runchain_22_03_23b.m that shows innovations with time, also
%allows for different ar order models. The following is an
%extensively-modified version of the code from Ex. 4 of Chapter 3 of
%Applied Bayesian Econometrics for Central Bankers, by Andrew Blake and
%Haroon Mumtaz.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Headers to modify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outDat.script=mfilename; %Save name of script

if ~exist('valM','var') || isempty(valM) %Load default observation array, otherwise load provided one
obsmatrix='obs_23_02_01'; %Load data array, with colLabels corresponding to observer source for each column
load(obsmatrix); %From makeobsmatrix.m
% outDat.mgScale=nanstd(valM(:,1));
% valM(:,1)=valM(:,1)./nanstd(valM(:,1));
%Seven columns of valM correspond to the following observers:
% "BremenMgII"
% "ERBS/ERBE"
% "NIMBUS-7/HF"
% "SILSO"
% "SMM/ACRIM1"
% "SOHO/VIRGO"
% "UARS/ACRIM2"

else
    obsmatrix='synthetic';
end
if ~isfield(opts,'randomizeChain') %Load default observation array, otherwise load provided one
    opts.randomizeChain=false;
end
if ~isfield(opts,'excludeFliers')
    opts.excludeFliers=false;%true to remove outlier observations from examined dataset
end
if ~isfield(opts,'reps')
    opts.reps=10500; %Total steps of Gibbs Sampler
end
if ~isfield(opts,'burn')
    opts.burn=500; %Steps excluded in burnin period
end
if ~isfield(opts,'logContributions')
    opts.logContributions=false; %Set to true to include record of innovation contributions
end
if ~isfield(opts,'dispProgress')
    opts.dispProgress=false; %Set to true to display progress of every 100 iterations
end
if ~isfield(opts,'lags')
    opts.lags=2; %Set process to AR(2)
end
opts.normalize=false; %Set to true to normalize data within. For getpriors call

%Develop a set of monthly observations from satellite and proxy
%observations
dateS=getdates;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of main script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.excludeFliers %Code to remove outliers using a past run of BTSI
    %load excludeMask_22_11_03.mat %from exclude_fliers_22_04_26.m
    valM(excludeMask) = NaN;
    oM(excludeMask) = false;
end
%Specify priors for H coefficients
oindex=[1 0 1 1 1 1 1]; %oindex=1 for observers with varying offset, 0 for fixed
tindex=[1 1 0 1 1 0 1]; %tindex=1 for observers with time dependent drift, 0 otherwise
sindex=[0 0 1 0 0 1 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
satindex=tindex; %satindex=1 for observers that are satellites, 0 otherwise

valMAll=valM;
tsi.data = valM;
iObs = nansum(oM,1);
tsi.T=size(tsi.data,1);
tsi.L=opts.lags;  %number of lags in the VAR
tsi.N=1; %Number of factors/variables in the transition equation
tsi.NN=size(tsi.data,2);% Number of observers

cx = ones(tsi.T,1); %ones vector in first row of z for offset
tau =repmat(linspace(0,tsi.T./120,tsi.T)',[1 tsi.NN]); %Make time rows for t
for ii=1:tsi.NN
    TM=mean(tau(oM(:,ii),ii));
    if isnan(TM)
        TM=0;
    end
    tau(:,ii)=tau(:,ii)-TM;
end
tsi.tau=tau;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 1: establish starting values and priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Load priors for observation model
[H0, Hsig, T0, th0,tsi.Xprior,tsi.Xsig] = getpriors(tsi.NN,oindex,tindex,sindex,valM,oM,...
    satindex,colLabels,opts);

%get an intial guess for the process
x0=zeros(tsi.T,1);
guess1 = nanmean(tsi.data(:,~sindex),2); %mean of satellite observations
x0(~isnan(guess1))=guess1(~isnan(guess1));
if opts.randomizeChain
    rng('shuffle')
    x0=x0+randn(length(x0),1).*(nanstd(x0)./5);
end
tsi.x0=x0;
tsi.X0=[x0(1) zeros(1,tsi.L-1)];%state vector X[t-1|t-1] of x and nth order lagged x
tsi.ns=size(tsi.X0,2);
ns=size(tsi.X0,2);
tsi.V00=eye(tsi.ns);  %v[t-1|t-1]
rmat=1E9.*ones(tsi.NN,1); %arbitrary starting value for the variance of process x
Sigma=eye(tsi.N);  %arbitrary starting value for the variance of transition model errors
%Save the records of contributions to innovation at each time i
contributionChain = NaN(tsi.T,tsi.NN);

%Prepare arrays to be saved
NS=opts.reps-opts.burn;
outDat.contributionChain=zeros(NS,tsi.T,tsi.NN);
outDat.m=zeros(NS,1);
outDat.b=zeros(NS,1);
xAll=zeros(tsi.T,NS);
A=zeros(tsi.NN,3,NS);
a=zeros(tsi.L,NS);
sigY=zeros(tsi.NN,NS);
sigX=zeros(tsi.T,NS);
theta=zeros(tsi.NN,NS);

mm=1;%index for saved chain post-burnin
tic
for m=1:opts.reps %Gibbs sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 2: sample loadings that compose H through Bayesian regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hload=[];
error=NaN(tsi.T,tsi.NN);
for ii=1:tsi.NN %loop over all observers
    infI(:,ii)=oM(:,ii);
    y=tsi.data(infI(:,ii),ii);
    alpha = [cx(infI(:,ii)) tsi.x0(infI(:,ii)) tau(infI(:,ii),ii)]; %This needs to be the full predictor matrix
    precY=1./rmat(ii); %Precision of observer i
    sH=diag(Hsig(ii,:));
    M=inv(inv(sH) + precY.*alpha'*alpha)*(inv(sH)*H0(ii,:)'+precY*alpha'*y);
    V=inv(inv(sH) + precY.*alpha'*alpha);
    %draw
    hf=M+(randn(1,size(alpha,2))*chol(V))';
    hload=[hload;hf']; %The factor dependent coefficients for H
    error(infI(:,ii),ii)=y-alpha*hf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 3: sample variance of the observers from inverse gamma distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmat=[];
for ii=1:tsi.NN
    rmati= IG(T0(ii),th0(ii),error(infI(:,ii),ii));
    th(ii)=error(infI(:,ii),ii)'*error(infI(:,ii),ii);
    rmat=[rmat rmati];
end
tsi.rmat=rmat;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 4: sample estimates of autoregressive parameters alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=tsi.x0;
X=[];
for iL=1:tsi.L %Create columns for lagged estimates of x
    X=[X lag0(Y,iL)];
end
X=[X ones(size(Y,1),1)];%To ensure a stable linear regression, estimate mean of TSI
Y=Y(2:end,:);
X=X(2:end,:);
tsi.X=X;

M=inv(X'*X)*(X'*Y);M=M(:);  %conditional mean (right now just the obs mean) 
V=mean(Sigma).*inv(X'*X); %conditional variance
chck=-1;                 %make sure VAR is stationary
while chck<0
alpha=M+(randn(1,tsi.N*(tsi.N*tsi.L+1))*chol(V))';  %draw for VAR coefficients
S=stability(alpha,tsi.N,tsi.L);
if S==0
    chck=10;
end
end
alpha1=reshape(alpha,tsi.N*tsi.L+1,tsi.N);
tsi.alpha=alpha1;

errorsx=Y-X*alpha1;

%sample VAR covariance for time-dependent X uncertainty
[Sigma,mSigma,bSigma]=epsilonmagdependent(tsi.x0,tsi.X,tsi.alpha,tsi.Xprior,tsi.Xsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 5: Run Kalman filter to estimate x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tsi.Sigma=Sigma;
tsi.hload=hload;
opts.nonLin=false;
opts.tDependence=true;
[tsi.x0,tsi.contributionChain,tsi.x2,tsi.F,tsi.H]=carterkohn(tsi,opts,m);

if opts.dispProgress
    if mod(m,100) == 0 %display progress of chain
        disp([num2str(m) ' reps completed'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 7: If burnin completed, store state and observation model estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if m>opts.burn
    %save
    xAll(:,mm)=tsi.x2(:,1); %Reconstructed TSI
    for ii=1:tsi.NN
        A(ii,:,mm)=tsi.H(ii,[1 2 ii+tsi.L+1])'; %observation coefficients
    end
    a(:,mm)=tsi.F(1,:); %VAR autoregressive coefficients [lag1 lag2]
    sigY(:,mm)=tsi.rmat; %observer noise estimate
    sigX(:,mm) = tsi.Sigma; %TSI noise estimate
    theta(:,mm)=th; %Theta parameter estimate
    if exist('mSigma','var') && exist('bSigma','var')
        outDat.m(mm)=mSigma;
        outDat.b(mm)=bSigma;
    end
    if opts.logContributions
        outDat.contributionChain(mm,:,:)=contributionChain; %Innovation contributions
    end
    mm=mm+1;
end
    

end
if opts.dispProgress
    toc
end
outDat.H0=H0;outDat.Hsig=Hsig;outDat.T0=T0;
outDat.th0=th0;outDat.oindex=oindex;outDat.tindex=tindex;outDat.sindex=sindex;
outDat.satindex=satindex;
outDat.obsmatrix=obsmatrix;
outDat.opts=opts; %Save the input info
outDat.runDate=datetime;
if isfield(opts,'saveFile')
    save(opts.saveFile,'xAll','sigY','sigX','theta','a','A','tau','outDat','-v7.3')
    %out1x=prctile(xAll',[10 20 30 40 50 60 70 80 90])';
end


end
