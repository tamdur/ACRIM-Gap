%runchain_22_04_25.m PRODUCE GIBBS SAMPLING OUTPUT OF BAYESIAN KALMAN
%FILTER MODEL USING SATLELLITE AND PROXY DATA.
%
% Ted Amdur 04/25/22
%Based upon runchain_22_03_23b.m that shows innovations with time, also
%allows for different ar order models. The following is an
%extensively-modified version of the code from Ex. 4 of Chapter 4 of
%Applied Bayesian Econometrics for Central Bankers, by Andrew Blake and
%Haroon Mumtaz.

clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Headers to modify
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveString = 'chain_output/ar2_23_01_14.mat';
outDat.script=mfilename; %Save name of script
obsmatrix='obs_23_01_13'; %Load data array, with colLabels corresponding to observer source for each column
excludeFliers=0;%1 to remove outlier observations from examined dataset

%Develop a set of monthly observations from satellite and proxy
%observations
dateS=getdates;
dateCycles=dateS.cycles;


load(obsmatrix); %From makeobsmatrix.m
%Twelve columns of valM correspond to the following observers:
%     "ACRIM1/SMM"
%     "ACRIM2/UARS"
%     "BremenMgII"
%     "ERBE/ERBS"
%     "NIMBUS-7/HF"
%     "SILSO"
%     "VIRGO/SOHO"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Beginning of main script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if excludeFliers %Code to remove outliers using a past run of BTSI
    load excludeMask_22_11_03.mat %from exclude_fliers_22_04_26.m
    valM(excludeMask) = NaN;
    oM(excludeMask) = false;
end
%Specify priors for H coefficients
oindex=[1 1 1 1 1 1 0]; %oindex=1 for observers with varying offset, 0 for fixed
tindex=[1 1 0 1 1 0 1]; %tindex=1 for observers with time dependent drift, 0 otherwise
sindex=[0 0 1 0 0 1 0]; %sindex=1 for observers with non-identity scaling to TSI, 0 otherwise
satindex=tindex; %satindex=1 for observers that are satellites, 0 otherwise

valMAll=valM;
data = valM;
iObs = nansum(oM,1);
T=size(data,1);
KK=1;  %number of factors
L=2;  %number of lags in the VAR
N=KK; %number of Variables in transition equation
NN=size(data,2);% Number of observers

cx = ones(T,1); %ones vector in first row of z for offset
t =repmat(linspace(0,T./120,T)',[1 size(data,2)]); %Make time rows for t
for ii=1:size(data,2)
    TM=mean(t(oM(:,ii),ii));
    if isnan(TM)
        TM=0;
    end
    t(:,ii)=t(:,ii)-TM;
end

%Load priors for observation model
opts.normalize=0; %Set to 1 to normalize data within
[H0, Hsig, T0, th0] = getpriors(NN,oindex,tindex,sindex,valM,oM,...
    satindex,colLabels,opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 1: establish starting values and priors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get an intial guess for the process
x0=zeros(T,1);
guess1 = nanmean(data(:,~sindex),2); %mean of satellite observations
x0(~isnan(guess1))=guess1(~isnan(guess1));
X0=[x0(1) zeros(1,L-1)];%state vector X[t-1|t-1] of x and nth order lagged x
ns=size(X0,2);
V00=eye(ns);  %v[t-1|t-1]
rmat=1E9.*ones(NN,1); %arbitrary starting value for the variance of process x
Sigma=eye(N);  %arbitrary starting value for the variance of transition model errors
%Save the records of contributions to innovation at each time i
contributionChain = NaN(T,size(data,2));

reps=10500; %Total steps of Gibbs Sampler
burn=500; %Steps excluded in burnin period
mm=1;%index for saved chain post-burnin
tic
for m=1:reps %Gibbs sampling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 2: sample loadings that compose H through Bayesian regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hload=[];
error=NaN(size(data));
for i=1:NN %loop over all observers
    infI(:,i)=oM(:,i);
    y=data(infI(:,i),i);
    alpha = [cx(infI(:,i)) x0(infI(:,i)) t(infI(:,i),i)]; %This needs to be the full predictor matrix
    precY=1./rmat(i); %Precision of observer i
    sH=diag(Hsig(i,:));
    M=inv(inv(sH) + precY.*alpha'*alpha)*(inv(sH)*H0(i,:)'+precY*alpha'*y);
    V=inv(inv(sH) + precY.*alpha'*alpha);
    %draw
    hf=M+(randn(1,size(alpha,2))*chol(V))';
    hload=[hload;hf']; %The factor dependent coefficients for H
    error(infI(:,i),i)=y-alpha*hf;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 3: sample variance of the observers from inverse gamma distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmat=[];
for i=1:NN
    rmati= IG(T0(i),th0(i),error(infI(:,i),i));
    th(i)=error(infI(:,i),i)'*error(infI(:,i),i);
    rmat=[rmat rmati];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 4: sample estimates of autoregressive parameters alpha
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=x0;
X=[];
for iL=1:L %Create columns for lagged estimates of x
    X=[X lag0(Y,iL)];
end
X=[X ones(size(Y,1),1)];%To ensure a stable linear regression, estimate mean of TSI
Y=Y(2:end,:);
X=X(2:end,:);

M=inv(X'*X)*(X'*Y);M=M(:);  %conditional mean (right now just the obs mean) 
V=mean(Sigma).*inv(X'*X); %conditional variance
chck=-1;                 %make sure VAR is stationary
while chck<0
alpha=M+(randn(1,N*(N*L+1))*chol(V))';  %draw for VAR coefficients
S=stability(alpha,N,L);
if S==0
    chck=10;
end
end
alpha1=reshape(alpha,N*L+1,N);

errorsx=Y-X*alpha1;

%sample VAR covariance for time-dependent X uncertainty
bSigma=IG(0,0,errorsx); %Estimate of baseline TSI innovation noise
precX=1./bSigma;
Mx=inv(precX.*Y'*Y)*(precX*Y'*(errorsx.^2)); %Estimate noise as fn of TSI magnitude
Vx=bSigma.*inv(Y'*Y);
mSigma=Mx+(randn*chol(Vx))'; %Estimate of magnitude-dependent innovation noise
Sigma=bSigma+mSigma*Y;
Sigma(Sigma<bSigma)=bSigma;
Sigma=[Sigma(1);Sigma];%lengthen to full observation interval

%Create matrix of factor loadings
H=zeros(NN,NN+L+1);
H(:,1:2)=hload(:,1:2);
for ii = 1:NN
    H(ii,ii+L+1)=hload(ii,3);
end
%matrix R of measurement uncertainty
R=diag(rmat);
%vector MU of state vector mean state
MU=[alpha1(end,:)';zeros(N*(L-1),1)]';
%matrix F, transition model used to estimate prior \hat{x}_i
F=[alpha1(1:N*L,:)';eye(N*(L-1),N*L)];
%matrix Q of TSI uncertainty \eta at time i
Q=zeros(size(F,1),size(F,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 5: Run Kalman filter to estimate x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xp=[];          %will hold the filtered state variable x'
vtt=zeros(T,ns,ns);    % will hold its variance
i=1;
Q(1:N,1:N)=Sigma(i);
h=H(oM(i,:),:);%x is now a matrix of coefficients
ht=h(:,2:(L+1)); %subset of coefficients scaling state vectors [x_t x_t-1]
%Prediction for the first step
z10=[1 MU+X0*F' t(i,:)];
v10=F*V00*F'+Q;
yhat=(h*(z10)')';
xi=data(i,oM(i,:))-yhat;
fxi=(ht*v10*ht')+diag(diag(R(oM(i,:),oM(i,:))));
%updating
K=(v10*ht')*inv(fxi);
contributionChain(i,oM(i,:)) = K(1,:).*xi; %Record contribution of each obs to innovation
z11=[1 (z10(2:(L+1))'+K*xi')' t(1,:)];
v11=v10-K*(ht*v10);
xp=[xp;z11];
vtt(i,:,:)=v11;
%Prediction for other steps
for i=2:T
    Q(1:N,1:N)=Sigma(i);
    h=H(oM(i,:),:); %subset observation model for observers with observations at i
    ht=h(:,2:L+1); %Part of observation model responsible for scaling x
    z10=[1 MU+z11(2:(L+1))*F' t(i,:)];
    v10=F*v11*F'+Q;
    yhat=(h*(z10)')';
    xi=data(i,oM(i,:))-yhat; 
    fxi=(ht*v10*ht')+diag(diag(R(oM(i,:),oM(i,:))));
    %updating
    K=(v10*ht')*inv(fxi); %Calculate Kalman gain
    contributionChain(i,oM(i,:)) = K(1,:).*xi;
    z11=[1 (z10(2:(L+1))'+K*xi')' t(i,:)];
    v11=v10-K*(ht*v10);
    vtt(i,:,:)=v11;
    xp=[xp;z11];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 6: Run backward recursion to determing x_i using Carter Kohn
%algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x2 = zeros(T,ns);   %this will hold the draw of the state variable
jv=2:(L+1); jv1=1; %index of state variables to extract
wa=randn(T,ns);
f=F(jv1,:);
q=Q(jv1,jv1);
mu=MU(jv1);
i=T;  %period t
p00=squeeze(vtt(i,jv1,jv1)); 
%draw for updated x in time i
x2(i,jv1)=xp(i,jv(jv1))+(wa(i,jv1)*chol(p00));   
%periods t-1..to .1
for i=T-1:-1:1
    Q(1:N,1:N)=Sigma(i);q=Q(jv1,jv1);
    pt=squeeze(vtt(i,:,:));
    bm=xp(i,jv)+(pt*f'*inv(f*pt*f'+q)*(x2(i+1,jv1)-mu-xp(i,jv)*f')')';
    pm=pt-pt*f'*inv(f*pt*f'+q)*f*pt;
    x2(i,jv1)=bm(jv1)+(wa(i,jv1)*chol(pm(jv1,jv1)));
end
x0=x2(:,jv1);   %update the state variable estimate
if mod(m,100) == 0 %display progress of chain
    disp([num2str(m) ' reps completed'])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%step 7: If burnin completed, store state and observation model estimates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if m>burn
    %save
    xAll(:,mm)=x2(:,1); %Reconstructed TSI
    for ii=1:NN
        A(ii,:,mm)=H(ii,[1 2 ii+L+1])'; %observation coefficients
    end
    a(:,mm)=F(1,:); %VAR autoregressive coefficients [lag1 lag2]
    sigY(:,mm)=rmat; %observer noise estimate
    sigX(:,mm) = Sigma; %TSI noise estimate
    theta(:,mm)=th; %Theta parameter estimate
    outDat.contributionChain(mm,:,:)=contributionChain; %Innovation contributions
    mm=mm+1;
end
    

end

toc
outDat.reps=reps;outDat.burn=burn;outDat.H0=H0;outDat.Hsig=Hsig;outDat.T0=T0;
outDat.th0=th0;outDat.oindex=oindex;outDat.tindex=tindex;outDat.sindex=sindex;
outDat.satindex=satindex;
outDat.excludeFliers=excludeFliers;
outDat.obsmatrix=obsmatrix;
save(saveString,'xAll','sigY','sigX','theta','a','A','t','outDat','-v7.3')
out1x=prctile(xAll',[10 20 30 40 50 60 70 80 90])';



