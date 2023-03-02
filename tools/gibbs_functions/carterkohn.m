function [x0,contributionChain,x2,F,H]=carterkohn(pS,opts,m)
% carterkohn
% This function applies the Carter and Kohn algorithm to run a Kalman filter
% and backward recursion to estimate the state variable of a time series.
%
% INPUTS:
%
% * |pS|: a structure containing various inputs required for running the
% Kalman filter and backward recursion. It contains the following fields:
% * |hload|: a matrix of factor loadings
% * |T|: the length of the time series
% * |L|: the lag length of the state vector
% * |N|: the number of factors
% * |NN|: the number of observers
% * |ns|: the number of state variables
% * |rmat|: a vector of observer error variance
% * |alpha|: a matrix of lag coefficients
% * |Sigma|: a vector of TSI uncertainty priors in time
% * |data|: a matrix of observations
% * |cmpStYr|: the year for which regression coefficients start being
% estimated
% * |V00|: a matrix of prior state vector variances
% * |glf|: a structure of non-linear function handles
% * |pindex|: a vector of indices of parameters subject to non-linear
% transformation
% * |xbar|: a vector of mean values for each variable used in the
% non-linear transformation
% * |t|: a matrix of time-dependent variables (used if opts.tDependence is true)
% * |floor|: a floor value for the state variable. Used in event of imposed
% minimum solar constant value
% * |m|: the number of MCMC draws
% * |burn|: the burn-in period
%
%
% * |opts|: a structure of options to tweak the function. It contains the
% following fields:
% * |restrictOsf|: a flag indicating whether to restrict the use of open solar
% flux (OSF) in determining TSI during the satellite era
% * |nonLin|: a flag indicating whether the model has non-linear relationships
% * |tDependence|: a flag indicating whether to incorporate time-dependent
% variables
%
% OUTPUTS:
%
% * |x0|: a matrix of estimated state variable values
% * |contributionChain|: a matrix of innovation contribution values
% * |x2|: a matrix of MCMC draws of state variable values
% * |F|: a matrix of transition models used to estimate prior x_i
% * |H|: a matrix of factor loadings

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare variables for Kalman filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract variables from input struct pS
hload=pS.hload; % matrix of factor loadings
T=pS.T; % number of time steps
L=pS.L; % number of lags in the state vector
N=pS.N; % number of variables in the state vector
NN=pS.NN; % number of variables in the observation vector
ns=pS.ns; % number of state variables
rmat=pS.rmat; % matrix of measurement uncertainty
alpha=pS.alpha; % matrix of state vector  autoregression coefficients
Sigma=pS.Sigma; % matrix of time-specific variances for the state vector
data=pS.data; % matrix of observations
X0=pS.X0; %initial expectation for the state vector
V00=pS.V00; % initial covariance matrix for the state vector

% Extract variables from input struct opts
burn=opts.burn;
if isfield(opts,'nonLin') && opts.nonLin % flag for non-linear relationships
glf=pS.glf; % struct containing non-linear functions
pindex=pS.pindex; % vector indicating which variables have non-linear relationships
xbar=pS.xbar; % mean of the state vector
end
if isfield(opts,'tDependence') && opts.tDependence % flag for time-dependent factor loadings
tau=pS.tau; % vector of time steps
end
if isfield(opts,'restrictOsf') && opts.restrictOsf
    oM=~isnan(data);
    oM(pS.dateM.Year>=pS.cmpStYr,pS.osfIndex)=false;
else
    oM=~isnan(data);
end

% Create matrix of factor loadings H
if opts.tDependence % if time-dependent factor loadings
H=zeros(NN,NN+L+1);
H(:,1:2)=hload(:,1:2); % first two columns are the same for all time steps
for ii = 1:NN % set the third column for each variable based on the time step
H(ii,ii+L+1)=hload(ii,3);
end
else % if not time-dependent factor loadings
H=zeros(NN,L+1);
H(:,1:2)=hload(:,1:2);
end

% Create matrix R of measurement uncertainty and vector MU of mean state vector
R=diag(rmat);
MU=[alpha(end,:)';zeros(N*(L-1),1)]';

% Create matrices F and Q for transition model and TSI uncertainty respectively
F=[alpha(1:N*L,:)';eye(N*(L-1),N*L)];
Q=zeros(size(F,1),size(F,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run Kalman filter to estimate x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize variables
xp=[]; % filtered state variable x
vtt=zeros(T,ns,ns); % variance of the state variable x
contributionChain = zeros(T,NN); % contribution of each observation to innovation

t=1;
Q(1:N,1:N)=Sigma(t);

% Prediction for the first step
if opts.tDependence % if time-dependent factor loadings
z10=[1 MU+X0*F' tau(t,:)];
else % if not time-dependent factor loadings
z10=[1 MU+X0*F'];
end

h=H(oM(t,:),:); % subset observation model for observers with observations at i

%subset of coefficients scaling state vectors [x_t x_t-1]
ht=h(:,2:(L+1)); 

%Calculate the predicted covariance matrix of the state variables, v10, 
%based on the previous covariance matrix V00, the state transition matrix F, 
%and the noise covariance matrix Q.
v10=F*V00*F'+Q;
%Calculate the predicted observation vector
yhat=(h*(z10)')';
%Perform update for nonlinear relationships using EKF method, perform a 
%nonlinear regression if the relationship is nonlinear
if opts.nonLin
    idx=find(oM(t,:));
    for iP=1:NN
        if ismember(iP,idx) && pindex(iP) && m>burn
            yhat(idx==iP)=glf(iP).glm(glf(iP).ff,z10(2));
            ht(idx==iP,1)=teglf(z10(2),glf(iP).ff,xbar);
        end
    end
end
xi=data(t,oM(t,:))-yhat;
fxi=(ht*v10*ht')+diag(diag(R(oM(t,:),oM(t,:))));
%updating
K=(v10*ht')*inv(fxi);
contributionChain(t,oM(t,:)) = K(1,:).*xi; %Record contribution of each obs to innovation

%Update the state vector, z11, by adding the Kalman gain times the innovation, 
%xi, to the predicted state vector, z10. If there is time dependence in the 
%model, the time vector, t, is also included.
if opts.tDependence
    z11=[1 (z10(2:(L+1))'+K*xi')' tau(1,:)];
else
    z11=[1 (z10(2:(L+1))'+K*xi')'];
end
v11=v10-K*(ht*v10);
xp=[xp;z11];
vtt(t,:,:)=v11;
%Prediction for other steps
for t=2:T
    Q(1:N,1:N)=Sigma(t);
    h=H(oM(t,:),:); %subset observation model for observers with observations at i
    ht=h(:,2:L+1); %subset of coefficients scaling state vectors [x_t x_t-1]
    %NOTE: AS OF 9/6/22 removed reference to MU in following line, as it's assumed to be
    %0
    if opts.tDependence
        z10=[1 z11(2:(L+1))*F' tau(t,:)];
    else
        z10=[1 z11(2:(L+1))*F'];
    end

    v10=F*v11*F'+Q;
    yhat=(h*(z10)')';
    %Perform update for nonlinear relationships using EKF method
    if opts.nonLin
        idx=find(oM(t,:));
        for iP=1:NN
            if ismember(iP,idx) && pindex(iP) && m>burn
                yhat(idx==iP)=glf(iP).glm(glf(iP).ff,z10(2));
                ht(idx==iP,1)=teglf(z10(2),glf(iP).ff,xbar);
            end
        end
    end
    xi=data(t,oM(t,:))-yhat; 
    fxi=(ht*v10*ht')+diag(diag(R(oM(t,:),oM(t,:))));
    %updating
    K=(v10*ht')*inv(fxi); %Calculate Kalman gain
    contributionChain(t,oM(t,:)) = K(1,:).*xi;
    if opts.tDependence
        z11=[1 (z10(2:(L+1))'+K*xi')' tau(t,:)];
    else
        z11=[1 (z10(2:(L+1))'+K*xi')'];
    end
    v11=v10-K*(ht*v10);
    vtt(t,:,:)=v11;
    xp=[xp;z11];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run backward recursion to determing x_i using Carter Kohn algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize x2, a matrix that will hold the draw of the state variable
% jv and jv1 are used to index the state variables
x2 = zeros(T, ns);
jv = 2:(L+1);
jv1 = 1;
wa = randn(T, ns);
f = F(jv1,:);
mu = 0;
% Start with the last time step
t = T;
p00 = squeeze(vtt(t,jv1,jv1));
% Draw the updated value of x at time i
x2(t, jv1) = xp(t, jv(jv1)) + (wa(t, jv1) * chol(p00));
% Iterate backwards from time T-1 to time 1
for t = T-1:-1:1
% Update the covariance matrices with values at time i
Q(1:N, 1:N) = Sigma(t);
q = Q(jv1, jv1);
% Update pt with value at time i
pt = squeeze(vtt(t,:,:));
% Calculate the updated expectation and variance for time i 
bm=xp(t,jv)+(pt*f'*inv(f*pt*f'+q)*(x2(t+1,jv1)-mu-xp(t,jv)*f')')';
pm=pt-pt*f'*inv(f*pt*f'+q)*f*pt;
% Draw the updated value of x at time i from posterior distribution
x2(t,jv1) = bm(jv1) + (wa(t, jv1) * chol(pm(jv1,jv1)));
end
x0 = x2(:,jv1);
% If the pS struct contains a 'floor' field, set any values in x2 and x0 below
% this floor to the value of the floor
if isfield(pS, 'floor')
x2(x2(:,jv1) < pS.floor, jv1) = pS.floor;
x0(x0 < pS.floor) = pS.floor;
end
end