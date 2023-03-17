function [rmat, mSigma, bSigma] = epsilonmagdependent(pS)
% This function estimates the innovation noise in a time series of solar
% irradiance observations (TSI) as a function of the magnitude of the
% solar cycle anomaly. It returns a vector of innovation noise values
% for each time step (rmat), as well as the coefficients of the linear
% relationship between the magnitude of the anomaly and the innovation
% noise (mSigma and bSigma).
%
% INPUTS:
% pS        : Structure containing the following fields:
% x0        : Nx1 vector of modeled TSI values
% X         : NxP matrix of lagged TSI, where P is lag
% alpha     : Px1 vector of regression coefficients
% eta0      : (optional) Prior value of eta
% Xsig      : (optional) Prior covariance of X
%
% OUTPUTS:
% rmat      : Nx1 vector of innovation noise values
% mSigma    : slope of the relationship between anomaly magnitude and noise
% bSigma    : intercept of the relationship between anomaly magnitude and noise

dateM=pS.dateM(2:end);
x0=pS.x0;
X=pS.X;
alpha=pS.alpha;
eta0=pS.eta0;
Xsig=pS.Xsig;
% Compute the solar cycle anomaly for each observation
Y=x0(2:end);
YCycle = x0 - movmean(x0, 12*11, 'omitnan');
YCycle = YCycle(2:end);
YCycleAll = x0 - movmean(x0, 12*11, 'omitnan');

if isfield(pS,'cmpStYr')
    dateI=dateM.Year>=pS.cmpStYr;
    Y=Y(dateI);
    X=X(dateI,:);
    YCycle=YCycle(dateI);
end
cycMin=min(YCycle);
YCycle = YCycle - cycMin;

% Construct the predictor matrix for the linear regression
pred = [ones(length(YCycle), 1) YCycle];

% Sample VAR covariance for time-dependent X uncertainty
errorsx = Y - X * alpha;
sig=IG(0,0,errorsx);
if isfield(pS,'magDependent') && pS.magDependent
    precX = 1 ./ sig;
    
    % If prior X information is provided, estimate noise as a function of TSI magnitude
    if exist('eta0', 'var') && exist('Xsig', 'var')
        sig0 = diag(Xsig.^2);
        Mx = inv(inv(sig0) + precX .* pred' * pred) * (inv(sig0) * eta0 + precX * pred' * errorsx.^2);
        Vx = inv(inv(sig0) + precX .* pred' * pred);
    else
        % Otherwise, estimate noise without prior X information
        Mx = inv(precX .* pred' * pred) * (precX * pred' * errorsx.^2);
        Vx = inv(precX .* pred' * pred);
    end
    
    % Estimate the magnitude-dependent innovation noise using a Gaussian process
    p = -1;
    ind=0;
    while p <= 0
        b = Mx + (randn(1, 2) * chol(Vx))';
        if ind>=10
            b(1)=sig;
        end
        p = b(1);
        ind=ind+1;
    end
    if b(2)<0
        b(2)=0;
    end
    
    % Compute the innovation noise for each observation as a function of the anomaly magnitude
    predAll = [ones(length(YCycleAll), 1) YCycleAll-cycMin];
    rmat = predAll * b;
    rmat(rmat < 0.1.*sig) = 0.1.*sig; %To avoid stability issues near zero
    
    % Output the slope and intercept of the linear relationship between anomaly magnitude and noise
    bSigma = b(1);
    mSigma = b(2);
else
    T0=ceil(length(x0).*(.25./(1-.25))); %Give prior estimate weight
    b = IG(T0, (T0-1).*eta0(1), errorsx); % Estimate of baseline TSI innovation noise
    rmat=b.*ones(length(YCycleAll),1);
    bSigma=b;
    mSigma=0;
end

end

