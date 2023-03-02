function [rmat, mSigma, bSigma] = epsilonmagdependent(x0, X, alpha, X0, Xsig)
% This function estimates the innovation noise in a time series of solar
% irradiance observations (TSI) as a function of the magnitude of the
% solar cycle anomaly. It returns a vector of innovation noise values
% for each time step (rmat), as well as the coefficients of the linear
% relationship between the magnitude of the anomaly and the innovation
% noise (mSigma and bSigma).
%
% INPUTS:
% x0        : Nx1 vector of modeled TSI values
% X         : NxP matrix of lagged TSI, where P is lag
% alpha     : Px1 vector of regression coefficients
% X0        : (optional) Prior value of X
% Xsig      : (optional) Prior covariance of X
%
% OUTPUTS:
% rmat      : Nx1 vector of innovation noise values
% mSigma    : slope of the relationship between anomaly magnitude and noise
% bSigma    : intercept of the relationship between anomaly magnitude and noise

% Compute the solar cycle anomaly for each observation
Y=x0(2:end);
YCycleAll = x0 - movmean(x0, 12*11, 'omitnan');
YCycleAll = YCycleAll - min(YCycleAll);
YCycle = YCycleAll(2:end);

% Construct the predictor matrix for the linear regression
pred = [ones(length(YCycle), 1) YCycle];

% Sample VAR covariance for time-dependent X uncertainty
errorsx = Y - X * alpha;
sig = IG(0, 0, errorsx); % Estimate of baseline TSI innovation noise
precX = 1 ./ sig;

% If prior X information is provided, estimate noise as a function of TSI magnitude
if exist('X0', 'var') && exist('Xsig', 'var')
    sig0 = diag(Xsig.^2);
    Mx = inv(inv(sig0) + precX .* pred' * pred) * (inv(sig0) * X0 + precX * pred' * abs(errorsx));
    Vx = inv(inv(sig0) + precX .* pred' * pred);
else
    % Otherwise, estimate noise without prior X information
    Mx = inv(precX .* pred' * pred) * (precX * pred' * abs(errorsx));
    Vx = inv(precX .* pred' * pred);
end

% Estimate the magnitude-dependent innovation noise using a Gaussian process
p = -1;
while p <= 0
    b = Mx + (randn(1, 2) * chol(Vx))';
    p = b(1);
end

% Compute the innovation noise for each observation as a function of the anomaly magnitude
rmat = (pred * b).^2;
rmat = [rmat(1); rmat];
rmat(rmat < sig) = sig;

% Output the slope and intercept of the linear relationship between anomaly magnitude and noise
bSigma = b(1);
mSigma = b(2);

end

