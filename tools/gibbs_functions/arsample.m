function [alpha,X,Y]=arsample(pS,cmpStYr)
% ARSAMPLE draws samples of autoregressive parameters from a conditional
% multivariate normal distribution given observations
%
%   [alpha,X,Y] = arsample(pS, cmpStYr) returns the VAR coefficients
%   sampled from a conditional multivariate normal distribution given
%   observations. 
%
%   Inputs:
%       - pS: structure containing the following fields:
%           - dateM: date vector of the input data
%           - x0: initial state vector
%           - N: number of variables
%           - L: number of lags in the VAR
%           - Sigma: covariance of the observation error
%           - T: number of observations
%       - cmpStYr: starting year for the regression
%
%   Outputs:
%       - alpha: VAR coefficients (matrix of size N*(L+1) by N)
%       - X: regressor matrix (matrix of size T-L by N*(L+1))
%       - Y: response variable matrix (matrix of size T-L by N)
%
% Ted Amdur, 2023/03/02
% Adapted from code by Andrew Blake and Haroon Mumtaz.


% Select observations from input data from starting year cmpStYr onwards
infI = pS.dateM.Year >= cmpStYr;

% Obtain the response variable matrix
Y = pS.x0(infI);

% Extract the number of variables and the number of lags in the VAR
N = pS.N;
L = pS.L;

% Create regressor matrix by lagging the response variable matrix
X = [];
for iL = 1:L
    X = [X lag0(Y, iL)]; %lag0 is a function that lags the data by iL periods
end

% Add a column of ones to X to ensure a stable linear regression, estimate
% mean of TSI
X = [X ones(size(Y, 1), 1)];

% Remove the first observation from Y and X (since there are no lags)
Y = Y(2:end, :);
X = X(2:end, :);

% Calculate the conditional mean and variance of the VAR coefficients
M = inv(X' * X) * (X' * Y);
M = M(:);  % columnize M
V = mean(pS.Sigma) .* inv(X' * X);

% Draw VAR coefficients from the conditional multivariate normal
% distribution until a stationary VAR is obtained
chck = -1;
while chck < 0
    alphaT = M + (randn(1, N * (N * L + 1)) * chol(V))';
    S = stability(alphaT, N, L);
    if S == 0
        chck = 10;
    end
end

% Reshape the vector of VAR coefficients into a matrix
alpha = reshape(alphaT, N * L + 1, N);
% Note that alpha is a matrix of size N*(L+1) by N, with each column
% corresponding to the VAR coefficients for a given variable

% The regressor matrix X is returned for potential use in subsequent
% calculations
end
