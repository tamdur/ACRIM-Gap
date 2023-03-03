function[hload,errorN,infI]=coeffsample(pS,opts)
% coeffsample - Sample loadings that compose H through Bayesian regression
%
% Syntax: [hload,errorN,infI]=coeffsample(pS,opts)
%
% Inputs:
%   pS - struct, contains the following fields:
%       .data - array of observations
%       .oM - logical array of same size as data, indicates which entries
%       have observations
%       .rmat - array of observer noise estimates
%       .dateM - datetime array corresponding to the observations
%       .cx - array of ones to be used for the offset coefficient in H
%       .tau - array of time rows for t (only used if pS.tDependence is true)
%       .x0 - array of initial guesses for the process
%       .tDependence - logical, true if drift predictor is used, false otherwise
%       .H0 - prior mean for observation model
%       .Hsig - prior variance for observation model
%       .NN - number of observers
%   opts - struct, contains the following fields:
%       .cmpStYr - the starting year of the comparison period used for drift predictor
%       .cmpStYrRegress - logical, true if drift predictor is used, false otherwise
%
% Outputs:
%   hload - array of the factor dependent coefficients for H
%   errorN - array of observation errors
%   infI - logical array indicating which observation entries were used

data=pS.data;
oM=pS.oM;
rmat=pS.rmat;
dateM=pS.dateM;
cx=pS.cx;
if pS.tDependence
    tau=pS.tau;
end
H0=pS.H0;
Hsig=pS.Hsig;
x0=pS.x0;

if isfield(opts,'cmpStYrRegress') && opts.cmpStYrRegress % Use drift predictor
    infI=logical((dateM.Year>=opts.cmpStYr).*oM);
else
    infI=oM;
end

hload=[];
errorN=NaN(size(data));

for i=1:pS.NN % Loop over all observers
    y=data(infI(:,i),i);
    
    % Z needs to be the full predictor matrix
    if pS.tDependence
        Z = [cx(infI(:,i)) x0(infI(:,i)) tau(infI(:,i),i)];
    else
        Z = [cx(infI(:,i)) x0(infI(:,i))];
    end
    
    precY=1./rmat(i); % Precision of observer i
    sH=diag(Hsig(i,:));
    
    M=inv(inv(sH) + precY.*Z'*Z)*(inv(sH)*H0(i,:)'+precY*Z'*y); % Posterior mean
    V=inv(inv(sH) + precY.*Z'*Z); % Posterior variance
    
    % Draw from multivariate normal distribution with mean M and variance V
    hf=M+(randn(1,size(Z,2))*chol(V))';
    
    % Calculate observation errors
    errorN(infI(:,i),i)=y-Z*hf;
    
    hload=[hload;hf']; % Factor dependent coefficients for H
end
end