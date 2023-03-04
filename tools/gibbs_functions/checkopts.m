function optsOut = checkopts(optsIn)
%CHECKOPTS Create default settings for all Gibbs sampling options that have
%not been explicitly specified by the user


if ~isfield(optsIn,'randomizeChain') %Load default observation array, otherwise load provided one
    optsIn.randomizeChain=false;
end
if ~isfield(optsIn,'excludeFliers')
    optsIn.excludeFliers=false;%true to remove outlier observations from examined dataset
end
if ~isfield(optsIn,'reps')
    optsIn.reps=10500; %Total steps of Gibbs Sampler
end
if ~isfield(optsIn,'burn')
    optsIn.burn=500; %Steps excluded in burnin period
end
if ~isfield(optsIn,'logContributions')
    optsIn.logContributions=false; %Set to true to include record of innovation contributions
end
if ~isfield(optsIn,'dispProgress')
    optsIn.dispProgress=false; %Set to true to display progress of every 100 iterations
end
if ~isfield(optsIn,'lags')
    optsIn.lags=2; %Set process to AR(2)
end
optsIn.normalize=false; %Set to true to normalize data within. For getpriors call
optsIn.cmpStYr=1978;
optsOut=optsIn;
end

