%For a given number of chains and a length of run, examine the time
%required to reach convergence
%
%Ted Amdur
%6/17/21, updated 7/23/21

clearvars
fsize = 16;
burnin = 200; %Number of burn-in reps assumed for chain length analysis
load observations_22_01_17.mat
rng(101)
pN = sIn.pN; %Number of realizations
reps = sIn.reps; %Number of Gibbs samples in chain
N = sIn.N;  %Length of dataset
chains = sIn.chains; %Number of chains
nObs = sIn.nObs;


%Goal: calculate the mixing ratios for the different variables of interest
%For each parameter, produce a single value giving the mixing ratio

%Cut out the burnin 
aC = xburnin(a,burnin, reps, chains, sIn.pN);
sigYC = xburnin(sigY,burnin, reps, chains, pN,nObs);
sigXC = xburnin(sigX,burnin, reps, chains, pN,5);
sP = xburnin(squeeze(A(:,:,2,:)),burnin, reps, chains, nObs);
bP = xburnin(squeeze(A(:,:,1,:)),burnin, reps, chains, pN,nObs);
mP = xburnin(squeeze(A(:,:,3,:)),burnin, reps, chains, pN,nObs);
aP = aC(:);
qP = sigXC';
rP = squeeze(sigYC)';

%FIRST: Find the mixing ratio of the AR(1) coefficient
[av,~,aw] = mpv(aP,chains,reps-burnin,1,reps-burnin);
ra = sqrt(av/aw);
aa = corr(aP(2:end),aP(1:end-1));

%SECOND: Find the mixing ratio for all the parameters in A
Av = NaN(size(A,[2 3])); Aw = NaN(size(A,[2 3]));
rA = NaN(size(A,[2 3]));
for iX = 1:size(A,2)
    %First get the ratios for the offsets (column 1)
    [v,~,w] = mpv(bP(iX,:),chains,reps-burnin,1,reps-burnin);
    Av(iX,1) = v; Aw(iX,1) = w;
    rA(iX,1) = sqrt(v/w);
    aA(iX,1) = corr(bP(iX,2:end)',bP(iX,1:end-1)');
    
    %Then get the ratios for the scalings (column 2)
    [v,~,w] = mpv(sP(iX,:),chains,reps-burnin,1,reps-burnin);
    Av(iX,2) = v; Aw(iX,2) = w;
    rA(iX,2) = sqrt(v/w);
    aA(iX,2) = corr(sP(iX,2:end)',sP(iX,1:end-1)');
    
    %Last, get the ratios for the linear drift parameters (column 3)
    [v,~,w] = mpv(mP(iX,:),chains,reps-burnin,1,reps-burnin);
    Av(iX,3) = v; Aw(iX,3) = w;
    rA(iX,3) = sqrt(v/w);
    aA(iX,3) = corr(mP(iX,2:end)',mP(iX,1:end-1)');
end

%THIRD: Find the mixing ratio for all the parameters for sigX
for iN = 1:size(qP,2)
    [qv,~,qw] = mpv(qP(:,iN),chains,reps-burnin,1,reps-burnin);
    rQ(iN) = sqrt(qv/qw);
    aQ(iN) = corr(qP(2:end,iN),qP(1:end-1,iN));
end

%FOURTH: Find the mixing ratio for all the parameters for sigY
for iN = 1:size(rP,2)
    [rv,~,rw] = mpv(rP(:,iN),chains,reps-burnin,1,reps-burnin);
    rR(iN) = sqrt(rv/rw);
    aR(iN) = corr(rP(2:end,iN),rP(1:end-1,iN));
end

    

function b = getb(p,chains,reps,c,c2)
    cp = reshape(p,[reps, chains]);
    cp = cp(c+1:c2,:);
    psij = mean(cp,1);
    psiall = mean(psij);
    b = ((c2-c)/(chains-1)).*(psij-psiall)*(psij-psiall)';
end
function w = getw(p,chains,reps,c,c2)
    cp = reshape(p,[reps, chains]);
    cp = cp(c+1:c2,:);
    psij = mean(cp,1);
    for ii = 1:chains
        sj2(ii) = (1/((c2-c)-1))*(cp(:,ii)-psij(ii))'*(cp(:,ii)-psij(ii));
    end
    w = inv(chains)*sum(sj2);
end
function [vout, b, w] = mpv(p,chains,reps,c,c2)
%Return the values needed to calculate 
%INPUTS:
%       p: posterior chain being input
%       chains: number of chains
%       reps: number of reps in each chain
%       c: index in an individual chain to begin observation
%       c2: index in an individual chain to end observation

%Get the relative variance between chains versus within chain
    b = getb(p,chains,reps,c,c2);
    w = getw(p,chains,reps,c,c2);
    vout = (c2-c-1)*w/(c2-c) + inv(c2-c)*b;
end
function out = xburnin(in,c, reps, chains, pN,numObs)
    %Revised as of 12/18/21 to reflect modified sigXchain structure
    if nargin < 6
        out = reshape(in,[pN,reps, chains]);
        %Remove the burn-in period from the variable record of interest
        out = out(:,c+1:end,:);
        out = reshape(out,[pN, (chains*reps-c*chains)]);
    else
        out = squeeze(reshape(in,[pN,numObs,reps, chains]));
        %Remove the burn-in period from the variable record of interest
        if ndims(out) > 3
            out = out(:,:,c+1:end,:);
        else
            out = out(:,c+1:end,:);
        end
        out = squeeze(reshape(out,[pN, numObs, (chains*reps-c*chains)]));
    end
end

