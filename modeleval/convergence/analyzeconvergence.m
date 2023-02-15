%Generate multiple chains using randomized starting parameters, analyze
%chains for intra and inter chain mixing using the potential scale
%reduction metric \hat{R}.
%
%Ted Amdur
%2023/02/02, based off analyze convergence.m in code_21_11_20 directory

clearvars
rng(1) %Set random number generator for reproducability
fsize = 16;
load obs_23_02_01.mat
opts.burnin = 500; %Number of burn-in reps assumed for chain length analysis
opts.chain=10500; %Total length of chain, including burn-in
opts.randomizeChain=true; %Disperse chains for independent comparison
opts.dispProgress=true;

%First, run chain three times for comparison
opts.saveFile='ar2_chain1_23_02_08.mat';
runchain_23_01_13(valM,oM,colLabels,opts);
opts.saveFile='ar2_chain2_23_02_08.mat';
runchain_23_01_13(valM,oM,colLabels,opts);
opts.saveFile='ar2_chain3_23_02_08.mat';
runchain_23_01_13(valM,oM,colLabels,opts);

%%%%%%%%%%%%%%%
%Goal: calculate the mixing ratios for the different variables of interest
%For each parameter, produce a single value giving the mixing ratio

%Select variable of interest 
chain=struct;
variableName='xAll(50,:)';
load ar2_chain1_23_02_04.mat
eval(['chain(1).x=',variableName,';']);
load ar2_chain2_23_02_04.mat
eval(['chain(2).x=',variableName,';']);
load ar2_chain3_23_02_04.mat
eval(['chain(3).x=',variableName,';']);

reps=opts.chain-opts.burnin; %Number of reps
chainN=3; %Number of chains
ind=1;
for ii=1:3
    chainMat(:,ind)=chain(ii).x(1:floor(reps/2)); %chainMat in dimension repsX(chainsx2)
    ind=ind+1;
    chainMat(:,ind)=chain(ii).x(floor(reps/2)+1:end);
    ind=ind+1;
end
% [vout, b, w] = mpv(chainMat,size(chainMat,2),size(chainMat,1),1,size(chainMat,1));
[vout, b, w] = mpv(chainMat,size(chainMat,2),size(chainMat,1),1);
Rhat=rhat(vout,w);



% %FIRST: Find the mixing ratio of the AR(1) coefficient
% [av,~,aw] = mpv(aP,chains,reps-burnin,1,reps-burnin);
% ra = sqrt(av/aw);
% aa = corr(aP(2:end),aP(1:end-1));
% 
% %SECOND: Find the mixing ratio for all the parameters in A
% Av = NaN(size(A,[2 3])); Aw = NaN(size(A,[2 3]));
% rA = NaN(size(A,[2 3]));
% for iX = 1:size(A,2)
%     %First get the ratios for the offsets (column 1)
%     [v,~,w] = mpv(bP(iX,:),chains,reps-burnin,1,reps-burnin);
%     Av(iX,1) = v; Aw(iX,1) = w;
%     rA(iX,1) = sqrt(v/w);
%     aA(iX,1) = corr(bP(iX,2:end)',bP(iX,1:end-1)');
%     
%     %Then get the ratios for the scalings (column 2)
%     [v,~,w] = mpv(sP(iX,:),chains,reps-burnin,1,reps-burnin);
%     Av(iX,2) = v; Aw(iX,2) = w;
%     rA(iX,2) = sqrt(v/w);
%     aA(iX,2) = corr(sP(iX,2:end)',sP(iX,1:end-1)');
%     
%     %Last, get the ratios for the linear drift parameters (column 3)
%     [v,~,w] = mpv(mP(iX,:),chains,reps-burnin,1,reps-burnin);
%     Av(iX,3) = v; Aw(iX,3) = w;
%     rA(iX,3) = sqrt(v/w);
%     aA(iX,3) = corr(mP(iX,2:end)',mP(iX,1:end-1)');
% end
% 
% %THIRD: Find the mixing ratio for all the parameters for sigX
% for iN = 1:size(qP,2)
%     [qv,~,qw] = mpv(qP(:,iN),chains,reps-burnin,1,reps-burnin);
%     rQ(iN) = sqrt(qv/qw);
%     aQ(iN) = corr(qP(2:end,iN),qP(1:end-1,iN));
% end
% 
% %FOURTH: Find the mixing ratio for all the parameters for sigY
% for iN = 1:size(rP,2)
%     [rv,~,rw] = mpv(rP(:,iN),chains,reps-burnin,1,reps-burnin);
%     rR(iN) = sqrt(rv/rw);
%     aR(iN) = corr(rP(2:end,iN),rP(1:end-1,iN));
% end

    

function b = getb(cp,m,reps,c)
%Estimate the between-sequence variance, from pg. 284 of BDA by Gelman
cp=cp(c:c+reps-1,:);
n=reps;
psij = mean(cp,1);
psiall = mean(psij);
b = (n/(m-1)).*(psij-psiall)*(psij-psiall)';
end
function w = getw(cp,m,reps,c)
%Estimate the within-sequence variance, from pg. 284 of BDA by Gelman
cp=cp(c:c+reps-1,:);
n=reps;
psij = mean(cp,1);
for ii = 1:m
    sj2(ii) = (1/(n-1))*(cp(:,ii)-psij(ii))'*(cp(:,ii)-psij(ii));
end
w = (1/m)*sum(sj2);
end
function [vout, b, w] = mpv(p,chains,reps,c,c2)
%Return the values needed to calculate the marginal posterior variance of
%the estimand
%INPUTS:
%       p: posterior chain being input
%       chains: number of chains
%       reps: number of reps in each chain
%       c: index in an individual chain to begin observation
%       c2: index in an individual chain to end observation

%Get the relative variance between chains versus within chain
    b = getb(p,chains,reps,c);
    w = getw(p,chains,reps,c);
    vout = (reps-1)*w/reps + (1/reps)*b;
end

function Rhat=rhat(mpv,W)
%Estimate the potential scale reduction, from pg. 285 of BDA by Gelman
    Rhat=sqrt(mpv./W);
end

