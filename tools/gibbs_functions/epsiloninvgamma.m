function [rmat,theta] = epsiloninvgamma(T0,th0,infI,errorN,NN)
%epsiloninvgamma.m calculates the variance of the observation error for each observer 
%using inverse gamma distribution
%Inputs:
% T0: vector of prior inverse gamma shape parameters for each observer
% th0: vector of prior inverse gamma scale parameters for each observer
% infI: logical matrix indicating observations that are included in the analysis
% NN: number of observers
%
%Outputs:
% rmat: variance of observation error for each observer
% theta: sum of squared error for each observer

rmat=zeros(NN,1); %initialize the output variable for variances of observation error
theta=zeros(NN,1);
for ii=1:NN %loop over all observers
    %calculate the variance of observation error using inverse gamma distribution
    rmat(ii) = IG(T0(ii), th0(ii), errorN(infI(:,ii),ii));
    %calculate the sum of squared error for each observer
    theta(ii) = errorN(infI(:,ii),ii)' * errorN(infI(:,ii),ii);
end

end


