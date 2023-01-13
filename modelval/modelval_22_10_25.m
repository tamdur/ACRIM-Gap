%Evaluate the predicted observer characteristics by excising half of each
%observer's record, seeing if the result can be consistently predicted
%Ted Amdur
%10/25/22


%Twelve columns of valM correspond to the following observers:
%     "ACRIM1/SMM"
%     "ACRIM2/UARS"
%     "ACRIM3"
%     "BremenMgII"
%     "ERBE/ERBS"
%     "NIMBUS-7/HF"
%     "PREMOS/PICARD"
%     "SILSO"
%     "SORCE/TIM"
%     "TCTE"
%     "TSIS-1"
%     "VIRGO/SOHO"
tic
rng(7);
obsTest=12; %Observer to run experiment for


obsmatrix='obs_22_7_12';
load(obsmatrix); %From makeobsmatrix.m
load excludeMask22_6_29.mat %from exclude_fliers_22_04_26.m
valM(excludeMask) = NaN;
oM(excludeMask) = false;
col=12;
fracEx=0.5;
[valE,oE,iX]=exciseobs(valM,oM,col,fracEx);


[obsPrior] = estimateobspriors(valE,oE,colLabels,1);
[A,xAll,sigY,t] = runevalchain_22_10_25(valE,oE,colLabels,obsPrior);

toc
plotpredictions(A,xAll,t,valM,iX,col,sigY,offsets,dateM);