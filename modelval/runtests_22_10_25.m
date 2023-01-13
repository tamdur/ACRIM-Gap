%Run tests of observer predictions in parallel computing environment
% Ted Amdur
% 10/25/22


parpool('local',str2num(getenv('SLURM_CPUS_PER_TASK')))
rng(1)
tN=100; %Number of experiments per observer
obsmatrix='obs_22_7_12.mat';
load(obsmatrix); %From makeobsmatrix.m
load excludeMask22_6_29.mat %from exclude_fliers_22_04_26.m
valM(excludeMask) = NaN;
oM(excludeMask) = false;
fracEx=0.5; %Exclude this fraction of results
reps=1000; %Length of kept chain
parfor col=1:length(colLabels)
    tic;
    oMi=oM;
    len=ceil(fracEx.*sum(oMi(:,col)))+1;
    iXM=NaN(tN,len);
    AM=zeros(tN,3,reps);
    sigYM=zeros(tN,reps);
    pM=NaN(tN,len,reps);
    xM=NaN(tN,len,reps);
    for ii=1:tN
        [valE,oE,iX]=exciseobs(valM,oMi,col,fracEx);
        len=length(iX);
        [obsPrior] = estimateobspriors(valE,oE,colLabels,1);
        [A,xAll,sigY,t] = runevalchain_22_10_25(valE,oE,colLabels,obsPrior);
         pred=[ones(len,1) xAll(iX,ii) t(iX,col)]; %make predictor matrix
         pO=zeros(len,reps);
         for iN=1:reps
             pred=[ones(len,1) xAll(iX,iN) t(iX,col)]; %make predictor matrix
             p=squeeze(A(col,:,iN))*pred'; %Predict from observation model
             p=p';
             p=p+randn(len,1).*sqrt(sigY(col,iN));
             pO(:,iN)=p;
         end
         iXM(ii,1:len)=iX;
         AM(ii,:,:)=squeeze(A(col,:,:));
         sigYM(ii,:)=sigY(col,:);
         pM(ii,1:len,:)=pO;
         xM(ii,1:len,:)=xAll(iX,:);
    end
    oTest(col).iX=iXM;
    oTest(col).A=AM;
    oTest(col).sigY=sigYM;
    oTest(col).p=pM;
    oTest(col).x=xM;
    oTest(col).tRun=toc;
end
save('test_22_10_25.mat','oTest')


