function plotpredictions(A,xAll,t,valM,iX,col,sigY,offsets,dateM)
reps = size(A,3);
len=length(iX);
pO=zeros(len,reps);
res=zeros(len,reps);
xm=mean(xAll,2);
for ii=1:reps
    pred=[ones(len,1) xAll(iX,ii) t(iX,col)]; %make predictor matrix
    p=squeeze(A(col,:,ii))*pred'; %Predict from observation model
    p=p';
    p=p+randn(len,1).*sqrt(sigY(col,ii));
    res(:,ii)=p-xm(iX);
    p=p;
    pO(:,ii)=p;
end

figure
plot(dateM(iX),prctile(res',[2.5 10 50 90 97.5]),'Color','k')
hold on
plot(dateM(iX),prctile(res',50),'Color','k','LineWidth',3)
hold on
plot(dateM,valM(:,col)-xm,'.','Color','r','MarkerSize',20)
xlabel('Year')
ylabel('TSI (W/m^{2})')
set(gca,'FontSize',16)
end

