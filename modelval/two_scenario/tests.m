%Sandbox for determining how to correctly infer autocorrelation from time
%series data

for ii=1:1000
    aa=0;
    bb=0;
    cc=aa+bb;
    aMat(1,ii)=aa;
    bMat(1,ii)=bb;
    cMat(1,ii)=cc;
    for iB=2:50
        aMat(iB,ii)=0.95.*aMat(iB-1,ii)+randn;
        bMat(iB,ii)=0.99.*bMat(iB-1,ii)+randn;
        cMat(iB,ii)=aMat(iB,ii)+bMat(iB,ii);
    end
end

for ii=1:1000
    a(ii)=aMat(2:end,ii)'*aMat(1:end-1,ii)/(aMat(1:end-1,ii)'*aMat(1:end-1,ii));
    b(ii)=bMat(2:end,ii)'*bMat(1:end-1,ii)/(bMat(1:end-1,ii)'*bMat(1:end-1,ii));
    c(ii)=cMat(2:end,ii)'*cMat(1:end-1,ii)/(cMat(1:end-1,ii)'*cMat(1:end-1,ii));
    acorr(ii)=corr(aMat(2:end,ii),aMat(1:end-1,ii));
    bcorr(ii)=corr(bMat(2:end,ii),bMat(1:end-1,ii));
    ccorr(ii)=corr(cMat(2:end,ii),cMat(1:end-1,ii));
end
    