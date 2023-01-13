function [valE,oE,iX]=exciseobs(valM,oM,col,fracEx)

len=sum(oM(:,col));
oI=find(oM(:,col));
rI=ceil(rand.*len);
nani=mod(rI:rI+ceil(fracEx.*len),len);
nani=nani(nani~=0);
iX=oI(nani);
valM(iX,col)=NaN;
oM(iX,col)=false;
valE=valM;
oE=oM;
end