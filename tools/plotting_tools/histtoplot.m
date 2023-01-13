function [x,y,hp] = histtoplot(counts,edges,xSmooth)
x = edges(1:end-1)+diff(edges)./2;
counts = smoothPH(counts,xSmooth);
hp = 1./sum(counts.*diff(edges));
y = counts.*hp;
end
