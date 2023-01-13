function [trends,offsets] = returntrend(A,t,ind)
reps = size(A,3);
offsets = NaN(size(A,1),reps);
trends = NaN([size(t,1),size(A,1),reps]);
for ii = 1:reps
    offsets(:,ii) = squeeze(A(:,1,ii));
    trends(:,:,ii) = squeeze(A(:,3,ii)').*t;
end
if nargin == 3
    trends=squeeze(trends(:,ind,:));
    offsets=offsets(ind,:);
end
end
