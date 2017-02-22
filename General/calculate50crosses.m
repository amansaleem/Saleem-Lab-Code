function [p50] = calculate50crosses(in, showplot)

if nargin<2
    showplot = false;
end

for n = 1:size(in,1)
    nIn(n,:) = normalise1var(in(n,:));
end
snIn = sort(nIn,2);
for n = 1:size(in,1)
    if sum(~isnan(snIn(n,:)))>3
        p50(n) = min(find(snIn(n,:)>0.5));
    else
        p50(n) = nan;
    end
end

if showplot
    [~, sidx] = sort(p50);
    subplot(121)
    imagesc(nIn); colorbar;
    subplot(122)
    imagesc(snIn(sidx,:)); colorbar;
    hold on
    plot(p50(sidx),1:size(in,1),'k.');
end
