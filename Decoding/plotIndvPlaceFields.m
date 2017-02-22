function plotIndvPlaceFields(P, animalIdx, icell)
bins = P(animalIdx).decoder.bins;

ax(1) = subplot(121);
hold off
plot(bins, 7.5*P(animalIdx).decoder_high.model.meanModel(icell,:),'r')
hold on;
plot(bins, 7.5*P(animalIdx).decoder_low.model.meanModel(icell,:),'b')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis tight;
ylims = ylim;
axis square
xlabel('Position')
ylabel('Firing rate (spikes/s)')
set(gca, 'XTick', [0 50 100]);

ax(2) = subplot(122);
hold off
plot(bins, 7.5*P(animalIdx).decoder_low_decPos.model.meanModel(icell,:),'b--')
hold on;
plot(bins, 7.5*P(animalIdx).decoder_low.model.meanModel(icell,:),'b')
% plot(bins, 7.5*P(animalIdx).decoder_high_decPos.model.meanModel(icell,:),'r--')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis tight;
axis square;
ymax = max(max(ylim),max(ylims));

xlabel('Position')
set(gca, 'XTick', [0 50 100]);
linkaxes(ax,'xy')
set(gca, 'xlim',[-0.001 100.001]);
set(gca, 'ylim',[0 ymax]);