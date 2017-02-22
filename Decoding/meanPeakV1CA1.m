close all;
clear;

% load Posterior_all_150226_trainedCV_125_noInter.mat
% % load Posterior_all_150226_trainedCV_thetaBins_noInter.mat
% % load Posterior_all_150624_trainedCV_thetaBins_noInterall_new

% load Posterior_all_V1_150601_trainedCV_250
load Posterior_all_151014_V1_trainedCV_125
samp_rate = 60;
start_lim = 100;
getMeanPeakCentre;
peakRate_V1 = peakRate;
meanRate_V1 = meanRate;
mRate_V1 = mRate;
pRate_V1 = pRate;
lRate_V1 = lRate;

clear meanRate peakRate mRate pRate
% load Posterior_all_150723_trainedCV_thetaBins_FineBin_EVover50
% load Posterior_all_150706_trainedCV_thetaBins_FineBin1000_FixedSmth20
load Posterior_all_150226_trainedCV_thetaBins_noInter
samp_rate = 8;
start_lim = 100;
getMeanPeakCentre;
peakRate_CA1 = peakRate;
meanRate_CA1 = meanRate;
mRate_CA1 = mRate;
pRate_CA1 = pRate;
lRate_CA1 = lRate;

contrasts = [0 0.18 0.6 0.72];
%% 
figure(1)
subplot(131)
errorbarxy(pRate_CA1.l,pRate_CA1.h,...
    pRate_CA1.l_s,pRate_CA1.h_s,{'ko','k','k'});
% hist(2*(peakRate_CA1.h-peakRate_CA1.l)./(peakRate_CA1.h+peakRate_CA1.l),30)
% axis equal; axis([0 5 0 5]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% line(xlim, ylim,'linestyle','--','color','k');
% axis equal; axis([0.75 1.3 0.75 1.3]);
% line(xlim, [1 1],'linestyle','--','color','k'); line([1 1],  ylim,'linestyle','--','color','k');
axis equal
axis([0 35 0 35])
line(xlim, ylim,'linestyle','--','color','k');
xlabel('Peak rate CA1, low')
ylabel('Peak rate CA1, high')

subplot(132)
errorbarxy(pRate_CA1.ld,pRate_CA1.hd,...
    pRate_CA1.ld_s,pRate_CA1.hd_s,{'ko','k','k'});
% hist(2*(peakRate_CA1.hd-peakRate_CA1.ld)./(peakRate_CA1.hd+peakRate_CA1.ld),30)
% axis equal; axis([0 5 0 5]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis equal
axis([0 35 0 35])
line(xlim, ylim,'linestyle','--','color','k');% line(xlim, ylim,'linestyle','--','color','k');
% axis equal; axis([0.75 1.3 0.75 1.3]);
% line(xlim, [1 1],'linestyle','--','color','k'); line([1 1],  ylim,'linestyle','--','color','k');
xlabel('Peak rate CA1 dec, low')
ylabel('Peak rate CA1 dec, high')

subplot(133)
errorbarxy(pRate_CA1.l_low,pRate_CA1.ld_low,...
    pRate_CA1.l_low_s,pRate_CA1.ld_low_s,{'ko','k','k'});
% hist(2*(peakRate_CA1.ld_low-peakRate_CA1.l_low)./(peakRate_CA1.ld_low+peakRate_CA1.l_low),30)
%     axis equal; axis([0 5 0 5]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis equal
axis([0 35 0 35])
line(xlim, ylim,'linestyle','--','color','k');% line(xlim, ylim,'linestyle','--','color','k');
% axis equal; axis([0.75 1.3 0.75 1.3]);
% line(xlim, [1 1],'linestyle','--','color','k'); line([1 1],  ylim,'linestyle','--','color','k');
xlabel('Peak rate CA1, low')
ylabel('Peak rate CA1, low dec')

%%
figure(2)
hold off
% subplot(121);
% rate_ratio_CA1 = [mRate_CA1.l'./mRate_CA1.n' mRate_CA1.n'./mRate_CA1.n' mRate_CA1.h'./mRate_CA1.n']';
rate_ratio_CA1 = [mRate_CA1.l' mRate_CA1.n' mRate_CA1.h']';
% rate_ratio_CA1 = [meanRate_CA1.l' meanRate_CA1.n' meanRate_CA1.h']';
% rate_ratio_CA1 = [meanRate_CA1.l'./meanRate_CA1.n' meanRate_CA1.n'./meanRate_CA1.n' meanRate_CA1.h'./meanRate_CA1.n']';
bar([ 1 2 3],mean(rate_ratio_CA1'))
hold on;
errorbar([ 1 2 3], mean(rate_ratio_CA1'), sem(rate_ratio_CA1'),'o-');
% errorbar([1 2 3], mean(rate_ratio_CA1'),'ro-','linewidth',3)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
line(xlim, [1 1], 'linestyle','--','color','k');

% subplot(122);
% rate_ratio_V1 = [mRate_V1.l'./mRate_V1.n' mRate_V1.n'./mRate_V1.n' mRate_V1.h'./mRate_V1.n']';
rate_ratio_V1 = [mRate_V1.l' mRate_V1.n' mRate_V1.h']';
% rate_ratio_V1 = [meanRate_V1.l' meanRate_V1.n' meanRate_V1.h']';
% rate_ratio_V1 = [meanRate_V1.l'./meanRate_V1.n' meanRate_V1.n'./meanRate_V1.n' meanRate_V1.h'./meanRate_V1.n']';
bar([4 5 6],nanmean(rate_ratio_V1'))
hold on;
errorbar([4 5 6], nanmean(rate_ratio_V1'), nansem(rate_ratio_V1'),'o-');
% errorbar([4 5 6], mean(rate_ratio_V1'),'ko-','linewidth',3)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
line(xlim, [1 1], 'linestyle','--','color','k');

set(gca,'YLim',[0.80 1.1])


%%
figure(3)
subplot(131)
errorbarxy(pRate_V1.l,pRate_V1.h,...
    pRate_V1.l_s,pRate_V1.h_s,{'ko','k','k'});
% hist(2*(peakRate_V1.h-peakRate_V1.l)./(peakRate_V1.h+peakRate_V1.l),30)
% axis equal; axis([0 5 0 5]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% line(xlim, ylim,'linestyle','--','color','k');
% axis equal; axis([0.75 1.3 0.75 1.3]);
% line(xlim, [1 1],'linestyle','--','color','k'); line([1 1],  ylim,'linestyle','--','color','k');
axis equal
axis([0 35 0 35])
line(xlim, ylim,'linestyle','--','color','k');
xlabel('Peak rate V1, low')
ylabel('Peak rate V1, high')

subplot(132)
errorbarxy(pRate_V1.ld,pRate_V1.hd,...
    pRate_V1.ld_s,pRate_V1.hd_s,{'ko','k','k'});
% hist(2*(peakRate_V1.hd-peakRate_V1.ld)./(peakRate_V1.hd+peakRate_V1.ld),30)
% axis equal; axis([0 5 0 5]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis equal
axis([0 35 0 35])
line(xlim, ylim,'linestyle','--','color','k');% line(xlim, ylim,'linestyle','--','color','k');
% axis equal; axis([0.75 1.3 0.75 1.3]);
% line(xlim, [1 1],'linestyle','--','color','k'); line([1 1],  ylim,'linestyle','--','color','k');
xlabel('Peak rate V1 dec, low')
ylabel('Peak rate V1 dec, high')

subplot(133)
errorbarxy(pRate_V1.l_low,pRate_V1.ld_low,...
    pRate_V1.l_low_s,pRate_V1.ld_low_s,{'ko','k','k'});
% hist(2*(peakRate_V1.ld_low-peakRate_V1.l_low)./(peakRate_V1.ld_low+peakRate_V1.l_low),30)
%     axis equal; axis([0 5 0 5]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis equal
axis([0 35 0 35])
line(xlim, ylim,'linestyle','--','color','k');% line(xlim, ylim,'linestyle','--','color','k');
% axis equal; axis([0.75 1.3 0.75 1.3]);
% line(xlim, [1 1],'linestyle','--','color','k'); line([1 1],  ylim,'linestyle','--','color','k');
xlabel('Peak rate V1, low')
ylabel('Peak rate V1, low dec')

% %%
% figure(1)
% hold off;
% % subplot(121);
% rate_ratio_CA1 = [mRate_CA1.l'./mRate_CA1.n' mRate_CA1.n'./mRate_CA1.n' mRate_CA1.h'./mRate_CA1.n']';
% plot([4 5 6], rate_ratio_CA1,'o-','color',[1 .5 .5],'linewidth',1);
% hold on;
% plot([4 5 6], mean(rate_ratio_CA1'),'ro-','linewidth',3)
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% line(xlim, [1 1], 'linestyle','--','color','k');
% 
% % subplot(122);
% rate_ratio_V1 = [mRate_V1.l'./mRate_V1.n' mRate_V1.n'./mRate_V1.n' mRate_V1.h'./mRate_V1.n']';
% plot([1 2 3], rate_ratio_V1,'o-','color',[.5 .5 .5],'linewidth',1);
% hold on;
% plot([1 2 3], mean(rate_ratio_V1'),'ko-','linewidth',3)
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% line(xlim, [1 1], 'linestyle','--','color','k');
% %%
% peakRate_CA1.l = peakRate_CA1.l   -   peakRate_CA1.n;
% peakRate_CA1.h = peakRate_CA1.h   -   peakRate_CA1.n;
% peakRate_CA1.n = peakRate_CA1.n   -   peakRate_CA1.n;
% 
% peakRate_V1.l = peakRate_V1.l   -   peakRate_V1.n;
% peakRate_V1.h = peakRate_V1.h   -   peakRate_V1.n;
% peakRate_V1.n = peakRate_V1.n   -   peakRate_V1.n;
% 
% meanRate_CA1.l = meanRate_CA1.l   -   meanRate_CA1.n;
% meanRate_CA1.g = meanRate_CA1.g   -   meanRate_CA1.n;
% meanRate_CA1.h = meanRate_CA1.h   -   meanRate_CA1.n;
% meanRate_CA1.n = meanRate_CA1.n   -   meanRate_CA1.n;
% 
% meanRate_V1.l = meanRate_V1.l   -   meanRate_V1.n;
% meanRate_V1.g = meanRate_V1.g   -   meanRate_V1.n;
% meanRate_V1.h = meanRate_V1.h   -   meanRate_V1.n;
% meanRate_V1.n = meanRate_V1.n   -   meanRate_V1.n;
% 
% cc    = contrasts;
% rr_V1 = [nanmean(meanRate_V1.g) nanmean(meanRate_V1.l) nanmean(meanRate_V1.n) nanmean(meanRate_V1.h)];%[meanRate_V1.g' meanRate_V1.l' meanRate_V1.n' meanRate_V1.h'];
% cc_CA1    = repmat(contrasts(2:4), size(meanRate_CA1.l'));
% rr_CA1 = [(meanRate_CA1.l)' (meanRate_CA1.n)' (meanRate_CA1.h)'];%[meanRate_CA1.l' meanRate_CA1.n' meanRate_CA1.h'];
% [err_V1, pars_V1] = fit_hyper_ratio(cc(2:4),rr_V1(2:4));
% % [err_CA1, pars_CA1] = fit_hyper_ratio(cc(2:4),rr_CA1);
% 
% %%
% [xData, yData] = prepareCurveData( cc_CA1, rr_CA1 );
% ft = fittype( 'poly1' );
% opts = fitoptions( ft );
% opts.Lower = [-Inf -Inf];
% opts.Robust = 'LAR';
% opts.Upper = [Inf Inf];
% % Fit model to data.
% [fitresult, gof] = fit( xData, yData, ft, opts );
% 
% %%
% V1mean_mean = [mean(meanRate_V1.g) mean(meanRate_V1.l) mean(meanRate_V1.n) mean(meanRate_V1.h)];
% V1sems_mean = [sem(meanRate_V1.g) sem(meanRate_V1.l) sem(meanRate_V1.n) sem(meanRate_V1.h)];
% V1mean_peak = [mean(peakRate_V1.l) mean(peakRate_V1.n) mean(peakRate_V1.h)];
% V1sems_peak = [sem(peakRate_V1.l) sem(peakRate_V1.n) sem(peakRate_V1.h)];
% 
% CA1mean_mean = [mean(meanRate_CA1.g) mean(meanRate_CA1.l) mean(meanRate_CA1.n) mean(meanRate_CA1.h)];
% CA1sems_mean = [sem(meanRate_CA1.g) sem(meanRate_CA1.l) sem(meanRate_CA1.n) sem(meanRate_CA1.h)];
% CA1mean_peak = [mean(peakRate_CA1.l) mean(peakRate_CA1.n) mean(peakRate_CA1.h)];
% CA1sems_peak = [sem(peakRate_CA1.l) sem(peakRate_CA1.n) sem(peakRate_CA1.h)];
% 
% figure;
% subplot(121)
% % bar(contrasts([2:4])-0.15, V1mean_mean([2:4]),0.25,'k');
% % hold on;
% % bar(contrasts([2:4])+0.15, CA1mean_mean([2:4]),0.25,'r');
% X = 0.15:0.01:0.75;
% errorbar(contrasts([2:4]), V1mean_mean([2:4]), V1sems_mean([2:4]),'k');
% hold on;
% errorbar(contrasts([2:4]), CA1mean_mean([2:4]), CA1sems_mean([2:4]),'r');
% plot(X, hyper_ratio(pars_V1, X),'k');
% plot(X, fitresult.p2 + fitresult.p1*X,'r');
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('Mean Rate');
% xlabel('Contrast');
% ylabel('Mean rate - mean_{medium contrast} (Spikes/s)');
% % axis tight
% set(gca, 'XLim', [-0.5 0.8]);
% 
% subplot(122)
% % bar(contrasts([2 4])-0.15, V1mean_peak([1 3]),0.25,'k');
% % hold on;
% % bar(contrasts([2 4])+0.15, CA1mean_peak([1 3]),0.25,'r');
% errorbar(contrasts([2:4]), V1mean_peak([1:3]), V1sems_peak([1:3]),'k');
% hold on;
% errorbar(contrasts([2:4]), CA1mean_peak([1:3]), CA1sems_peak([1:3]),'r','linewidth',2);
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('Peak Rate')
% xlabel('Contrast');
% ylabel('Peak rate - peak_{medium contrast} (Spikes/s)');
% % axis tight
% set(gca, 'XLim', [-0.5 0.8]);
% % 
% % subplot(223)
% % % bar(contrasts([2:4])-0.15, V1mean_mean([2:4]),0.25,'k');
% % % hold on;
% % bar(contrasts([2:4])+0.15, CA1mean_mean([2:4]),0.25,'r');
% % % errorbar(contrasts([2:4])-0.15, V1mean_mean([2:4]), V1sems_mean([2:4]),'k','linewidth',2);
% % hold on;
% % set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('Mean Rate');
% % xlabel('Contrast');
% % ylabel('Mean rate - mean_{medium contrast} (Spikes/s)');
% % axis tight
% % set(gca, 'XLim', [-1 3.5]);
% % 
% % subplot(224)
% % % bar(contrasts([2 4])-0.15, V1mean_peak([1 3]),0.25,'k');
% % % hold on;
% % bar(contrasts([2 4])+0.15, CA1mean_peak([1 3]),0.25,'r');
% % % errorbar(contrasts([2 4])-0.15, V1mean_peak([1 3]), V1sems_peak([1 3]),'k','linewidth',2);
% % hold on;
% % errorbar(contrasts([2 4])+0.15, CA1mean_peak([1 3]), CA1sems_peak([1 3]),'r','linewidth',2);
% % set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% % title('Peak Rate')
% % xlabel('Contrast');
% % ylabel('Peak rate - peak_{medium contrast} (Spikes/s)');
% % axis tight
% % set(gca, 'XLim', [-0.5 3.5]);
% 
% %%
% figure;
% [pC, pX] =  hist(-peakRate_CA1.l ,-7:1:7);
% [pV, pX] =  hist(-peakRate_V1.l ,-7:1:7);
% pC = pC  ./  sum(pC);
% pV = pV  ./ sum(pV);
% 
% [mC, mX] =  hist(-meanRate_CA1.l , -5:1:5);
% [mV, mX] =  hist(-meanRate_V1.l , -5:1:5);
% mV = mV  ./  sum(mV);
% mC = mC  ./  sum(mC);
% 
% subplot(121)
% plot(mX, mV, 'bo-', mX, mC, 'ro-');
% set(gca, 'YLim',[0 0.35]);
% line([0 0], ylim, 'linestyle','--', 'color','k')
% line(mean(-meanRate_CA1.l)*[1 1], ylim, 'color','r')
% line(mean(-meanRate_V1.l)*[1 1], ylim, 'color','b')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% set(gca, 'XTick',[-4 -2 0 2 4], 'YTick', [0 0.1 0.2 0.3]);
% xlabel('Difference in mean rates (spikes/sec)');
% ylabel('Fraction of cells')
% 
% subplot(122)
% plot(pX, pV, 'bo-', pX, pC, 'ro-');
% set(gca, 'YLim',[0 0.35]);
% line([0 0], ylim, 'linestyle','--', 'color','k')
% line(mean(-peakRate_CA1.l)*[1 1], ylim, 'color','r')
% line(mean(-peakRate_V1.l)*[1 1], ylim, 'color','b')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% set(gca, 'XTick',[-4 -2 0 2 4], 'YTick', [0 0.1 0.2 0.3]);
% xlabel('Difference in peak rates (spikes/sec)');
% ylabel('Fraction of cells')
% 
% %%
% figure;
% % subplot(121)
% cc2 = [1 2];
% errorbar(cc2, V1mean_mean([2 4]), V1sems_mean([2 4]),'k');
% hold on;
% errorbar(cc2+3, CA1mean_mean([2 4]), CA1sems_mean([2 4]),'r','linewidth',2);
% bar(cc2, V1mean_mean([2 4]),0.25,'k');
% bar(cc2+3, CA1mean_mean([2 4]),0.25,'r');
% hold on;
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('Mean Rate');
% xlabel('Contrast');
% ylabel('Mean rate - mean_{medium contrast} (Spikes/s)');
% axis tight
% set(gca, 'XLim', [-1 5.5]);
% 
% figure;
% % subplot(121)
% cc2 = [1 2];
% errorbar(cc2, V1mean_peak([1 3]), V1sems_peak([1 3]),'k');
% hold on;
% errorbar(cc2+3, CA1mean_peak([1 3]), CA1sems_peak([1 3]),'r','linewidth',2);
% bar(cc2, V1mean_peak([1 3]),0.25,'k');
% bar(cc2+3, CA1mean_peak([1 3]),0.25,'r');
% hold on;
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% title('Peak Rate');
% xlabel('Contrast');
% ylabel('Peak rate - peak_{medium contrast} (Spikes/s)');
% axis tight
% set(gca, 'XLim', [-1 5.5]);