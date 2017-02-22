%% Getting all the data into a nice format
start = 33;
stop  = 37;
alpha_val = 0.2;
numBins = size(Posterior_all(1).Posterior_norm,2);
try
    stopthere;
meanPostTrend(1) = Posterior_predictions_safe(Posterior_530_smooth_250, start, stop);
meanPostTrend(2) = Posterior_predictions_safe(Posterior_531_smooth_250, start, stop);
meanPostTrend(3) = Posterior_predictions_safe(Posterior_601_smooth_250, start, stop);
meanPostTrend(4) = Posterior_predictions_safe(Posterior_602_smooth_250, start, stop);
meanPostTrend(5) = Posterior_predictions_safe(Posterior_603_smooth_250, start, stop);
meanPostTrend(6) = Posterior_predictions_safe(Posterior_1030_smooth_250, start, stop);
meanPostTrend(7) = Posterior_predictions_safe(Posterior_1025_smooth_250, start, stop);
% meanPostTrend(8) = Posterior_predictions_safe(Posterior_604_smooth_250, start, stop);
catch
    for idx = 1:length(Posterior_all)
        meanPostTrend(idx) = Posterior_predictions_safe(Posterior_all(idx), start, stop);
    end
end
%% Average across all the experimets, posterior in reward zone with position
close all
allMeanPostTrend.correct = [];
allMeanPostTrend.miss = [];
allMeanPostTrend.late = [];
allMeanPostTrend.wrong = [];
for exptIdx = 1:length(Posterior_all) %[1:3 5:8]
    normvals = minmax(nanmean(meanPostTrend(exptIdx).all.correct,1));
    allMeanPostTrend.correct =  [allMeanPostTrend.correct (meanPostTrend(exptIdx).all.correct'-normvals(1))./(normvals(2)-normvals(1))];
    allMeanPostTrend.miss    =  [allMeanPostTrend.miss (meanPostTrend(exptIdx).all.miss'-normvals(1))./(normvals(2)-normvals(1))];
    allMeanPostTrend.late    =  [allMeanPostTrend.late (meanPostTrend(exptIdx).all.late'-normvals(1))./(normvals(2)-normvals(1))];
    allMeanPostTrend.wrong   =  [allMeanPostTrend.wrong (meanPostTrend(exptIdx).all.wrong'-normvals(1))./(normvals(2)-normvals(1))];
    
%     allMeanPostTrend.correct =  [allMeanPostTrend.correct meanPostTrend(exptIdx).all.correct'];
%     allMeanPostTrend.miss =     [allMeanPostTrend.miss meanPostTrend(exptIdx).all.miss'];
%     allMeanPostTrend.wrong =    [allMeanPostTrend.wrong meanPostTrend(exptIdx).all.wrong'];
end
figure;
disp('Warning! check is this is valid on the data!');
allMeanPostTrend.late(allMeanPostTrend.late<0) = NaN;
allMeanPostTrend.correct(allMeanPostTrend.correct<0) = NaN;
allMeanPostTrend.wrong(allMeanPostTrend.wrong<0) = NaN;
allMeanPostTrend.correct(allMeanPostTrend.correct>5) = NaN;
allMeanPostTrend.wrong(allMeanPostTrend.wrong>5) = NaN;
try
    for n = 1:100
        [~, sig_early(n)] = ttest2(allMeanPostTrend.correct(n,:), allMeanPostTrend.wrong(n,:));
        [~, sig_late(n)] = ttest2(allMeanPostTrend.correct(n,:), allMeanPostTrend.late(n,:));
    end
    errorarea_as(1:100, nanmean(allMeanPostTrend.wrong,2)',...
        nansem(allMeanPostTrend.wrong'),'r', alpha_val)
    hold on;
    errorarea_as(1:100, nanmean(allMeanPostTrend.correct,2)',...
        nansem(allMeanPostTrend.correct'),'c', alpha_val);
    errorarea_as(1:100, nanmean(allMeanPostTrend.late,2)',...
        nansem(allMeanPostTrend.late'),'b', alpha_val)
    errorarea_as(1:100, nanmean(allMeanPostTrend.miss,2)',...
        nansem(allMeanPostTrend.miss'),'k', alpha_val)
    plot(find(sig_late <0.05),0.6,'*b')
    plot(find(sig_early<0.05),0.55,'*r')
catch
    plot(1:100, nanmean(allMeanPostTrend.wrong,2)','r')
    hold on;
    plot(1:100, nanmean(allMeanPostTrend.correct,2)','c');
    plot(1:100, nanmean(allMeanPostTrend.late,2)','b')
    plot(1:100, nanmean(allMeanPostTrend.miss,2)','k')
    
end
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
set(gca, 'XLim',[-5 105])
% set(gca, 'YLim', [-1.5 1.25]);
set(gca, 'YLim', [0 1.1]);
%     xlabel('Position (cm)')
%     ylabel(['Posterior in ' num2str(2*start) ' to ' num2str(2*stop)]);
line([70-5 70-5],ylim,'color','c');
line([70+5 70+5],ylim,'color','c');
% line(xlim,[0 0],'color','k', 'linestyle','--');
ylims = ylim;
text(1, ylims(1)+0.1,[num2str(size(allMeanPostTrend.wrong,2)) ','],   'color','r','fontsize',14);
text(20,ylims(1)+0.1,[num2str(size(allMeanPostTrend.correct,2)) ','], 'color','c','fontsize',14);
text(40,ylims(1)+0.1,[num2str(size(allMeanPostTrend.late,2)) ','],        'color','b','fontsize',14);
text(60,ylims(1)+0.1,[num2str(size(allMeanPostTrend.miss,2))],        'color','k','fontsize',14);
% legend('Early','Correct','Late','Miss')
hold off;
%% Individual experiment posteriors in reward zone

figure('Position',[ 496   451   944   420])
for exptIdx = 1:length(Posterior_all) %[1:3 5:8]
    es = meanPostTrend(exptIdx).es;
    numWrong   = size(meanPostTrend(exptIdx).all.wrong,1);
    numCorrect = size(meanPostTrend(exptIdx).all.correct,1);
    numLate    = size(meanPostTrend(exptIdx).all.late,1);
    numMiss    = size(meanPostTrend(exptIdx).all.miss,1);
    subplot(2,4,exptIdx)
    if numWrong>4
        errorarea_as(1:100, nanmean(meanPostTrend(exptIdx).all.wrong,1), nansem(meanPostTrend(exptIdx).all.wrong),'r', alpha_val)
        %         plot(1:100, (meanPostTrend(exptIdx).all.wrong),'r')
        hold on;
        LickPos.wrong = es.traj(es.outcome==2 & es.lick);
%         plot(es.traj(es.outcome==2 & es.lick), es.posterior(es.outcome==2 & es.lick)./es.norm_factor, 'r+');
%         plot(nanmean(LickPos.wrong)-nansem(LickPos.wrong), 1.5, 'r>');
%         plot(nanmean(LickPos.wrong)+nansem(LickPos.wrong), 1.5, 'r<');
    elseif numWrong>0
        plot(1:100, (meanPostTrend(exptIdx).all.wrong),'r')
        hold on; 
%         plot(es.traj(es.outcome==2 & es.lick), es.posterior(es.outcome==2 & es.lick)./es.norm_factor, 'r+');
    end;
    hold on;
    if numCorrect>4
        errorarea_as(1:100, nanmean(meanPostTrend(exptIdx).all.correct,1), nansem(meanPostTrend(exptIdx).all.correct),'c', alpha_val);
        hold on;
        LickPos.correct = es.traj(es.outcome==1 & es.lick);
        plot(nanmean(LickPos.correct)-nansem(LickPos.correct), 1.5, 'k>');
        plot(nanmean(LickPos.correct)+nansem(LickPos.correct), 1.5, 'k<');
    elseif numCorrect>0
        plot(1:100, (meanPostTrend(exptIdx).all.correct),'c')
        hold on;
%         plot(es.traj(es.outcome==1 & es.lick), es.posterior(es.outcome==1 & es.lick)./es.norm_factor, 'k+');
    end
    if numLate>4
        errorarea_as(1:100, nanmean(meanPostTrend(exptIdx).all.late,1), nansem(meanPostTrend(exptIdx).all.late),'b', alpha_val)
        hold on;
        LickPos.late = es.traj(es.outcome==0 & es.lick);
%         plot(es.traj(es.outcome==0 & es.lick), es.posterior(es.outcome==0 & es.lick)./es.norm_factor, 'b+');
%         plot(nanmean(LickPos.late)-nansem(LickPos.late), 1.5, 'b>');
%         plot(nanmean(LickPos.late)+nansem(LickPos.late), 1.5, 'b<');
    elseif numMiss>0
        plot(1:100, (meanPostTrend(exptIdx).all.miss),'k')
        hold on;
%         plot(es.traj(es.outcome==0 & es.lick), es.posterior(es.outcome==0 & es.lick)./es.norm_factor, 'b+');
    end;
    if numMiss>4
        errorarea_as(1:100, nanmean(meanPostTrend(exptIdx).all.miss,1), nansem(meanPostTrend(exptIdx).all.miss),'k', alpha_val)
        hold on;
        LickPos.miss = es.traj(es.outcome==0 & es.lick);
%         plot(es.traj(es.outcome==0 & es.lick), es.posterior(es.outcome==0 & es.lick)./es.norm_factor, 'b+');
%         plot(nanmean(LickPos.miss)-nansem(LickPos.miss), 1.5, 'b>');
%         plot(nanmean(LickPos.miss)+nansem(LickPos.miss), 1.5, 'b<');
    elseif numMiss>0
        plot(1:100, (meanPostTrend(exptIdx).all.miss),'k')
        hold on;
%         plot(es.traj(es.outcome==0 & es.lick), es.posterior(es.outcome==0 & es.lick)./es.norm_factor, 'b+');
    end;
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    set(gca, 'XLim',[-5 105])
    set(gca, 'YLim', [-2.5 2]);
    %     xlabel('Position (cm)')
    %     ylabel(['Posterior in ' num2str(2*start) ' to ' num2str(2*stop)]);
    line([70-5 70-5],ylim,'color','c');
    line([70+5 70+5],ylim,'color','c');
    ylims = ylim;
    text(1, ylims(1)+0.2,[num2str(numWrong) ','],   'color','r','fontsize',14);
    text(20,ylims(1)+0.2,[num2str(numCorrect) ','], 'color','c','fontsize',14);
    text(40,ylims(1)+0.2,[num2str(numLate) ','],        'color','k','fontsize',14);
    text(60,ylims(1)+0.2,[num2str(numMiss)],        'color','k','fontsize',14);
    hold off;
end
% figure;
% subplot(2,4,8)
% try
%     errorarea_as(1:100, nanmean(meanPostTrend(exptIdx).all.wrong,1), nansem(meanPostTrend(exptIdx).all.wrong),'r', alpha_val)
% catch; end;
% hold on;
% errorarea_as(1:100, nanmean(meanPostTrend(exptIdx).all.correct,1), nansem(meanPostTrend(exptIdx).all.correct),'c', alpha_val);
% try
%     errorarea_as(1:100, nanmean(meanPostTrend(exptIdx).all.miss,1), nansem(meanPostTrend(exptIdx).all.miss),'k', alpha_val)
% catch; end;
% axis([0 0.1 0 0.1])
% %     hold off
% plot(1,1)
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% legend('Early','Correct', 'Late', 'Location','EastOutside');
% 
% % subplot(2,4,1)
% % ylabel(['Normalized posterior in ' num2str(2*start) ' to ' num2str(2*stop)]);
% subplot(2,4,5)
% ylabel({'Normalized posterior';  ['in reward region']});% ' num2str(2*start) ' to ' num2str(2*stop) ' cm']});
% % for idx = 5:7
% xlabel('Position (cm)')
% % end
% subplot(2,4,8)
% axis off

%% Individual experiment, licks
% figure;
% for n = 1:8
%     subplot(2,4,n)
%     patch(2*[start stop stop start],2*[start start stop stop],'c');
%     hold on;
%     axis equal; axis square;
%     set(gca,'XLim',[0 100],'YLim',[0 100])
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     line(xlim, ylim, 'linestyle','--', 'color',[.5 .5 .5]);
%     line(xlim, 2*[start start ],'color','c');
%     line(xlim, 2*[stop stop ],'color','c');
%     line(2*[start start],ylim, 'color','c');
%     line(2*[stop stop],ylim, 'color','c');
%     
%     plot(meanPostTrend(n).licks.correct.traj  +randn(1)*0.1, meanPostTrend(n).licks.correct.X_ML +randn(1)*0.1,'.k');
%     xlabel('Actual position (cm)');
%     ylabel('MLE decoded position (cm)');
% end
figure('Position',[100 100 944 488]);
for n = 1:length(Posterior_all)
    subplot(2,4,n)
%     patch(2*[start stop stop start],2*[start start stop stop],'c');
    hold on;
    axis equal; axis square;
    set(gca,'XLim',[0 100],'YLim',[0 100])
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line(xlim, ylim, 'linestyle','--', 'color',[.5 .5 .5]);
    line(xlim, 2*[start start ],'color','c');
    line(xlim, 2*[stop stop ],'color','c');
    line(2*[start start],ylim, 'color','c');
    line(2*[stop stop],ylim, 'color','c');
    
    plot(meanPostTrend(n).licks.wrong.traj +randn(1)*0.1, meanPostTrend(n).licks.wrong.X_ML+randn(1)*0.1,'r+');
    plot(meanPostTrend(n).licks.late.traj  +randn(1)*0.1, meanPostTrend(n).licks.late.X_ML +randn(1)*0.1,'k+');
    xlabel('Actual position (cm)');
    ylabel('MLE decoded position (cm)');
% end
% 
% %%
% figure;
% for n = 1:8
    subplot(2,4,n)
%     patch(2*[start stop stop start],2*[start start stop stop],'c');
    hold on;
    axis equal; axis square;
    set(gca,'XLim',[0 100],'YLim',[0 100])
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line(xlim, ylim, 'linestyle','--', 'color',[.5 .5 .5]);
    line(xlim, 2*[start start ],'color','c');
    line(xlim, 2*[stop stop ],'color','c');
    line(2*[start start],ylim, 'color','c');
    line(2*[stop stop],ylim, 'color','c');
    
    a=nanstd(meanPostTrend(n).licks.correct.traj); %horizontal radius
    b=nanstd(meanPostTrend(n).licks.correct.X_ML); %vertical radius
    x0=nanmean(meanPostTrend(n).licks.correct.traj); % x0,y0 ellipse centre coordinates
    y0=nanmean(meanPostTrend(n).licks.correct.X_ML);
    t=-pi:0.01:pi;
    x=x0+a*cos(t);
    y=y0+b*sin(t);
    plot(x,y,'color',[.5 .5 .5],'linewidth',1.5)
    
    a=nanstd(meanPostTrend(n).licks.wrong.traj); %horizontal radius
    b=nanstd(meanPostTrend(n).licks.wrong.X_ML); %vertical radius
    x0=nanmean(meanPostTrend(n).licks.wrong.traj); % x0,y0 ellipse centre coordinates
    y0=nanmean(meanPostTrend(n).licks.wrong.X_ML);
    t=-pi:0.01:pi;
    x=x0+a*cos(t);
    y=y0+b*sin(t);
%     plot(x,y,'r','linewidth',1.5)
    patch(x,y,'r','linewidth',1.5,'FaceAlpha',alpha_val)
    
    a=nanstd(meanPostTrend(n).licks.late.traj); %horizontal radius
    b=nanstd(meanPostTrend(n).licks.late.X_ML); %vertical radius
    x0=nanmean(meanPostTrend(n).licks.late.traj); % x0,y0 ellipse centre coordinates
    y0=nanmean(meanPostTrend(n).licks.late.X_ML);
    t=-pi:0.01:pi;
    x=x0+a*cos(t);
    y=y0+b*sin(t);
%     plot(x,y,'color',[0 .5 0],'linewidth',1.5)
    patch(x,y,'k','linewidth',1.5,'FaceAlpha',alpha_val)
    
    xlabel('Actual position (cm)');
    ylabel('MLE decoded position (cm)');
end

%% All experiments, licks
figure('Position',[9         400        1194         717]);

subplot(121)
traj_vals_c = [];
ML_vals_c   = [];
traj_vals_w = [];
ML_vals_w   = [];
traj_vals_m = [];
ML_vals_m   = [];
for n = 1:length(Posterior_all) %[1:3 5:8]
       
%     a=nanstd(meanPostTrend(n).licks.correct.traj); %horizontal radius
%     b=nanstd(meanPostTrend(n).licks.correct.X_ML); %vertical radius
%     x0=nanmean(meanPostTrend(n).licks.correct.traj); % x0,y0 ellipse centre coordinates
%     y0=nanmean(meanPostTrend(n).licks.correct.X_ML);
%     t=-pi:0.01:pi;
%     x=x0+a*cos(t);
%     y=y0+b*sin(t);
%     plot(x,y,'color',[.5 .5 .5],'linewidth',1.5)
    traj_vals_w = [traj_vals_w meanPostTrend(n).licks.wrong.traj'];
    ML_vals_w   = [ML_vals_w meanPostTrend(n).licks.wrong.X_ML'];
    traj_vals_m = [traj_vals_m meanPostTrend(n).licks.late.traj'];
    ML_vals_m   = [ML_vals_m meanPostTrend(n).licks.late.X_ML'];
    traj_vals_c = [traj_vals_c meanPostTrend(n).licks.correct.traj'];
    ML_vals_c   = [ML_vals_c meanPostTrend(n).licks.correct.X_ML'];
    
    a=nanstd(meanPostTrend(n).licks.wrong.traj); %horizontal radius
    b=nanstd(meanPostTrend(n).licks.wrong.X_ML*50/numBins); %vertical radius
    x0=nanmean(meanPostTrend(n).licks.wrong.traj); % x0,y0 ellipse centre coordinates
    y0=nanmean(meanPostTrend(n).licks.wrong.X_ML*50/numBins);
    t=-pi:0.01:pi;
    x=x0+a*cos(t);
    y=y0+b*sin(t);
    patch(x,y,'r','linewidth',1.5,'FaceAlpha',alpha_val)
    hold on;
    
    a=nanstd(meanPostTrend(n).licks.late.traj); %horizontal radius
    b=nanstd(meanPostTrend(n).licks.late.X_ML*50/numBins); %vertical radius
    x0=nanmean(meanPostTrend(n).licks.late.traj); % x0,y0 ellipse centre coordinates
    y0=nanmean(meanPostTrend(n).licks.late.X_ML*50/numBins);
    t=-pi:0.01:pi;
    x=x0+a*cos(t);
    y=y0+b*sin(t);
    patch(x,y,'k','linewidth',1.5,'FaceAlpha',alpha_val)
    
    if n == 4
        legend('Early','Late');
    end
    
    patch(2*[start stop stop start],2*[start start stop stop],'c');
    hold on;
    axis equal; axis square;
    set(gca,'XLim',[0 100],'YLim',[0 100])
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line(xlim, ylim, 'linestyle','--', 'color',[.5 .5 .5]);
    line(xlim, 2*[start start ],'color','c');
    line(xlim, 2*[stop stop ],'color','c');
    line(2*[start start],ylim, 'color','c');
    line(2*[stop stop],ylim, 'color','c');
    
    xlabel('Actual position (cm)');
    ylabel('Decoded position, MAP (cm)');
    title('Lick distributions');
    
end

ML_vals_c = ML_vals_c*50/numBins;
ML_vals_w = ML_vals_w*50/numBins;
ML_vals_m = ML_vals_m*50/numBins;

X = 1:6:100;
% for exptIdx = 1:8
%     trajDistr_c(exptIdx,:) = hist(meanPostTrend(n).licks.correct.traj, X)./sum(~isnan(meanPostTrend(n).licks.correct.traj));
%     trajDistr_w(exptIdx,:) = hist(meanPostTrend(n).licks.wrong.traj, X)./sum(~isnan(meanPostTrend(n).licks.wrong.traj));
%     trajDistr_m(exptIdx,:) = hist(meanPostTrend(n).licks.late.traj, X)./sum(~isnan(meanPostTrend(n).licks.late.traj));
%     
%     ML_Distr_c(exptIdx,:) = hist(meanPostTrend(n).licks.correct.X_ML, X)./sum(~isnan(meanPostTrend(n).licks.correct.X_ML));
%     ML_Distr_w(exptIdx,:) = hist(meanPostTrend(n).licks.wrong.X_ML, X)./sum(~isnan(meanPostTrend(n).licks.wrong.X_ML));
%     ML_Distr_m(exptIdx,:) = hist(meanPostTrend(n).licks.late.X_ML, X)./sum(~isnan(meanPostTrend(n).licks.late.X_ML));
% end
% This is the older version
trajDistr_c = hist(traj_vals_c, X)./length(traj_vals_c);
ML_Distr_c = hist(ML_vals_c, X)./sum(~isnan(ML_vals_c));
trajDistr_w = hist(traj_vals_w, X)./length(traj_vals_w);
ML_Distr_w = hist(ML_vals_w, X)./sum(~isnan(ML_vals_w));
trajDistr_m = hist(traj_vals_m, X)./length(traj_vals_m);
ML_Distr_m = hist(ML_vals_m, X)./sum(~isnan(ML_vals_m));


% ML_Distr = ML_Distr./sum(ML_Distr);
% trajDistr = trajDistr./sum(trajDistr);
subplot(222)
plot(X, nanmean(trajDistr_c,1), 'c',X, nanmean(trajDistr_w,1), 'r', X, nanmean(trajDistr_m,1),'k','linewidth',2);
line([2*start 2*start], ylim, 'color','c');
line([2*stop 2*stop], ylim, 'color','c');
axis square;
set(gca,'XLim',[0 100])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
legend('Correct','Licked early', 'Licked late','Location','Best');
xlabel('Actual Position (cm)');
ylabel('Density of licks');
subplot(224)
plot(X, nanmean(ML_Distr_c,1), 'c',X, nanmean(ML_Distr_w,1), 'r', X, nanmean(ML_Distr_m,1),'k','linewidth',2);
line([2*start 2*start], ylim, 'color','c');
line([2*stop 2*stop], ylim, 'color','c');
axis square;
set(gca,'XLim',[0 100])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
legend('Correct','Licked early', 'Licked late','Location','Best');
xlabel('Decoded Position (cm)');
ylabel('Density of licks');

% %%
% figure;
% % traj_vals_c
% 
% a=nanstd(traj_vals_c); %horizontal radius
% b=nanstd(ML_vals_c); %vertical radius
% x0=nanmean(traj_vals_c); % x0,y0 ellipse centre coordinates
% y0=nanmean(ML_vals_c);
% t=-pi:0.01:pi;
% x=x0+a*cos(t);
% y=y0+b*sin(t);
% patch(x,y,'c','linewidth',1.5,'FaceAlpha',alpha_val)
% hold on;
% a=nanstd(traj_vals_w); %horizontal radius
% b=nanstd(ML_vals_w); %vertical radius
% x0=nanmean(traj_vals_w); % x0,y0 ellipse centre coordinates
% y0=nanmean(ML_vals_w);
% t=-pi:0.01:pi;
% x=x0+a*cos(t);
% y=y0+b*sin(t);
% patch(x,y,'r','linewidth',1.5,'FaceAlpha',alpha_val)
% hold on;
% a=nanstd(traj_vals_m); %horizontal radius
% b=nanstd(ML_vals_m); %vertical radius
% x0=nanmean(traj_vals_m); % x0,y0 ellipse centre coordinates
% y0=nanmean(ML_vals_m);
% t=-pi:0.01:pi;
% x=x0+a*cos(t);
% y=y0+b*sin(t);
% patch(x,y,'k','linewidth',1.5,'FaceAlpha',alpha_val)
% hold on;
% axis square;
% set(gca,'XLim',[0 100],'YLim',[0 100])
% line(xlim, ylim, 'linestyle','--','color','k');
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% legend('Correct','Licked early', 'Licked late','Location','Best');
% line([2*start 2*start], ylim, 'color','c');
% line([2*stop 2*stop], ylim, 'color','c');
% line(xlim, [2*start 2*start], 'color','c');
% line(xlim, [2*stop 2*stop], 'color','c');
% xlabel('Decoded Position (cm)');
% ylabel('Density of licks');

%% To get the licks from reward position, marginals
X2 = -70:10:30;

figure;

trajDistr_fromRew_c = hist(traj_vals_c-70, X2)./length(traj_vals_c);
ML_Distr_fromRew_c = hist(ML_vals_c-70, X2)./sum(~isnan(ML_vals_c));

trajDistr_fromRew_e = hist([traj_vals_w traj_vals_m]-70, X2)./length([traj_vals_w traj_vals_m]);
ML_Distr_fromRew_e = hist([ML_vals_w ML_vals_m]-70, X2)./sum(~isnan([ML_vals_w ML_vals_m]));

trajDistr_fromRew_w = hist([traj_vals_w]-70, X2)./length([traj_vals_w]);
ML_Distr_fromRew_w = hist([ML_vals_w]-70, X2)./sum(~isnan([ML_vals_w]));

trajDistr_fromRew_m = hist([traj_vals_m]-70, X2)./length([traj_vals_m]);
ML_Distr_fromRew_m = hist([ML_vals_m]-70, X2)./sum(~isnan([ML_vals_m]));

% % Newer version
% trajDistr_fromRew_c = nanmean(trajDistr_c,1);
% ML_Distr_fromRew_c  = nanmean(ML_Distr_c,1);
% 
% trajDistr_fromRew_e = nanmean(trajDistr_m,1);
% ML_Distr_fromRew_e  = nanmean(ML_Distr_m,1);
subplot(221)
hold off
% plot(X2, trajDistr_fromRew_c, 'k--','linewidth',1.5)
% hold on; plot(X2, ML_Distr_fromRew_c, 'k-','linewidth',1.5);
patch([X2 X2(end) X2(1)], [trajDistr_fromRew_c 0 0], [0.5 0.5 0.5])
hold on; 
patch([X2 X2(end) X2(1)], [ML_Distr_fromRew_c 0 0], [0 0 0])
axis tight;
ylims = ylim;
patch([-5 5 5 -5], [ylims(1) ylims(1) ylims(2) ylims(2)], [0 1 1], 'EdgeColor','none');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% legend('Position','Decoded position','Location','Best');
xlabel('Distance from reward location (cm)');
ylabel('Density of licks');
set(gca, 'XLim', [-50 30]);
% axis square
title('Correct trials')

subplot(222)
% plot(X2, trajDistr_fromRew_w, 'r--','linewidth',1.5)
patch([X2 X2(end) X2(1)], [trajDistr_fromRew_w 0 0], [0.5 0 0])
hold on; 
patch([X2 X2(end) X2(1)], [ML_Distr_fromRew_w 0 0], [1 0 0])
% plot(X2, ML_Distr_fromRew_w, 'r-','linewidth',1.5);
% plot(X2, trajDistr_fromRew_e, 'k--','linewidth',1.5)
% hold on; plot(X2, ML_Distr_fromRew_e, 'k-','linewidth',1.5);
axis tight;
ylims = ylim;
patch([-5 5 5 -5], [ylims(1) ylims(1) ylims(2) ylims(2)], [0 1 1], 'EdgeColor','none');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
set(gca, 'XLim', [-50 30])
% legend('Position','Decoded position','Location','Best');
xlabel('Distance from reward location (cm)');
ylabel('Density of licks');
% axis square
title('Early trials')
% title('Error trials')

subplot(224)
% plot(X2, trajDistr_fromRew_m, 'b--','linewidth',1.5)
% hold on; plot(X2, ML_Distr_fromRew_m, 'b-','linewidth',1.5);
patch([X2 X2(end) X2(1)], [trajDistr_fromRew_m 0 0], [0 0 0.5])
hold on; 
patch([X2 X2(end) X2(1)], [ML_Distr_fromRew_m 0 0], [0 0 1])
axis tight;
ylims = ylim;
patch([-5 5 5 -5], [ylims(1) ylims(1) ylims(2) ylims(2)], [0 1 1], 'EdgeColor','none');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
set(gca, 'XLim', [-50 30])
% legend('Position','Decoded position','Location','Best');
xlabel('Distance from reward location (cm)');
ylabel('Density of licks');
% axis square
title('Late trials')