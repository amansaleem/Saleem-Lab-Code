


% function plot_summary_figures
for idx = 1:length(Posterior_all)
    Posterior_all(idx).confidence.high  = [nanmax(2.^Posterior_all(idx).Posterior_high,[],2) - nanmin(2.^Posterior_all(idx).Posterior_high,[],2)]'/50;
    Posterior_all(idx).confidence.norm  = [nanmax(2.^Posterior_all(idx).Posterior_norm,[],2) - nanmin(2.^Posterior_all(idx).Posterior_norm,[],2)]'/50;
    Posterior_all(idx).confidence.low   = [nanmax(2.^Posterior_all(idx).Posterior_low,[],2)  - nanmin(2.^Posterior_all(idx).Posterior_low,[],2)]'/50;
    Posterior_all(idx).confidence.gray  = [nanmax(2.^Posterior_all(idx).Posterior_gray,[],2) - nanmin(2.^Posterior_all(idx).Posterior_gray,[],2)]'/50;
end
%% all marginals (actual)
figure

for idx = 1:length(Posterior_all)
    subplot(2,4,idx)
    plot(Posterior_all(idx).marginals.lowX*2*sqrt(2),Posterior_all(idx).marginals.low)
    hold on;
    plot(Posterior_all(idx).marginals.normX*2*sqrt(2),Posterior_all(idx).marginals.norm,'k')
    plot(Posterior_all(idx).marginals.highX*2*sqrt(2),Posterior_all(idx).marginals.high,'r');
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    xlabel('Distance from original position (cm)');
    ylabel('Decoded posterior');
    line(xlim, [0 0], 'linestyle','--','color','k');
    line([0 0], ylim, 'linestyle','--','color','k');
end

%% All Marginals (normalised)

% figure;
% for idx = 1:length(Posterior_all)
%     subplot(2,4,idx)
%     nfactor = max(Posterior_all(idx).marginals.norm);
%     plot(Posterior_all(idx).marginals.lowX*2*sqrt(2),Posterior_all(idx).marginals.low./nfactor)
%     hold on;
%     plot(Posterior_all(idx).marginals.normX*2*sqrt(2),Posterior_all(idx).marginals.norm./nfactor,'k')
%     plot(Posterior_all(idx).marginals.highX*2*sqrt(2),Posterior_all(idx).marginals.high./nfactor,'r');
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     xlabel('Distance from original position (cm)');
%     ylabel('Decoded posterior');
%     line(xlim, [0 0], 'linestyle','--','color','k');
%     line([0 0], ylim, 'linestyle','--','color','k');
% end

%% Summary marginals
for idx = 1:length(Posterior_all)
    nfactor = max(Posterior_all(idx).marginals.norm);
    allhigh(idx,:) = Posterior_all(idx).marginals.high./nfactor;
    allnorm(idx,:) = Posterior_all(idx).marginals.norm./nfactor;
    alllow(idx,:) = Posterior_all(idx).marginals.low./nfactor;
end

figure
subplot(121)
plot(Posterior_all(1).marginals.lowX*2*sqrt(2), allhigh,'r')
hold on; plot(Posterior_all(1).marginals.lowX*2*sqrt(2), allnorm,'k')
hold on; plot(Posterior_all(1).marginals.lowX*2*sqrt(2), alllow,'b')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlabel('Distance from original position (cm)');
ylabel('Decoded posterior');
line(xlim, [0 0], 'linestyle','--','color','k');
line([0 0], ylim, 'linestyle','--','color','k');
set(gca, 'YLim', [-2 1.5])

subplot(122)
errorarea_as(Posterior_all(1).marginals.lowX*2*sqrt(2), mean(allhigh,1),sem(allhigh),'r')
hold on;
errorarea_as(Posterior_all(1).marginals.lowX*2*sqrt(2), mean(allnorm),sem(allnorm),'k')
errorarea_as(Posterior_all(1).marginals.lowX*2*sqrt(2), mean(alllow), sem(alllow),'b')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlabel('Distance from original position (cm)');
ylabel('Decoded posterior');
line(xlim, [0 0], 'linestyle','--','color','k');
line([0 0], ylim, 'linestyle','--','color','k');
set(gca, 'YLim', [-2 1.5])


% %% to plot all mean firing rates;
% figure
%
% for idx = 1:8
%     subplot(2,4,idx)
%     plot(Posterior_all(idx).meanrates.low, Posterior_all(idx).meanrates.high,'ko');
%     axis equal; axis square;
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     maxLim = max(max(xlim), max(ylim));
%     set(gca, 'xlim', [0 maxLim], 'ylim', [0 maxLim]);
%     line(xlim, ylim, 'color','k', 'linestyle','--');
% end

%%
colour_list = [...
    1 0 0;
    0 1 0;
    0 0 1;
    0 1 1;
    0 0 0;
    1 1 0;
    1 0 1;
    0.5 1 0;
    0 0.5 1;
    0.5 0.5 0.5;
    1 0.5 0;
    1 0 0.5;
    0 1 0.5...
    ];

symbol_list = [...
    '.','o','s','d',...
    '.','o','s','d',...
    '.','o','s','d',...
    '.','o','s','d',...
    ];

%% all mean firing rates
figure;
thres = 0.00;
diffRate_all = [];
for idx = 1:length(Posterior_all)
    subplot(231)
    
    subsetCells = nanmean(Posterior_all(idx).decoder_high.model.EV)>thres &...
        nanmean(Posterior_all(idx).decoder_low.model.EV)>thres; 
    
    plot(Posterior_all(idx).meanrates.low, Posterior_all(idx).meanrates.high,...
        'color',colour_list(idx,:),'Marker', symbol_list(idx),'linestyle','none');
    hold on;
    axis equal; axis square;
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    maxLim = max(max(xlim), max(ylim));
    maxLim = 2.5;
    set(gca, 'xlim', [0 maxLim], 'ylim', [0 maxLim]);
    line(xlim, ylim, 'color','k', 'linestyle','--');
    ylabel('Mean Rate @ high contrast');
    xlabel('Mean Rate @ low contrast');
    subplot(234)
    errorbarxy(...
        nanmean(60*Posterior_all(idx).meanrates.low(subsetCells)),...
        nanmean(60*Posterior_all(idx).meanrates.high(subsetCells)),...
        nansem (60*Posterior_all(idx).meanrates.low(subsetCells)),...
        nansem (60*Posterior_all(idx).meanrates.high(subsetCells)),{'ko','k','k'})%...
    p_rates(idx) = signrank(Posterior_all(idx).meanrates.high(subsetCells),Posterior_all(idx).meanrates.low(subsetCells));
    diffRate(idx) = nanmean(Posterior_all(idx).meanrates.high(subsetCells)) - nanmean(Posterior_all(idx).meanrates.low(subsetCells));
    diffRate_all = [diffRate_all Posterior_all(idx).meanrates.high(subsetCells) - Posterior_all(idx).meanrates.low(subsetCells)];
    
    hold on;
    axis equal; axis tight;
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    maxLim = max(max(xlim), max(ylim));
%     maxLim = 1.5;
    axis square;
    set(gca, 'xlim', [0 maxLim], 'ylim', [0 maxLim]);
    line(xlim, ylim, 'color','k', 'linestyle','--');
    ylabel('Mean Rate @ high contrast');
    xlabel('Mean Rate @ low contrast');
%     axis([0 1 0 1]);
end

%% differences in explained variance
subplot(233)
thres = 0.01;

for idx = 1:length(Posterior_all)
    lowEV = nanmean(Posterior_all(idx).decoder_low.model.EV,1);
    highEV = nanmean(Posterior_all(idx).decoder_high.model.EV,1);
    normEV = nanmean(Posterior_all(idx).decoder.model.EV,1);
    t = (lowEV>thres & highEV>thres) & normEV>thres;
    try
        lowEV_decP = nanmean(Posterior_all(idx).decoder_low_decPos.model.EV,1);
        highEV_decP = nanmean(Posterior_all(idx).decoder_high_decPos.model.EV,1);
        normEV_decP = nanmean(Posterior_all(idx).decoder_norm_decPos.model.EV,1);
        t_decP = (lowEV_decP>thres & highEV_decP>thres) & normEV_decP>thres;
    catch; end
    diffEVs(idx) = nanmean(highEV(t)./normEV(t)-lowEV(t)./normEV(t));
    
    diffEV(idx) = nanmean(highEV-lowEV);
    try
    ps(idx) = signrank(lowEV(t)-highEV(t));
    p(idx) = signrank(lowEV-highEV);
    catch; end
    try
        diffEVs_decP(idx) = nanmean(highEV_decP(t_decP)./normEV_decP(t_decP)-lowEV_decP(t_decP)./normEV_decP(t_decP));
        diffEV_decP(idx) = nanmean(highEV_decP(t_decP)-lowEV_decP(t_decP));
        p_decP(idx) = signrank(lowEV_decP-highEV_decP);
        ps_decP(idx) = signrank(lowEV_decP(t_decP)-highEV_decP(t_decP));
    catch; end
    %     plot(lowEV(t), highEV(t),...
    %         'color',colour_list(idx,:),'Marker', symbol_list(idx),'linestyle','none');
    %     plot(nanmean(lowEV(t)./normEV(t)), nanmean(highEV(t)./normEV(t)),'ko')%...
%     plot(nanmean(lowEV(t)), nanmean(highEV(t)),'ko')%...
        errorbarxy(...
        nanmean(lowEV(t)./normEV(t)),...
        nanmean(highEV(t)./normEV(t)),...
        nansem(lowEV(t)./normEV(t)),...
        nansem(highEV(t)./normEV(t)),...
        {'ko','k','k'})%...
%             'color',colour_list(idx,:),'Marker', symbol_list(idx),'linestyle','none');
    hold on;
    if idx ==4
        errorbarxy(...
        nanmean(lowEV(t)./normEV(t)),...
        nanmean(highEV(t)./normEV(t)),...
        nansem(lowEV(t)./normEV(t)),...
        nansem(highEV(t)./normEV(t)),...
        {'ro','k','k'})%...
    end
    hold on;
    try
%         plot(nanmean(lowEV_decP(t_decP)), nanmean(highEV_decP(t_decP)),'ro')%...
%         errora
        errorbarxy(...
            nanmean(lowEV_decP(t_decP)./normEV_decP(t_decP)),...
            nanmean(highEV_decP(t_decP)./normEV_decP(t_decP)),...
            nansem(lowEV_decP(t_decP)./normEV_decP(t_decP)),...
            nansem(highEV_decP(t_decP)./normEV_decP(t_decP)),...
            {'ro','r','r'})%...
        hold on;
        if idx ==4
        errorbarxy(...
            nanmean(lowEV_decP(t_decP)./normEV_decP(t_decP)),...
            nanmean(highEV_decP(t_decP)./normEV_decP(t_decP)),...
            nansem(lowEV_decP(t_decP)./normEV_decP(t_decP)),...
            nansem(highEV_decP(t_decP)./normEV_decP(t_decP)),...
            {'ko','r','r'})%...
        hold on;
    end
    catch; end;
    axis equal; axis square;
    axis tight
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    maxLim = max(max(xlim), max(ylim));
    minLim = min(min(xlim), min(ylim));
    set(gca, 'xlim', [minLim maxLim], 'ylim', [minLim maxLim]);
end
set(gca, 'xlim', [0 maxLim+0.01], 'ylim', [0 maxLim+0.01]);
line(xlim, ylim, 'color','k', 'linestyle','--');
xlabel('EV @ low contrast');
ylabel('EV @ high contrast');

%% EV actual vs. estimated
subplot(236)
hold off
thres = 0.01;

try
    diffdiffEV = [];
    diffdiffEVlow = [];
    diffEVall = [];
    diffEVall_decP = [];
    
for idx = 1:length(Posterior_all)
    normEV = nanmean(Posterior_all(idx).decoder_norm_decPos.model.EV_orig,1);
    lowEV = nanmean(Posterior_all(idx).decoder_low_decPos.model.EV_orig,1);
    highEV = nanmean(Posterior_all(idx).decoder_high_decPos.model.EV_orig,1);
    
    normEV_decP = nanmean(Posterior_all(idx).decoder_norm_decPos.model.EV,1);
    lowEV_decP = nanmean(Posterior_all(idx).decoder_low_decPos.model.EV,1);
    highEV_decP = nanmean(Posterior_all(idx).decoder_high_decPos.model.EV,1);
    
    t_decP = (lowEV_decP>thres & highEV_decP>thres);% & normEV_decP>thres;
    t      = (lowEV>thres & highEV>thres);% & normEV>thres;
    t_all  = t & t_decP;
    
%     lowEV = lowEV./normEV;
%     lowEV_decP = lowEV_decP./normEV_decP;
%     highEV = highEV./normEV;
%     highEV_decP = highEV_decP./normEV_decP;
       
    lowEV(isinf(lowEV)) = nan;
    lowEV_decP(isinf(lowEV_decP)) = nan;
    highEV(isinf(highEV)) = nan;
    highEV_decP(isinf(highEV_decP)) = nan;
    
    numCells(idx) = sum(t_all);
    diff_dec_act_high(idx) = nanmean(highEV_decP(t_all) - highEV(t_all)); 
    diff_dec_act_low(idx)  = nanmean(lowEV_decP(t_all) - lowEV(t_all)); 
%     diffdiffEV = [diffdiffEV ((highEV_decP(t_all)-lowEV_decP(t_all)) - (highEV(t_all)-lowEV(t_all)))];
    diffEVs(idx) = nanmean((highEV(t)./normEV(t) - lowEV(t))./normEV(t));
    ps(idx) = signrank(lowEV(t)-highEV(t));
    
    diffEVall = [diffEVall ((highEV(t_all)-lowEV(t_all)))];
    diffEV(idx) = nanmean(highEV-lowEV);
    p(idx) = signrank(lowEV-highEV);
    
    diffEVs_decP(idx) = nanmean((highEV_decP(t_decP)./normEV_decP(t_decP)-lowEV_decP(t_decP)./normEV_decP(t_decP)));
    diffEV_decP(idx) = nanmean((highEV_decP-lowEV_decP));
    diffEVall_decP = [diffEVall_decP (highEV_decP(t_decP)-lowEV_decP(t_decP))];
    
    p_decP(idx) = signrank(lowEV_decP-highEV_decP);
    ps_decP(idx) = signrank(lowEV_decP(t_decP)-highEV_decP(t_decP));
    %     plot(lowEV(t), highEV(t),...
    %         'color',colour_list(idx,:),'Marker', symbol_list(idx),'linestyle','none');
    %     plot(nanmean(lowEV(t)./normEV(t)), nanmean(highEV(t)./normEV(t)),'ko')%...
    %     plot(nanmean(lowEV(t)), nanmean(highEV(t)),'ko')%...
    %     t_all = t & t_decP;
    
%     diffEVs_decP(idx) = nanmean((highEV_decP(t_all)-lowEV_decP(t_all)));
%     t_all = lowEV>thres | lowEV_decP>thres;
    errorbarxy(nanmean(lowEV(t_all)), nanmean(lowEV_decP(t_all)), nansem(lowEV(t_all)), nansem(lowEV_decP(t_all)),{'bo','b','b'})%...
    %         'color',colour_list(idx,:),'Marker', symbol_list(idx),'linestyle','none');
    hold on;
    diffdiffEVlow = [diffdiffEVlow (lowEV_decP(t_all)-lowEV(t_all))];
    diffEVlow(idx) = nanmean((lowEV_decP(t_all)-lowEV(t_all)));
    p_diffEVlow(idx) = signrank((lowEV_decP(t_all)-lowEV(t_all)));
    
%         plot(nanmean(lowEV_decP(t_decP)), nanmean(highEV_decP(t_decP)),'ro')%...
%     t_all = highEV>thres | highEV_decP>thres;
    diffEVhigh(idx) = nanmean((highEV_decP(t_all)-highEV(t_all)));
%     errorbarxy(nanmean(highEV(t_all)), nanmean(highEV_decP(t_all)),nansem(highEV(t_all)), nansem(highEV_decP(t_all)),{'ro','r','r'})%...
    hold on;
%     t_all = normEV>thres | normEV_decP>thres;
    diffdiffEV = [diffdiffEV (normEV_decP(t_all)-normEV(t_all))];
    diffEVnorm(idx) = nanmean((normEV_decP(t_all)-normEV(t_all)));
%         errorbarxy(nanmean(normEV(t_all)), nanmean(normEV_decP(t_all)),nansem(normEV(t_all)), nansem(normEV_decP(t_all)),{'ko','k','k'})%...
    hold on;
    axis equal; axis square;
    axis tight
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    maxLim = max(max(xlim), max(ylim));
    minLim = min(min(xlim), min(ylim));
    set(gca, 'xlim', [minLim maxLim], 'ylim', [minLim maxLim]);
end
set(gca, 'xlim', [0 maxLim+0.01], 'ylim', [0 maxLim+0.01]);
line(xlim, ylim, 'color','k', 'linestyle','--');
xlabel('EV, actual position');
ylabel('EV, estimated contrast');
catch; end
%% differences in widths
subplot(235)
hold off
subplot(232)
hold off
for idx = 1:length(Posterior_all)
    lowEV = nanmean(Posterior_all(idx).decoder_low.model.EV,1);
    highEV = nanmean(Posterior_all(idx).decoder_high.model.EV,1);
    normEV = nanmean(Posterior_all(idx).decoder.model.EV,1);
    try
        lowEV_decP = nanmean(Posterior_all(idx).decoder_low_decPos.model.EV,1);
        highEV_decP = nanmean(Posterior_all(idx).decoder_high_decPos.model.EV,1);
        normEV_decP = nanmean(Posterior_all(idx).decoder_norm_decPos.model.EV,1);
        lowW_decP = Posterior_all(idx).fieldWidths.low_decPos;
        normW_decP = Posterior_all(idx).fieldWidths.norm_decPos;
        highW_decP = Posterior_all(idx).fieldWidths.high_decPos;
    catch; end
    
    lowW = Posterior_all(idx).fieldWidths.low;
    highW = Posterior_all(idx).fieldWidths.high;
    normW = Posterior_all(idx).fieldWidths.norm;
    
    thres = 0.01;
    cthres_high= 51;
    cthres_low = 0;
    [mv,c] = max(Posterior_all(idx).decoder.model.bestModel');
%     [mvl,c] = max(Posterior_all(idx).decoder_low.model.bestModel');
%     [mvh,c] = max(Posterior_all(idx).decoder_high.model.bestModel');
%     mv(idx) = nanmean(mvh-mvl);
    
    t = (lowEV>=thres) & (highEV>=thres) ...
        & c<cthres_high & c>cthres_low ...
        & (highW<50 & lowW<50);
    try
        t_decP = ((lowEV_decP>=thres) & (highEV_decP>=thres));%  & c<cthres_high & c>cthres_low;% & (normEV>=thres);
        t_all = t | t_decP;
        diffW_decP(idx) = nanmean(highW_decP(t_all)-lowW_decP(t_all));
        pW_decP(idx) = signrank(lowW_decP(t_all)-highW_decP(t_all));
    catch
        t_all = t;
    end
    %     display(num2str([idx sum(t_all) pW(idx)]))
    
%     diffEV(idx) = nanmean(highEV(t_all)-lowEV(t_all));
    p(idx) = signrank(lowEV(t_all)-highEV(t_all));
%     t_all = t;
    diffW(idx) = nanmean(highW(t)-lowW(t));
%     pW(idx) = signrank(lowW(t)-highW(t));
    
    subplot(232)
%     if pW(idx)<0.05
%         scatter(nanmean(lowW(t_all)), nanmean(highW(t_all)),'ko','filled');
%     else
%         scatter(nanmean(lowW(t_all)), nanmean(highW(t_all)),'ko');%...
%     end
    errorbarxy(nanmean(lowW(t_all)), nanmean(highW(t_all)),nansem(lowW(t_all)), nansem(highW(t_all)),{'ko','k','k'})
    %         'color',colour_list(idx,:),'Marker', symbol_list(idx),'linestyle','none');
    %     plot((lowW), (highW),'ko')%...
    %     plot((lowW_decP(t_all)), (highW_decP(t_all)),'ro')%...
    hold on;
    try
%         if pW_decP(idx)<0.05
%             scatter(nanmean(lowW_decP(t_all)), nanmean(highW_decP(t_all)),'ro','filled')%...
%         else
%             scatter(nanmean(lowW_decP(t_all)), nanmean(highW_decP(t_all)),'ro')%...
%         end
        errorbarxy(nanmean(lowW_decP(t_all)), nanmean(highW_decP(t_all)),nansem(lowW_decP(t_all)), nansem(highW_decP(t_all)),{'ro','r','r'})%...
        hold on;
    catch; end
    %     plot(nanmean(lowW_decP), nanmean(highW_decP),'ro')%...
    axis equal; axis square;
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    axis tight
    maxLim = max(max(xlim), max(ylim));
    minLim = min(min(xlim), min(ylim));
    set(gca, 'xlim', [minLim maxLim], 'ylim', [minLim maxLim]);
    %     pause
    %     hold off
    try
        subplot(235)
        pW_decP_P_low(idx) = signrank(lowW_decP(t_all)-lowW(t_all));
        pW_decP_P_high(idx) = signrank(highW_decP(t_all)-highW(t_all));
%         if pW_decP_P_low(idx)<0.05
%             scatter(nanmean(lowW(t_all)), nanmean(lowW_decP(t_all)),[],[0 0.25 0],'o','filled');
%         else
%             scatter(nanmean(lowW(t_all)), nanmean(lowW_decP(t_all)),[],[0 0.25 0],'o');%...
%         end
        errorbarxy(nanmean(lowW(t_all)), nanmean(lowW_decP(t_all)),nansem(lowW(t_all)), nansem(lowW_decP(t_all)),{'ko','k','k'}); %...
        %         'color',colour_list(idx,:),'Marker', symbol_list(idx),'linestyle','none');
    %     plot((lowW), (highW),'ko')%...
    %     plot((lowW_decP(t_all)), (highW_decP(t_all)),'ro')%...
        hold on;
%         if pW_decP_P_high(idx)<0.05
%             scatter(nanmean(highW(t_all)), nanmean(highW_decP(t_all)),[],[0.5 0 0.5],'o','filled')%...
%         else
%             scatter(nanmean(highW(t_all)), nanmean(highW_decP(t_all)),[],[0.5 0 0.5],'o')%...
%         end
        errorbarxy(nanmean(highW(t_all)), nanmean(highW_decP(t_all)),nansem(highW(t_all)), nansem(highW_decP(t_all)),{'ro','r','r'}); %...
        hold on;
        errorbarxy(nanmean(normW(t_all)), nanmean(normW_decP(t_all)),nansem(normW(t_all)), nansem(normW_decP(t_all)),{'bo','b','b'}); %...
        hold on;
        %     plot(nanmean(lowW_decP), nanmean(highW_decP),'ro')%...
        axis equal; axis square;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        axis tight
        maxLim = max(max(xlim), max(ylim));
        minLim = min(min(xlim), min(ylim));
        set(gca, 'xlim', [minLim maxLim], 'ylim', [minLim maxLim]);
    subplot(235)
    set(gca, 'xlim', [minLim-2 maxLim+2], 'ylim', [minLim-2 maxLim+2]);
    line(xlim, ylim, 'color','k', 'linestyle','--');
    xlabel('Widths with position');
    ylabel('Widths with decoded position');
    catch
    end
end
% signrank(diffW)
subplot(232)
set(gca, 'xlim', [minLim-2 maxLim+2], 'ylim', [minLim-2 maxLim+2]);
line(xlim, ylim, 'color','k', 'linestyle','--');
xlabel('Widths @ low contrast');
ylabel('Widths @ high contrast');

%% confidence, error, accuracy

err_all.low = [];
err_all.high = [];
err_all.norm = [];

conf_all.low = [];
conf_all.high = [];
conf_all.norm = [];
conf_all.gray = [];




figure;
for idx = 1:length(Posterior_all)
    es = Posterior_all(idx).data;
    
    minXLim = 1;
    maxXLim = 32;
    t_l = (Posterior_all(idx).X_low>minXLim & Posterior_all(idx).X_low<maxXLim);
    t_h = (Posterior_all(idx).X_high>minXLim & Posterior_all(idx).X_high<maxXLim);
    t_n = (Posterior_all(idx).X_norm>minXLim & Posterior_all(idx).X_norm<maxXLim);
    
    Posterior_all(idx).X_low_orig = Posterior_all(idx).X_low_orig(1:length(Posterior_all(idx).accuracy.low_orig));
    Posterior_all(idx).X_high_orig = Posterior_all(idx).X_high_orig(1:length(Posterior_all(idx).accuracy.high_orig));
    t_lo = (Posterior_all(idx).X_low_orig>minXLim & Posterior_all(idx).X_low_orig<maxXLim);
    t_ho = (Posterior_all(idx).X_high_orig>minXLim & Posterior_all(idx).X_high_orig<maxXLim);
    
    acc_high(idx) = nanmean(Posterior_all(idx).accuracy.high(t_h));
    acc_norm(idx) = nanmean(Posterior_all(idx).accuracy.norm(t_n));
    acc_low(idx)  = nanmean(Posterior_all(idx).accuracy.low(t_l));
    
    err_high(idx) = nanmean(abs(Posterior_all(idx).error.high(t_h)));
    err_norm(idx) = nanmean(abs(Posterior_all(idx).error.norm(t_n)));
    err_low(idx)  = nanmean(abs(Posterior_all(idx).error.low(t_l)));
    
    err_all.high = [err_all.high (abs(Posterior_all(idx).error.high(t_h)))];
    err_all.norm = [err_all.norm (abs(Posterior_all(idx).error.norm(t_n)))];
    err_all.low  = [err_all.low (abs(Posterior_all(idx).error.low(t_l)))];
    
    conf_all.high = [conf_all.high (abs(Posterior_all(idx).confidence.high(t_h)))];
    conf_all.norm = [conf_all.norm (abs(Posterior_all(idx).confidence.norm(t_n)))];
    conf_all.low  = [conf_all.low (abs(Posterior_all(idx).confidence.low(t_l)))];
    conf_all.gray = [conf_all.gray (abs(Posterior_all(idx).confidence.gray))];
    
    conf_high(idx) = nanmean(Posterior_all(idx).confidence.high(t_h));
    conf_norm(idx) = nanmean(Posterior_all(idx).confidence.norm(t_n));
    conf_low(idx)  = nanmean(Posterior_all(idx).confidence.low(t_l));
    conf_gray(idx)  = nanmean(Posterior_all(idx).confidence.gray);
    
    drift_low(idx)  = nanmean(abs(diff(Posterior_all(idx).MAP.low)));
    drift_norm(idx) = nanmean(abs(diff(Posterior_all(idx).MAP.norm)));
    drift_high(idx) = nanmean(abs(diff(Posterior_all(idx).MAP.high)));
    drift_gray(idx) = nanmean(abs(diff(Posterior_all(idx).MAP.gray)));
    
    acc_higho(idx) = nanmean(Posterior_all(idx).accuracy.high_orig(t_ho));
    acc_lowo(idx)  = nanmean(Posterior_all(idx).accuracy.low_orig(t_lo));
    err_higho(idx) = nanmean(abs(Posterior_all(idx).error.high_orig(t_ho)));
    err_lowo(idx)  = nanmean(abs(Posterior_all(idx).error.low_orig(t_lo)));
    conf_higho(idx) = nanmean(Posterior_all(idx).confidence.high_orig(t_ho));
    conf_lowo(idx)  = nanmean(Posterior_all(idx).confidence.low_orig(t_lo));
    
    [~,p_acc(idx)] = ttest2((Posterior_all(idx).accuracy.high(t_h)),(Posterior_all(idx).accuracy.low(t_l)));
    [~,p_err(idx)] = ttest2((Posterior_all(idx).error.high(t_h)),(Posterior_all(idx).error.low(t_l)));
    [~,p_conf(idx)]= ttest2((Posterior_all(idx).confidence.high(t_h)),(Posterior_all(idx).confidence.low(t_l)));
    [~,p_drift(idx)]= ttest2(abs(diff(Posterior_all(idx).MAP.high)),abs(diff(Posterior_all(idx).MAP.low)));
    
    [~,p_acco(idx)] = ttest2((Posterior_all(idx).accuracy.high_orig(t_ho)),(Posterior_all(idx).accuracy.low_orig(t_lo)));
    [~,p_erro(idx)] = ttest2((Posterior_all(idx).error.high_orig(t_ho)),(Posterior_all(idx).error.low_orig(t_lo)));
    [~,p_confo(idx)]= ttest2((Posterior_all(idx).confidence.high_orig(t_ho)),(Posterior_all(idx).confidence.low_orig(t_lo)));
    
    X1 = 0:0.0001:1;
    X2 = 0:2:100;
    X3 = 0:50;
    subplot(4,8,idx)
    
    c_low(idx,:) = hist(Posterior_all(idx).confidence.low(t_l),X1);
    c_high(idx,:) = hist(Posterior_all(idx).confidence.high(t_h),X1);
    c_norm(idx,:) = hist(Posterior_all(idx).confidence.norm(t_n),X1);
    c_gray(idx,:) = hist(Posterior_all(idx).confidence.gray,X1);
    
    c_low(idx,:) = cumsum(c_low(idx,:)./sum(c_low(idx,:)));
    c_norm(idx,:) = cumsum(c_norm(idx,:)./sum(c_norm(idx,:)));
    c_high(idx,:) = cumsum(c_high(idx,:)./sum(c_high(idx,:)));
    c_gray(idx,:) = cumsum(c_gray(idx,:)./sum(c_gray(idx,:)));
    
    xlims = max([min(find(c_low(idx,:)>0.99)) min(find(c_norm(idx,:)>0.99)) min(find(c_high(idx,:)>0.99))]);
    plot(X1,c_low(idx,:),'b',X1,c_norm(idx,:),'k',X1,c_high(idx,:),'m',X1,c_gray(idx,:),'c')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    set(gca,'YLim',[0 1], 'XLim', [0 X1(xlims)])
    title(['Confidence ' num2str(Posterior_all(idx).iseries)]);
    
    e_low(idx,:) = hist(abs(Posterior_all(idx).error.low(t_l)),X2);
    e_high(idx,:) = hist(abs(Posterior_all(idx).error.high(t_h)),X2);
    e_norm(idx,:) = hist(abs(Posterior_all(idx).error.norm(t_n)),X2);
    
    e_low(idx,:) = cumsum(e_low(idx,:)./sum(e_low(idx,:)));
    e_norm(idx,:) = cumsum(e_norm(idx,:)./sum(e_norm(idx,:)));
    e_high(idx,:) = cumsum(e_high(idx,:)./sum(e_high(idx,:)));
    
    subplot(4,8,8+idx)
    plot(X2,e_low(idx,:),'b',X2,e_norm(idx,:),'k',X2,e_high(idx,:),'m')
    xlims = max([min(find(e_low(idx,:)>0.99)) min(find(e_norm(idx,:)>0.99)) min(find(e_high(idx,:)>0.99))]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    set(gca,'YLim',[0 1], 'XLim', [0 X2(xlims)])
    title(['Error ' num2str(Posterior_all(idx).iseries) ' ' num2str(sum(e_high(idx,:)-e_low(idx,:)))]);
    
    subplot(4,8,24+idx)
    plot(X2, cumsum(e_high(idx,:)-e_low(idx,:)),'k','linewidth',2)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    set(gca, 'XLim', [0 X2(xlims)])
    line(xlim, [0 0], 'color','k','linestyle','--');
    title(['\Delta Error ' num2str(Posterior_all(idx).iseries)]);
    
    subplot(4,8,16+idx)
    d_low(idx,:) = hist(abs(diff(Posterior_all(idx).MAP.low)),X3);
    d_high(idx,:) = hist(abs(diff(Posterior_all(idx).MAP.high)),X3);
    d_norm(idx,:) = hist(abs(diff(Posterior_all(idx).MAP.norm)),X3);
    d_gray(idx,:) = hist(abs(diff(Posterior_all(idx).MAP.gray)),X3);
    
    d_low(idx,:) = cumsum(d_low(idx,:)./sum(d_low(idx,:)));
    d_norm(idx,:) = cumsum(d_norm(idx,:)./sum(d_norm(idx,:)));
    d_high(idx,:) = cumsum(d_high(idx,:)./sum(d_high(idx,:)));
    d_gray(idx,:) = cumsum(d_gray(idx,:)./sum(d_gray(idx,:)));
    
    plot(X3,d_low(idx,:),'b',X3,d_norm(idx,:),'k',X3,d_high(idx,:),'m',X3,d_gray(idx,:),'c')
    xlims = max([min(find(d_low(idx,:)>0.99)) min(find(d_norm(idx,:)>0.99)) min(find(d_high(idx,:)>0.99))]);
    ylims_l = min([min(d_low(idx,:)) min(d_norm(idx,:)) min(d_high(idx,:))]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    set(gca,'YLim',[ylims_l 1], 'XLim', [0 X3(xlims)])
    title(['Drift ' num2str(Posterior_all(idx).iseries)]);
    
    for n = 1:50
        t_pts = (Posterior_all(idx).X_norm == n) & (es.outcome(Posterior_all(idx).t_norm)==1);
        confByPos.early(idx,n) = nanmean([abs(Posterior_all(idx).confidence.norm(t_pts)) nan]);
        erroByPos.early(idx,n) = nanmean([abs(Posterior_all(idx).error.norm(t_pts)) nan]);
        
        t_pts = (Posterior_all(idx).X_norm == n) & (es.outcome(Posterior_all(idx).t_norm)==2);
        confByPos.correct(idx,n) = nanmean([abs(Posterior_all(idx).confidence.norm(t_pts)) nan]);
        erroByPos.correct(idx,n) = nanmean([abs(Posterior_all(idx).error.norm(t_pts)) nan]);
        
        t_pts = (Posterior_all(idx).X_norm == n) & (es.outcome(Posterior_all(idx).t_norm)==3);
        confByPos.late(idx,n) = nanmean([abs(Posterior_all(idx).confidence.norm(t_pts)) nan]);
        erroByPos.late(idx,n) = nanmean([abs(Posterior_all(idx).error.norm(t_pts)) nan]);
        
        t_pts = (Posterior_all(idx).X_norm == n) & (es.outcome(Posterior_all(idx).t_norm)==4);
        confByPos.miss(idx,n) = nanmean([abs(Posterior_all(idx).confidence.norm(t_pts)) nan]);
        erroByPos.miss(idx,n) = nanmean([abs(Posterior_all(idx).error.norm(t_pts)) nan]);
    end
    confByPos.early(idx,:) = normalise1var(confByPos.early(idx,:));
    %     erroByPos.early(idx,:) = normalise1var(erroByPos.early(idx,:));
    confByPos.correct(idx,:) = normalise1var(confByPos.correct(idx,:));
    %     erroByPos.correct(idx,:) = normalise1var(erroByPos.correct(idx,:));
    confByPos.late(idx,:) = normalise1var(confByPos.late(idx,:));
    %     erroByPos.late(idx,:) = normalise1var(erroByPos.late(idx,:));
    confByPos.miss(idx,:) = normalise1var(confByPos.miss(idx,:));
    %     erroByPos.miss(idx,:) = normalise1var(erroByPos.miss(idx,:));
end

err =   [err_low' err_norm' err_high'];
conf =  [conf_low' conf_norm' conf_high'];
acc =   [acc_low' acc_norm' acc_high'];

erro =   [err_lowo' err_higho'];
confo =  [conf_lowo' conf_higho'];
acco =   [acc_lowo' acc_higho'];

err_remap_l =   [err_low' err_lowo'];
err_remap_h =   [err_high' err_higho'];
conf_remap_l =  [conf_low' conf_lowo'];
conf_remap_h =  [conf_high' conf_higho'];

figure;
subplot(221)
k_l = hist(abs(err_all.low),X2);
k_h = hist(abs(err_all.high),X2);
k_n = hist(abs(err_all.norm),X2);
k_l = cumsum(k_l./sum(k_l));
k_n = cumsum(k_n./sum(k_n));
k_h = cumsum(k_h./sum(k_h));

% plot(X2,k_l,'b',X2,k_n,'k',X2,k_h,'m')
errorarea_as(X2,nanmean(e_low),nansem(e_low),'b',1);
errorarea_as(X2,nanmean(e_high),nansem(e_high),'m',1);
plot(X2, nanmean(e_norm),'k--')
% ,X2,k_n,'k',X2,k_h,'m')
xlims = max([min(find(k_l>0.99)) min(find(k_n>0.99)) min(find(k_h>0.99))]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'YLim',[0 1], 'XLim', [0 X2(xlims)])
title(['Error Avg']);
ylabel('Cumulative fraction')
axis square

subplot(222)
k_l = hist(conf_all.low,X1);
k_h = hist(conf_all.high,X1);
k_n = hist(conf_all.norm,X1);
k_g = hist(conf_all.gray,X1);
k_l = cumsum(k_l./sum(k_l));
k_n = cumsum(k_n./sum(k_n));
k_h = cumsum(k_h./sum(k_h));
k_g = cumsum(k_g./sum(k_g));
xlims = max([min(find(k_l>0.99)) min(find(k_n>0.99)) min(find(k_h>0.99))]);
% plot(X1,k_l,'b',X1,k_n,'k',X1,k_h,'m',X1,k_g,'c')
errorarea_as(X1,nanmean(c_low),nansem(c_low),'b',1);
errorarea_as(X1,nanmean(c_high),nansem(c_high),'m',1);
plot(X1, nanmean(c_norm),'k--')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'YLim',[0 1], 'XLim', [0 X1(xlims)])
title(['Confidence Avg']);
axis square

subplot(223)
plot(err(:,1), err(:,3), 'ko');
axis equal;
axis([5 20 5 20])
line(xlim, ylim, 'color', 'k', 'linestyle','--')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlabel('Error @ low contrast')
ylabel('Error @ high contrast')

subplot(224)
plot(conf(:,1), conf(:,3), 'ko');
axis equal;
lims = max(conf(:));
axis([0 lims 0 lims])
line(xlim, ylim, 'color', 'k', 'linestyle','--')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlabel('Confidence @ low contrast')
ylabel('Confidence @ high contrast')

%%
% commenting out for V1
% display('commenting out for V1')

figure;
subplot(245)
bar([1 2], mean(err_remap_h,1),'b');
hold on;
plot(([1 2]')*ones(1,size(err_remap_h,1)),err_remap_h','o-','color',[.5 .5 .5])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'XTick',[1 2], 'XTickLabel', [{'Train base'}, {'Train high'}])
ylabel('Mean absolute error (cm)')
subplot(246)
bar([1 2], mean(err_remap_l,1),'b');
hold on;
plot(([1 2]')*ones(1,size(err_remap_l,1)),err_remap_l','o-','color',[.5 .5 .5])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'XTick',[1 2], 'XTickLabel', [{'Train base'}, {'Train low'}])
ylabel('Mean absolute error (cm)')

subplot(247)
bar([1 2], mean(erro,1),'b');
hold on;
plot(([1 2]')*ones(1,size(erro,1)),erro','o-','color',[.5 .5 .5])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'XTick',[1 2], 'XTickLabel', [{'Low'}, {'High'}])
ylabel('Mean absolute error (cm)')
subplot(248)
bar([1 2], mean(confo,1),'b');
hold on;
plot(([1 2]')*ones(1,size(confo,1)),confo','o-','color',[.5 .5 .5])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'XTick',[1 2], 'XTickLabel', [{'Low'}, {'High'}])
ylabel('Confidence')

subplot(231)
bar([1 2 3], mean(err,1),'b');
hold on;
plot(([1 2 3]')*ones(1,size(err,1)),err','o-','color',[.5 .5 .5])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'XTick',[1 2 3], 'XTickLabel', [{'Low'}, {'Base'}, {'High'}])
ylabel('Mean absolute error (cm)')

subplot(232)
bar([1 2 3], mean(conf,1),'b');
hold on;
plot(([1 2 3]')*ones(1,size(conf,1)),conf','o-','color',[.5 .5 .5])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'XTick',[1 2 3], 'XTickLabel', [{'Low'}, {'Base'}, {'High'}])
ylabel('Confidence')

subplot(233)
bar([1 2 3], mean(acc,1),'b');
hold on;
plot(([1 2 3]')*ones(1,size(acc,1)),acc','o-','color',[.5 .5 .5])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca,'XTick',[1 2 3], 'XTickLabel', [{'Low'}, {'Base'}, {'High'}])
ylabel('Accuracy')


%% Decoding performance
rescale = 1;

if ~rescale
for n = 1:length(Posterior_all)
    SummaryMaps.low(n,  :,:) = Posterior_all(n).meanPost.low;
    SummaryMaps.high(n, :,:) = Posterior_all(n).meanPost.high;
    SummaryMaps.norm(n, :,:) = Posterior_all(n).meanPost.norm;
    SummaryMaps.low_orig(n,  :,:) = Posterior_all(n).meanPost.low_orig;
    SummaryMaps.high_orig(n, :,:) = Posterior_all(n).meanPost.high_orig;
    
    SummaryMaps.allcond.low.early(n,:,:)    = Posterior_all(n).meanPost_new.low.early;
    SummaryMaps.allcond.low.correct(n,:,:)  = Posterior_all(n).meanPost_new.low.correct;
    SummaryMaps.allcond.low.late(n,:,:)     = Posterior_all(n).meanPost_new.low.late;
    SummaryMaps.allcond.low.miss(n,:,:)     = Posterior_all(n).meanPost_new.low.miss;
    
    SummaryMaps.allcond.high.early(n,:,:)    = Posterior_all(n).meanPost_new.high.early;
    SummaryMaps.allcond.high.correct(n,:,:)  = Posterior_all(n).meanPost_new.high.correct;
    SummaryMaps.allcond.high.late(n,:,:)     = Posterior_all(n).meanPost_new.high.late;
    SummaryMaps.allcond.high.miss(n,:,:)     = Posterior_all(n).meanPost_new.high.miss;
    
    SummaryMaps.allcond.norm.early(n,:,:)    = Posterior_all(n).meanPost_new.norm.early;
    SummaryMaps.allcond.norm.correct(n,:,:)  = Posterior_all(n).meanPost_new.norm.correct;
    SummaryMaps.allcond.norm.late(n,:,:)     = Posterior_all(n).meanPost_new.norm.late;
    SummaryMaps.allcond.norm.miss(n,:,:)     = Posterior_all(n).meanPost_new.norm.miss;
    
    SummaryMaps.allcond.all.early(n,:,:)    = Posterior_all(n).meanPost_new.all.early;
    SummaryMaps.allcond.all.correct(n,:,:)  = Posterior_all(n).meanPost_new.all.correct;
    SummaryMaps.allcond.all.late(n,:,:)     = Posterior_all(n).meanPost_new.all.late;
    SummaryMaps.allcond.all.miss(n,:,:)     = Posterior_all(n).meanPost_new.all.miss;
end
else
for n = 1:length(Posterior_all)
    SummaryMaps.low(n,  :,:) = (Posterior_all(n).meanPost.low);
    SummaryMaps.high(n, :,:) = (Posterior_all(n).meanPost.high);
    SummaryMaps.norm(n, :,:) = (Posterior_all(n).meanPost.norm);
    SummaryMaps.low_orig(n,  :,:) = (Posterior_all(n).meanPost.low_orig);
    SummaryMaps.high_orig(n, :,:) = (Posterior_all(n).meanPost.high_orig);
    
    SummaryMaps.allcond.low.early(n,:,:)    = (Posterior_all(n).meanPost_new.low.early);
    SummaryMaps.allcond.low.correct(n,:,:)  = (Posterior_all(n).meanPost_new.low.correct);
    SummaryMaps.allcond.low.late(n,:,:)     = (Posterior_all(n).meanPost_new.low.late);
    SummaryMaps.allcond.low.miss(n,:,:)     = (Posterior_all(n).meanPost_new.low.miss);
    
    SummaryMaps.allcond.high.early(n,:,:)    = (Posterior_all(n).meanPost_new.high.early);
    SummaryMaps.allcond.high.correct(n,:,:)  = (Posterior_all(n).meanPost_new.high.correct);
    SummaryMaps.allcond.high.late(n,:,:)     = (Posterior_all(n).meanPost_new.high.late);
    SummaryMaps.allcond.high.miss(n,:,:)     = (Posterior_all(n).meanPost_new.high.miss);
    
    SummaryMaps.allcond.norm.early(n,:,:)    = (Posterior_all(n).meanPost_new.norm.early);
    SummaryMaps.allcond.norm.correct(n,:,:)  = (Posterior_all(n).meanPost_new.norm.correct);
    SummaryMaps.allcond.norm.late(n,:,:)     = (Posterior_all(n).meanPost_new.norm.late);
    SummaryMaps.allcond.norm.miss(n,:,:)     = (Posterior_all(n).meanPost_new.norm.miss);
    
    SummaryMaps.allcond.all.early(n,:,:)    = (Posterior_all(n).meanPost_new.all.early);
    SummaryMaps.allcond.all.correct(n,:,:)  = (Posterior_all(n).meanPost_new.all.correct);
    SummaryMaps.allcond.all.late(n,:,:)     = (Posterior_all(n).meanPost_new.all.late);
    SummaryMaps.allcond.all.miss(n,:,:)     = (Posterior_all(n).meanPost_new.all.miss);
end
end
bins = Posterior_all(1).decoder.bins;


figure;
subplot(131)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.low,1))));
axis xy; colorbar; axis tight; axis equal; axis tight
title(['Low Avg'], 'fontsize', 14);
subplot(133)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.norm,1))));
axis xy; colorbar; axis tight; axis equal; axis tight
title('Normal Avg', 'fontsize', 14);
subplot(132)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.high,1))));
axis xy; colorbar; axis tight; axis equal; axis tight
title(['High Avg'], 'fontsize', 14);
RedWhiteBlue;
for n = 1:3
    subplot(1,3,n)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
end
for n = 1:3
    subplot(1,3,n)
    xlabel('Original Position')
    ylabel('Decoded Posterior')
    set(gca,'CLim',[-0.5 0.5])
    hcb = colorbar('YTick',[-0.5 0 0.5],'YTickLabel',{'2^-0.5 x', 'Chance', '2^0.5 x'});
    set(hcb,'YTickMode','manual')
end

figure;
subplot(121)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.low_orig,1))));
axis xy; colorbar; axis tight; axis equal; axis tight
title(['Low Avg'], 'fontsize', 14);

subplot(122)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.high_orig,1))));
axis xy; colorbar; axis tight; axis equal; axis tight
title(['High Avg'], 'fontsize', 14);
RedWhiteBlue;
for n = 1:2
    subplot(1,2,n)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
end
for n = 1:2
    subplot(1,2,n)
    xlabel('Original Position')
    ylabel('Decoded Posterior')
    set(gca,'CLim',[-0.5 0.5])
    hcb = colorbar('YTick',[-0.5 0 0.5],'YTickLabel',{'2^-0.5 x', 'Chance', '2^0.5 x'});
    set(hcb,'YTickMode','manual')
end

%%
figure;

subplot(4,4,1)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.low.early,1))));
subplot(4,4,2)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.low.correct,1))));
subplot(4,4,3)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.low.late,1))));
subplot(4,4,4)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.low.miss,1))));

subplot(4,4,5)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.norm.early,1))));
subplot(4,4,6)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.norm.correct,1))));
subplot(4,4,7)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.norm.late,1))));
subplot(4,4,8)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.norm.miss,1))));

subplot(4,4,9)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.high.early,1))));
subplot(4,4,10)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.high.correct,1))));
subplot(4,4,11)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.high.late,1))));
subplot(4,4,12)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.high.miss,1))));

subplot(4,4,13)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.all.early,1))));
subplot(4,4,14)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.all.correct,1))));
subplot(4,4,15)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.all.late,1))));
subplot(4,4,16)
imagesc(bins, bins, real(squeeze(nanmean(SummaryMaps.allcond.all.miss,1))));

RedWhiteBlue;
for n = 1:16
    subplot(4,4,n)
    axis xy; axis equal; axis tight
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
    %     xlabel('Original Position')
    %     ylabel('Decoded Posterior')
    set(gca,'CLim',[-0.5 0.5])
    %     axis off;
    %     hcb = colorbar('YTick',[-0.5 0 0.5],'YTickLabel',{'2^-0.5 x', 'Chance', '2^0.5 x'});
    %     set(hcb,'YTickMode','manual')
end
subplot(4,4,1); title('EARLY');
ylabel('Low')
subplot(4,4,2); title('CORRECT');
subplot(4,4,3); title('LATE');
subplot(4,4,4); title('MISS');
subplot(4,4,5)
ylabel('Medium')
subplot(4,4,9)
ylabel('High')
subplot(4,4,13)
ylabel('All')
% figure;
% for n = 1:8
%     subplot(2,8,n)
%     imagesc(bins, bins, Posterior_all(n).meanPost.low);
%     axis xy; colorbar; axis tight; axis equal; axis tight
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
%     line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
%     set(gca,'CLim',[-0.5 0.5])
%     subplot(2,8,n+8)
%     imagesc(bins, bins, Posterior_all(n).meanPost.high);
%     axis xy; colorbar; axis tight; axis equal; axis tight
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
%     line(xlim, ylim, 'color','k','linewidth',1, 'linestyle','--')
%     set(gca,'CLim',[-0.5 0.5])
% end

%% errors, beh & decoding
try
    numBins = size(Posterior_all(1).Posterior_norm,2);
    figure
    subplot(132)
    bar([1 2 3], mean(err*50/numBins,1),'r');
    hold on;
%     plot(([1 2 3]')*ones(1,8),err','o-','color',[.5 .5 .5])
    errorbar([1 2 3], mean(err*50/numBins,1), sem(err*50/numBins),'k')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    set(gca,'XTick',[1 2 3], 'XTickLabel', [{'Low'}, {'Base'}, {'High'}])
    ylabel('Mean absolute error (cm)')
    set(gca,'YLim',[0 18])
    
    err2 = [min([erro(:,1) err(:,1)]')', min([erro(:,2) err(:,3)]')'];
    subplot(133)
    bar([1 2], mean(err2*50/numBins,1),'r');
    hold on;
%     plot(([1 2 3]')*ones(1,8),err','o-','color',[.5 .5 .5])
    errorbar([1 2], mean(err2*50/numBins,1), sem(err2*50/numBins),'k')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    set(gca,'XTick',[1 2], 'XTickLabel', [{'Low'}, {'High'}])
    ylabel('Mean absolute error (cm)')
    set(gca,'YLim',[0 18])
    
    subplot(131)
    bar([1 2 3], mean(100-beh,1),'b');
    hold on;
%     plot(([1 2 3]')*ones(1,8),100-beh','o-','color',[.5 .5 .5])
    errorbar([1 2 3], mean(100-beh,1), sem(beh),'k')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
    set(gca,'XTick',[1 2 3], 'XTickLabel', [{'Low'}, {'Base'}, {'High'}])
    ylabel('Behavioral errors (%)')
    set(gca,'YLim',[0 40])
    
catch
end

figure;
imagesc(bins, bins, (real(log2(50*(squeeze(nanmean(SummaryMaps.norm,1)))))));
axis xy; axis equal; axis tight
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
set(gca, 'XTick', [0 50 100], 'YTick', [0 50 100])
colorbar; set(gca, 'clim', [-.5 0.5])
line(xlim, ylim,'linestyle','--','color','k');
xlabel('Actual Position (cm)');
ylabel('Decoded Position (cm)');
RedWhiteBlue;
% figure;
% thres = 0.01;
% peakRate_all = [];
% for animalIdx = 1:length(Posterior_all)
%     subsetCells = nanmean(Posterior_all(animalIdx).decoder_high.model.EV)>thres &...
%         nanmean(Posterior_all(animalIdx).decoder_low.model.EV)>thres; 
%     
%     peakRate_high = 60*max(Posterior_all(animalIdx).decoder_high.model.meanModel(subsetCells,:));
%     peakRate_low  = 60*max(Posterior_all(animalIdx).decoder_low.model.meanModel(subsetCells,:));
% 
%     peakRate_all = [peakRate_all peakRate_high-peakRate_low];
%     
%     errorbarxy(nanmean(peakRate_low),nanmean(peakRate_high),...
%         nansem(peakRate_low),nansem(peakRate_high),...
%         {'ko','k','k'});
%     hold on;axis equal;
%     axis([0 150 0 150]);
% %     pause
% end
% line(xlim, ylim, 'color','k', 'linestyle','--');
    