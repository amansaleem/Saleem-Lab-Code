function [o, numOver, numUnder] = varianceWithLick(Posterior_all, numBins)

Posterior_all = checkDrift_slide(Posterior_all, numBins, 1, 1, 0);

lowLim = 30;
highLim = 50;

allVar.nl = [];
allVar.la = [];
allVar.lc = [];
allVar.le = [];
allVar.lr = [];

for ianimal = 1:length(Posterior_all)
    es = Posterior_all(ianimal).data;
    P = Posterior_all(ianimal);
    rewZonePost = nan*ones(size(es.traj));
    
    rewZonePost(P.t_norm) = (nanmean((P.Posterior_norm(:,33:38)),2));
    normPost = nanmean(rewZonePost(P.t_norm & es.outcome==2));
    rewZonePost(P.t_norm) = (nanmean((P.Posterior_norm(:,33:38)),2))./abs(normPost);
    rewZonePost(P.t_low)  = (nanmean((P.Posterior_low(:,33:38)),2))./abs(normPost);
    rewZonePost(P.t_high) = (nanmean((P.Posterior_high(:,33:38)),2))./abs(normPost);
    
    decPos = nan*ones(size(es.traj));
    decPos(P.t_norm) = (P.MAP.norm);
    decPos(P.t_low)  = (P.MAP.low);
    decPos(P.t_high) = (P.MAP.high);
    
    varAll = nan*ones(size(es.traj));
    varAll(P.t_norm) = P.slideVar.norm;
    varAll(P.t_low) = P.slideVar.low;
    varAll(P.t_high) = P.slideVar.high;
    
    o.varNoLick{ianimal}    = (varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick==0));
    o.varWithLickC{ianimal} = (varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0 & es.outcome==2));
    o.varWithLickE{ianimal} = (varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0 & es.outcome~=2));
    o.varWithLick{ianimal} = (varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0));
    o.varWithLickR{ianimal} = (varAll(es.traj>65 & es.traj<75 ...
        & es.lick>0 & es.outcome==2));
    numOver(ianimal) = 0;
    numUnder(ianimal) = 0;
    posBins = 1:10:101;
    for iPos = 1:length(posBins)-1
        k = es.traj>=posBins(iPos) & es.traj<posBins(iPos+1) ...
            & es.outcome>0;
        meanPosVar(ianimal,iPos) = nanmedian(varAll(k & ...
            es.lick==0 & es.outcome>0));
        semPosVar(ianimal,iPos) = nansem(varAll(k & ...
            es.lick==0 & es.outcome>0));
        if posBins(iPos)>0 && posBins(iPos+1)<60 %  &&
                        numOver(ianimal)  = numOver(ianimal) + nansum(es.lick(varAll(k)>meanPosVar(ianimal,iPos)));
                        numUnder(ianimal) = numUnder(ianimal) + nansum(es.lick(varAll(k)<meanPosVar(ianimal,iPos)));
%             numOver(ianimal)  = numOver(ianimal) + nansum((es.lick(k) & varAll(k)>meanPosVar(ianimal,iPos)));
%             numUnder(ianimal) = numUnder(ianimal) + nansum((es.lick(k) & varAll(k)<meanPosVar(ianimal,iPos)));
        end
    end
    posBins = posBins(1:end-1);
    
    allVar.nl = [allVar.nl o.varNoLick{ianimal}'];
    allVar.la = [allVar.la o.varWithLick{ianimal}'];
    allVar.lc = [allVar.lc o.varWithLickC{ianimal}'];
    allVar.le = [allVar.le o.varWithLickE{ianimal}'];
    allVar.lr = [allVar.lr o.varWithLickR{ianimal}'];
    
    o.meanvarNoLick(ianimal) = nanmean(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick==0));
    o.stdvarNoLick(ianimal) = nanstd(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick==0));
    o.semvarNoLick(ianimal) = nansem(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick==0));
    
    o.meanvarWithLickC(ianimal) = nanmean(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0 & es.outcome==2));
    o.stdvarWithLickC(ianimal) = nanstd(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0 & es.outcome==2));
    o.semvarWithLickC(ianimal) = nansem(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0 & es.outcome==2));
    
    o.meanvarWithLickE(ianimal) = nanmean(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0 & es.outcome~=2));
    o.stdvarWithLickE(ianimal) = nanstd(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0 & es.outcome~=2));
    o.semvarWithLickE(ianimal) = nansem(varAll(es.traj>lowLim & es.traj<highLim ...
        & es.lick>0 & es.outcome~=2));
    
    figure(1)
    subplot(3,3,ianimal)
    hold off;
    semilogx(rewZonePost(es.traj>lowLim & es.traj<highLim & es.lick==0),...
        varAll(es.traj>lowLim & es.traj<highLim & es.lick==0),'k');
    hold on;
    semilogx(rewZonePost(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome~=2), ...
        varAll(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome~=2),'r.');
    semilogx(rewZonePost(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome==2), ...
        varAll(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome==2),'b.');
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    axis tight;
    
    figure(4)
    subplot(3,3,ianimal)
    hold off;
    errorarea_as(posBins, meanPosVar(ianimal,:), semPosVar(ianimal,:),'k')
    hold on;
    earlyLicks = es.lick>0 & es.traj<50;
    plot(es.traj(earlyLicks),varAll(earlyLicks),'r.');
    axis tight;
    ylims = ylim;
    axis([0 50 0 ylims(2)]); 
    %     plot(2*decPos(es.traj>lowLim & es.traj<highLim & es.lick==0),...
    %         varAll(es.traj>lowLim & es.traj<highLim & es.lick==0),'k.');
    %     hold on;
    %     plot(2*decPos(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome~=2), ...
    %         varAll(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome~=2),'r.');
    %     plot(2*decPos(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome==2), ...
    %         varAll(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome==2),'b.');
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    
%     figure(2)
%     subplot(3,3,ianimal)
%     hold off;
%     %     errorarea_as(posBins, meanPosVar(ianimal,:), semPosVar(ianimal,:),'k')
%     %     hold on;
%     %     earlyLicks = es.lick>0 & es.traj<100;
%     %     plot(es.traj(earlyLicks),varAll(earlyLicks),'.');
%     plot(2*decPos(es.traj>lowLim & es.traj<highLim & es.lick==0),...
%         varAll(es.traj>lowLim & es.traj<highLim & es.lick==0),'k.');
%     hold on;
%     plot(2*decPos(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome~=2), ...
%         varAll(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome~=2),'r.');
%     plot(2*decPos(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome==2), ...
%         varAll(es.traj>lowLim & es.traj<highLim & es.lick>0 & es.outcome==2),'b.');
%     set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
%     axis tight;
    
    figure(3)
    subplot(3,3,ianimal)
    [~,X] = hist([o.varNoLick{ianimal}' o.varWithLickC{ianimal}' o.varWithLickE{ianimal}'],10);
    distr_nl = hist(o.varNoLick{ianimal}, X);
    distr_la = hist(o.varWithLick{ianimal}, X);
    distr_lc = hist(o.varWithLickC{ianimal}, X);
    distr_le = hist(o.varWithLickE{ianimal}, X);
    distr_lr = hist(o.varWithLickR{ianimal}, X);
    distr_nl = distr_nl./sum(distr_nl);
    distr_lc = distr_lc./sum(distr_lc);
    distr_la = distr_la./sum(distr_la);
    distr_le = distr_le./sum(distr_le);
    distr_lr = distr_lr./sum(distr_lr);
    
    hold off;
%     plot(X, distr_nl, 'k', X, distr_lc, 'b', X, distr_le, 'r', 'linewidth',2)
    plot(X, distr_nl, 'k', X, distr_la, 'r', 'linewidth',2)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
end

figure(3)
subplot(3,3,9)
[~,X] = hist([allVar.nl allVar.lc allVar.le allVar.lr],10);
distr_nl = hist(allVar.nl, X);
distr_lc = hist(allVar.lc, X);
distr_le = hist(allVar.le, X);
distr_lr = hist(allVar.lr, X);
distr_nl = distr_nl./sum(distr_nl);
distr_lc = distr_lc./sum(distr_lc);
distr_le = distr_le./sum(distr_le);
distr_lr = distr_lr./sum(distr_lr);

hold off;
% plot(X, distr_nl, 'k', X, distr_lc, 'b',...
%     X, distr_le, 'r', X, distr_lr, 'c', 'linewidth',2)
plot(X, distr_nl, 'k', X, distr_la, 'r',...
    X, distr_lr, 'c', 'linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Across all animals');
xlabel('Jitter in time')
ylabel('Number of occurances')
legend('No Lick', 'Correct trial lick', 'Error trial Lick', 'Rewarded licks')
legend('No Lick', 'Early Licks', 'Rewarded licks')