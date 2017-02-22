function checkDrift_onAggregate(Aggregate, slideBins, cal_var)

if nargin<2
    slideBins = 5;
end
if nargin<3
    cal_var = 1;
end

figure('Position',[ 1  41  1920  1084])
X = -50:50;
lowD = 1;
highD = 20;
drange = [(50-highD):51-lowD 51+lowD:(51+highD)];
drange2= [1:51-lowD 51+lowD:100];
smth_win = 0;
gradient_n = 1;
diff_amount = 1;

low_lim = 1;
high_lim = 32;


clear gradients diffs
% for diff_amount = 1:30
clear varOnSpace;
%% main calculations
for animalIdx = 1:length(Posterior_all)
    maxLim = min([size(Posterior_all(animalIdx).Posterior_low,1) ...
        size(Posterior_all(animalIdx).Posterior_norm,1) ...
        size(Posterior_all(animalIdx).Posterior_high,1)]);
    
    gradients(animalIdx).low  = gradient(Posterior_all(animalIdx).MAP.low, gradient_n);
    gradients(animalIdx).norm = gradient(Posterior_all(animalIdx).MAP.norm,gradient_n);
    gradients(animalIdx).high = gradient(Posterior_all(animalIdx).MAP.high,gradient_n);
    
    [~, MAP.low_orig] = max(Posterior_all(animalIdx).Posterior_low_orig,[],2);
    [~, MAP.high_orig] = max(Posterior_all(animalIdx).Posterior_high_orig,[],2);
    
    gradients(animalIdx).low_o  = gradient(MAP.low_orig, gradient_n);
    gradients(animalIdx).high_o = gradient(MAP.high_orig,gradient_n);
    
    diffs(animalIdx).low  = Posterior_all(animalIdx).MAP.low;%  - circshift(Posterior_all(animalIdx).MAP.low', diff_amount)';
    diffs(animalIdx).norm = Posterior_all(animalIdx).MAP.norm;% - circshift(Posterior_all(animalIdx).MAP.norm',diff_amount)';
    diffs(animalIdx).high = Posterior_all(animalIdx).MAP.high;% - circshift(Posterior_all(animalIdx).MAP.high',diff_amount)';
    
    diffs(animalIdx).low(Posterior_all(animalIdx).MAP.low==1) = nan;%  - circshift(Posterior_all(animalIdx).MAP.low', diff_amount)';
    diffs(animalIdx).norm(Posterior_all(animalIdx).MAP.norm==1) = nan;% - circshift(Posterior_all(animalIdx).MAP.norm',diff_amount)';
    diffs(animalIdx).high(Posterior_all(animalIdx).MAP.high==1) = nan;% - circshift(Posterior_all(animalIdx).MAP.high',diff_amount)';
    
    diffs(animalIdx).low(Posterior_all(animalIdx).MAP.low==50) = nan;%  - circshift(Posterior_all(animalIdx).MAP.low', diff_amount)';
    diffs(animalIdx).norm(Posterior_all(animalIdx).MAP.norm==50) = nan;% - circshift(Posterior_all(animalIdx).MAP.norm',diff_amount)';
    diffs(animalIdx).high(Posterior_all(animalIdx).MAP.high==50) = nan;% - circshift(Posterior_all(animalIdx).MAP.high',diff_amount)';
    
    diffs(animalIdx).X_low  = Posterior_all(animalIdx).X_low;%  - circshift(Posterior_all(animalIdx).X_low', diff_amount)';
    diffs(animalIdx).X_norm = Posterior_all(animalIdx).X_norm;% - circshift(Posterior_all(animalIdx).X_norm',diff_amount)';
    diffs(animalIdx).X_high = Posterior_all(animalIdx).X_high;% - circshift(Posterior_all(animalIdx).X_high',diff_amount)';

    diffs(animalIdx).low_o  = MAP.low_orig - circshift(MAP.low_orig, diff_amount);
    diffs(animalIdx).high_o = MAP.high_orig - circshift(MAP.high_orig, diff_amount);
    
    slideVar(animalIdx).low_o  = (abs(diffs(animalIdx).low_o));
    slideVar(animalIdx).high_o = (abs(diffs(animalIdx).high_o));
    
    if cal_var
        slideVar(animalIdx).low  = ((diffs(animalIdx).low));
        slideVar(animalIdx).high = ((diffs(animalIdx).high));
        slideVar(animalIdx).norm = ((diffs(animalIdx).norm));
        
        slideVar(animalIdx).X_low  = ((diffs(animalIdx).X_low));
        slideVar(animalIdx).X_high = ((diffs(animalIdx).X_high));
        slideVar(animalIdx).X_norm = ((diffs(animalIdx).X_norm));
        
        for iSlide = 1:slideBins
            slideVar(animalIdx).low  = [slideVar(animalIdx).low'   circshift((diffs(animalIdx).low') , iSlide)]';
            slideVar(animalIdx).norm = [slideVar(animalIdx).norm'  circshift((diffs(animalIdx).norm'), iSlide)]';
            slideVar(animalIdx).high = [slideVar(animalIdx).high'  circshift((diffs(animalIdx).high'), iSlide)]';
            
            slideVar(animalIdx).X_low  = [slideVar(animalIdx).X_low   circshift((diffs(animalIdx).X_low) , iSlide)];
            slideVar(animalIdx).X_norm = [slideVar(animalIdx).X_norm  circshift((diffs(animalIdx).X_norm), iSlide)];
            slideVar(animalIdx).X_high = [slideVar(animalIdx).X_high  circshift((diffs(animalIdx).X_high), iSlide)];
            
            slideVar(animalIdx).low_o  = slideVar(animalIdx).low_o  + circshift(abs(diffs(animalIdx).low_o') , iSlide)';
            slideVar(animalIdx).high_o = slideVar(animalIdx).high_o + circshift(abs(diffs(animalIdx).high_o'), iSlide)';
        end
        slideVar(animalIdx).low  = nanstd(slideVar(animalIdx).low);
        slideVar(animalIdx).norm = nanstd(slideVar(animalIdx).norm);
        slideVar(animalIdx).high = nanstd(slideVar(animalIdx).high);
        
        slideVar(animalIdx).X_low  = nanstd(slideVar(animalIdx).X_low');
        slideVar(animalIdx).X_norm = nanstd(slideVar(animalIdx).X_norm');
        slideVar(animalIdx).X_high = nanstd(slideVar(animalIdx).X_high');
    else
        slideVar(animalIdx).low  = (abs(diffs(animalIdx).low));
        slideVar(animalIdx).norm = (abs(diffs(animalIdx).norm));
        slideVar(animalIdx).high = (abs(diffs(animalIdx).high));
        
        
        for iSlide = 1:slideBins
            slideVar(animalIdx).low  = slideVar(animalIdx).low  + circshift(abs(diffs(animalIdx).low') , iSlide)';
            slideVar(animalIdx).norm = slideVar(animalIdx).norm + circshift(abs(diffs(animalIdx).norm'), iSlide)';
            slideVar(animalIdx).high = slideVar(animalIdx).high + circshift(abs(diffs(animalIdx).high'), iSlide)';
            
            slideVar(animalIdx).low_o  = slideVar(animalIdx).low_o  + circshift(abs(diffs(animalIdx).low_o') , iSlide)';
            slideVar(animalIdx).high_o = slideVar(animalIdx).high_o + circshift(abs(diffs(animalIdx).high_o'), iSlide)';
        end
        slideVar(animalIdx).low  = slideVar(animalIdx).low./slideBins;
        slideVar(animalIdx).norm = slideVar(animalIdx).norm./slideBins;
        slideVar(animalIdx).high = slideVar(animalIdx).high./slideBins;
    end
    slideVar(animalIdx).low_o  = slideVar(animalIdx).low_o./slideBins;
    slideVar(animalIdx).high_o = slideVar(animalIdx).high_o./slideBins;
    
    
    dots(animalIdx).low  = sum(Posterior_all(animalIdx).Posterior_low .* circshift(Posterior_all(animalIdx).Posterior_low, diff_amount),2) ...
        ./ 	(sqrt(sum(Posterior_all(animalIdx).Posterior_low.^2,2)) .* sqrt(sum(circshift(Posterior_all(animalIdx).Posterior_low, diff_amount).^2,2)));
    dots(animalIdx).norm = sum(Posterior_all(animalIdx).Posterior_norm .* circshift(Posterior_all(animalIdx).Posterior_norm,diff_amount),2) ...
        ./ 	(sqrt(sum(Posterior_all(animalIdx).Posterior_norm.^2,2)) .* sqrt(sum(circshift(Posterior_all(animalIdx).Posterior_norm, diff_amount).^2,2)));
    dots(animalIdx).high = sum(Posterior_all(animalIdx).Posterior_high .* circshift(Posterior_all(animalIdx).Posterior_high,diff_amount),2) ...
        ./ 	(sqrt(sum(Posterior_all(animalIdx).Posterior_high.^2,2)) .* sqrt(sum(circshift(Posterior_all(animalIdx).Posterior_high, diff_amount).^2,2)));
    
    %     dots(n).low_o  = MAP.low_orig  .* circshift(MAP.low_orig, diff_amount);
    %     dots(n).high_o = MAP.high_orig .* circshift(MAP.high_orig, diff_amount);
    
    %     diffs(n).low(abs(diffs(n).low)>25)   =  50-diffs(n).low(abs(diffs(n).low)>25);
    %     diffs(n).norm(abs(diffs(n).norm)>25) =  50-diffs(n).norm(abs(diffs(n).norm)>25);
    %     diffs(n).high(abs(diffs(n).high)>25) =  50-diffs(n).high(abs(diffs(n).high)>25);
    %
    %     diffs(n).low_o(abs(diffs(n).low_o)>25)   =  50-diffs(n).low_o(abs(diffs(n).low_o)>25);
    %     diffs(n).high_o(abs(diffs(n).high_o)>25) =  50-diffs(n).high_o(abs(diffs(n).high_o)>25);
    
    subplot(3,5,[1 2 3])
    hold off;
    imagesc(Posterior_all(animalIdx).Posterior_low(1:maxLim,:)')
    hold on; axis xy
    %     plot(Posterior_all(n).X_low(1:maxLim),'k')
    plot(smthInTime(Posterior_all(animalIdx).MAP.low(1:maxLim),60,smth_win),'m')
    %     set(gca,'CLim',[-0.5 0.5]);
    title(['Animal #:' num2str(animalIdx) ', Low Contrast'])
    
    t_low = Posterior_all(animalIdx).X_low>low_lim & Posterior_all(animalIdx).X_low<high_lim;
    t_high = Posterior_all(animalIdx).X_high>low_lim & Posterior_all(animalIdx).X_high<high_lim;
    
    subplot(3,5,[4])
    low_hist = hist(gradient(Posterior_all(animalIdx).MAP.low(Posterior_all(animalIdx).X_low<high_lim),gradient_n),X);
    low_hist = low_hist./sum(low_hist);
    %     plot(X,low_hist);
    plot(Posterior_all(animalIdx).X_low(t_low), diffs(animalIdx).low(t_low), '.')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    xlabel('Position')
    ylabel('Diff from prev time')
    title('Low Contrast')
    
    subplot(3,5,[6 7 8])
    hold off;
    imagesc(Posterior_all(animalIdx).Posterior_norm(1:maxLim,:)')
    hold on; axis xy
    %     plot(Posterior_all(n).X_norm(1:maxLim),'k')
    plot(smthInTime(Posterior_all(animalIdx).MAP.norm(1:maxLim),60,smth_win),'m')
    %     set(gca,'CLim',[-0.5 0.5]);
    title(['Animal #:' num2str(animalIdx) ', Medium Contrast'])
    subplot(3,5,[9 ])%9
    norm_hist = hist(gradient(Posterior_all(animalIdx).MAP.norm(Posterior_all(animalIdx).X_norm<high_lim),gradient_n),X);
    norm_hist = norm_hist./sum(norm_hist);
    %     plot(X,norm_hist);
    plot(Posterior_all(animalIdx).X_norm, diffs(animalIdx).norm, '.')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    xlabel('Position')
    ylabel('Diff from prev time')
    title('Normal Contrast')
    
    subplot(3,5,[11 12 13])
    hold off;
    imagesc(Posterior_all(animalIdx).Posterior_high(1:maxLim,:)')
    hold on; axis xy
    %     plot(Posterior_all(n).X_high(1:maxLim),'k')
    plot(smthInTime(Posterior_all(animalIdx).MAP.high(1:maxLim),60,smth_win),'m')
    %     set(gca,'CLim',[-0.5 0.5]);
    title(['Animal #:' num2str(animalIdx) ', High Contrast'])
    
    subplot(3,5,[5]) % 15
    high_hist = hist(gradient(Posterior_all(animalIdx).MAP.high(Posterior_all(animalIdx).X_high<high_lim),gradient_n),X);
    high_hist = high_hist./sum(high_hist);
    %     plot(X,high_hist);
    plot(Posterior_all(animalIdx).X_high(t_high), diffs(animalIdx).high(t_high), '.')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    xlabel('Position')
    ylabel('Diff from prev time')
    title('High Contrast')
    
    subplot(3,5,10)
    plot(low_hist(drange), X(drange),'b',norm_hist(drange),X(drange), 'k', ...
        high_hist(drange),X(drange), 'r');
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    %     RedWhiteBlue;
    %     low(n)  = sum(low_hist(drange2));
    %     norm(n) = sum(norm_hist(drange2));
    %     high(n) = sum(high_hist(drange2));
    
    X_low  = Posterior_all(animalIdx).X_low;
    X_norm = Posterior_all(animalIdx).X_norm;
    X_high = Posterior_all(animalIdx).X_high;
    
    X_low_o  = Posterior_all(animalIdx).X_low_orig;
    X_high_o = Posterior_all(animalIdx).X_high_orig;
    
    low_v(animalIdx)  = nanvar(gradients(animalIdx).low( ...
        X_low<high_lim & X_low>low_lim));
    norm_v(animalIdx)  = nanvar(gradients(animalIdx).norm( ...
        X_norm<high_lim & X_norm>low_lim));
    high_v(animalIdx)  = nanvar(gradients(animalIdx).high( ...
        X_high<high_lim & X_high>low_lim));
    
    low_vo(animalIdx)  = nanvar(gradients(animalIdx).low_o( ...
        Posterior_all(animalIdx).X_low_orig>low_lim  & Posterior_all(animalIdx).X_low_orig<high_lim ));
    high_vo(animalIdx)  = nanvar(gradients(animalIdx).high_o( ...
        Posterior_all(animalIdx).X_high_orig>low_lim & Posterior_all(animalIdx).X_high_orig<high_lim ));
    
    
    low_vd(animalIdx)  = nanvar(gradients(animalIdx).low( ...
        X_low<high_lim & X_low>low_lim));
    norm_vd(animalIdx)  = nanvar(gradients(animalIdx).norm( ...
        X_norm<high_lim & X_norm>low_lim));
    high_vd(animalIdx)  = nanvar(gradients(animalIdx).high( ...
        X_high<high_lim & X_high>low_lim));
    
    %     low_v(n)  = nanvar(gradient(Posterior_all(n).MAP.low(Posterior_all(n).X_low<high_lim & Posterior_all(n).X_low   >1),gradient_n));
    %     norm_v(n) = nanvar(gradient(Posterior_all(n).MAP.norm(Posterior_all(n).X_norm<high_lim & Posterior_all(n).X_norm>1),gradient_n));
    %     high_v(n) = nanvar(gradient(Posterior_all(n).MAP.high(Posterior_all(n).X_high<high_lim & Posterior_all(n).X_high>1),gradient_n));
    
    %     low_vd(n)  = nanvar(diff(Posterior_all(n).MAP.low(Posterior_all(n).X_low<high_lim)));
    %     norm_vd(n) = nanvar(diff(Posterior_all(n).MAP.norm(Posterior_all(n).X_norm<high_lim)));
    %     high_vd(n) = nanvar(diff(Posterior_all(n).MAP.high(Posterior_all(n).X_high<high_lim)));
    low_vd(animalIdx)  = nanvar(diffs(animalIdx).low( ...
        Posterior_all(animalIdx).X_low<high_lim & Posterior_all(animalIdx).X_low>low_lim));
    norm_vd(animalIdx)  = nanvar(diffs(animalIdx).norm( ...
        Posterior_all(animalIdx).X_norm<high_lim & Posterior_all(animalIdx).X_norm>low_lim));
    high_vd(animalIdx)  = nanvar(diffs(animalIdx).high( ...
        Posterior_all(animalIdx).X_high<high_lim & Posterior_all(animalIdx).X_high>low_lim));
    
    low_vdo(animalIdx)  = nanvar(diffs(animalIdx).low_o( ...
        Posterior_all(animalIdx).X_low_orig<high_lim & Posterior_all(animalIdx).X_low_orig>low_lim));
    high_vdo(animalIdx)  = nanvar(diffs(animalIdx).high_o( ...
        Posterior_all(animalIdx).X_high_orig<high_lim & Posterior_all(animalIdx).X_high_orig>low_lim));
    
    posBins = [1:2:50];
    for posIdx = 1:length(posBins)-1
        meanSlide.low (animalIdx,posIdx) = 2*nanmean(slideVar(animalIdx).low(...
            X_low  >=posBins(posIdx)& X_low<posBins(posIdx+1)));
        meanSlide.norm(animalIdx,posIdx) = 2*nanmean(slideVar(animalIdx).norm(...
            X_norm>=posBins(posIdx) & X_norm<posBins(posIdx+1)));
        meanSlide.high(animalIdx,posIdx) = 2*nanmean(slideVar(animalIdx).high(...
            X_high>=posBins(posIdx) & X_high<posBins(posIdx+1)));
        
        meanSlide.X_low (animalIdx,posIdx) = 2*nanmean(slideVar(animalIdx).X_low(...
            X_low  >=posBins(posIdx)& X_low<posBins(posIdx+1)));
        meanSlide.X_norm(animalIdx,posIdx) = 2*nanmean(slideVar(animalIdx).X_norm(...
            X_norm>=posBins(posIdx) & X_norm<posBins(posIdx+1)));
        meanSlide.X_high(animalIdx,posIdx) = 2*nanmean(slideVar(animalIdx).X_high(...
            X_high>=posBins(posIdx) & X_high<posBins(posIdx+1)));
        
        meanSlide.low_o (animalIdx,posIdx) = 2*nanmean(slideVar(animalIdx).low_o(...
            X_low_o  >=posBins(posIdx)& X_low_o<posBins(posIdx+1)));
        meanSlide.high_o(animalIdx,posIdx) = 2*nanmean(slideVar(animalIdx).high_o(...
            X_high_o>=posBins(posIdx) & X_high_o<posBins(posIdx+1)));
    end
    
    %     posBins = [1 32 38 50];
    for posIdx = 1:length(posBins)-1
        varOnSpace.low(animalIdx,posIdx) = nanvar(gradients(animalIdx).low( ...
            X_low>=posBins(posIdx) & X_low<posBins(posIdx+1)));
        varOnSpace.norm(animalIdx,posIdx) = nanvar(gradients(animalIdx).norm( ...
            X_norm>=posBins(posIdx) & X_norm<posBins(posIdx+1)));
        varOnSpace.high(animalIdx,posIdx) = nanvar(gradients(animalIdx).high( ...
            X_high>=posBins(posIdx) & X_high<posBins(posIdx+1)));
        
        meanconf.low(animalIdx,posIdx) = nanmean(Posterior_all(animalIdx).confidence.low( ...
            X_low>=posBins(posIdx) & X_low<posBins(posIdx+1)));
        meanconf.norm(animalIdx,posIdx) = nanmean(Posterior_all(animalIdx).confidence.norm( ...
            X_norm>=posBins(posIdx) & X_norm<posBins(posIdx+1)));
        meanconf.high(animalIdx,posIdx) = nanmean(Posterior_all(animalIdx).confidence.high( ...
            X_high>=posBins(posIdx) & X_high<posBins(posIdx+1)));
        
        meanerr.low(animalIdx,posIdx) = nanmean(abs(Posterior_all(animalIdx).error.low( ...
            X_low>=posBins(posIdx) & X_low<posBins(posIdx+1))));
        meanerr.norm(animalIdx,posIdx) = nanmean(abs(Posterior_all(animalIdx).error.norm( ...
            X_norm>=posBins(posIdx) & X_norm<posBins(posIdx+1))));
        meanerr.high(animalIdx,posIdx) = nanmean(abs(Posterior_all(animalIdx).error.high( ...
            X_high>=posBins(posIdx) & X_high<posBins(posIdx+1))));
    end
    
    low_dot(animalIdx) = mean(dots(animalIdx).low);
    norm_dot(animalIdx) = mean(dots(animalIdx).norm);
    high_dot(animalIdx) = mean(dots(animalIdx).high);
    axis tight
    set(gca,'XLim',[0 0.1])
    
    subplot(3,5,14)
    plot(Posterior_all(animalIdx).X_low, slideVar(animalIdx).low, '.')
    hold on;
    plot(posBins(1:end-1), meanSlide.low(animalIdx,:),'r')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    hold off
    xlabel('Position')
    ylabel('Sliding mean diff')
    title('Low Contrast')
    axis([1 50 0 max([slideVar(animalIdx).low slideVar(animalIdx).high])]);
    
    subplot(3,5,15)
    plot(Posterior_all(animalIdx).X_high, slideVar(animalIdx).high, '.')
    hold on;
    plot(posBins(1:end-1), meanSlide.high(animalIdx,:),'r')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    hold off
    xlabel('Position')
    ylabel('Sliding mean diff')
    title('High Contrast')
    axis([1 50 0 max([slideVar(animalIdx).low slideVar(animalIdx).high])]);
    
    
    %     figure('Position',[ 1          41        1920        1084])
    %     pause
%     [hists(animalIdx).all, varBins] = hist([slideVar(animalIdx).low slideVar(animalIdx).norm slideVar(animalIdx).high],10);
%     hists(animalIdx).low  = hist(slideVar(animalIdx).low,varBins);
%     hists(animalIdx).norm = hist(slideVar(animalIdx).norm,varBins);
%     hists(animalIdx).high = hist(slideVar(animalIdx).high,varBins);
%     
%     lickVar = [slideVar(animalIdx).low(Posterior_all(animalIdx).t_low(Posterior_all(animalIdx).data.lick>0)) ...
%         slideVar(animalIdx).high(Posterior_all(animalIdx).t_high(Posterior_all(animalIdx).data.lick>0)) ...
%         slideVar(animalIdx).norm(Posterior_all(animalIdx).t_norm(Posterior_all(animalIdx).data.lick>0))];
%    
%     hists(animalIdx).allLicks = hist(lickVar,varBins);
%     hists(animalIdx).varBins  = varBins;
%     
%     subplot(3,5,9)
%     plot(varBins, hists(animalIdx).all./sum(hists(animalIdx).all), 'k', ...
%         varBins, hists(animalIdx).allLicks./sum(hists(animalIdx).allLicks), 'c');
    drawnow
%     pause
end
% close
% figure;

%% figure;
hold off
subplot(221)
hold off
plot(low_v, high_v, 'ko')
hold on;
% plot(low_v([4 6]), high_v([4 6]), 'k*')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis square
axis equal
lims = [0 nanmax([low_v high_v])];
axis([lims lims]);
line(lims, lims, 'linestyle','--','color','k');
xlabel('Variance of gradient (low contrast)');
ylabel('Variance of gradient (high contrast)');
title(['Gradient over ' num2str(round(16.67*gradient_n)) 'ms']);

subplot(222)
plot(low_vd, high_vd, 'ko')
hold on;
% plot(low_vd([4 6]), high_vd([4 6]), 'k*')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis square
axis equal
lims = [0 nanmax([low_vd high_vd])];
axis([lims lims]);
line(lims, lims, 'linestyle','--','color','k');
xlabel('Variance of diff (low contrast)');
ylabel('Variance of diff (high contrast)');
title(['Circular assumption; diff over ' num2str(round(diff_amount*16.67)) 'ms']  );
% title(['Gradient over +- ' num2str(gradient_n) ' points']);

subplot(224)
plot(low_vdo, high_vdo, 'ko')
hold on;
% plot(low_vd([4 6]), high_vd([4 6]), 'k*')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis square
axis equal
lims = [0 nanmax([low_vdo high_vdo])];
axis([lims lims]);
line(lims, lims, 'linestyle','--','color','k');
xlabel('Variance of diff orig (low contrast)');
ylabel('Variance of diff orig (high contrast)');
title(['Circular assumption; diff over ' num2str(round(diff_amount*16.67)) 'ms']  );
% title(['Gradient over +- ' num2str(gradient_n) ' points']);

subplot(223)
hold off
plot(low_vo, high_vo, 'ko')
hold on;
% plot(low_v([4 6]), high_v([4 6]), 'k*')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis square
axis equal
lims = [0 nanmax([low_vo high_vo])];
axis([lims lims]);
line(lims, lims, 'linestyle','--','color','k');
xlabel('Variance of gradient orig (low contrast)');
ylabel('Variance of gradient orig (high contrast)');
title(['Gradient over ' num2str(round(16.67*gradient_n)) 'ms']);

% pause(1/2)
% end

figure;
plot(low_dot, high_dot, 'ko')
hold on;
% plot(low_v([4 6]), high_v([4 6]), 'k*')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis square
axis equal
lims = [min([low_dot high_dot])-0.1 max([low_dot high_dot])+0.1];
axis([lims lims]);
line(xlim, ylim, 'linestyle','--','color','k');
xlabel('Dot prod (low contrast)');
ylabel('Dot prod (high contrast)');
title(['Over ' num2str(round(16.67*diff_amount)) 'ms']);

%% Variance over time in position
figure
for animalIdx = 1:size(meanSlide.low,1)
    subplot(3,3,animalIdx)
    hold off
    plot(2*posBins(1:end-1), meanSlide.high(animalIdx,:),'r')
    hold on;
    plot(2*posBins(1:end-1), meanSlide.low(animalIdx,:),'b')
    line([20 20],ylim, 'color','k', 'linestyle','--');
    line([40 40],ylim, 'color','k', 'linestyle','--');
    line([60 60],ylim, 'color','k', 'linestyle','--');
    line([80 80],ylim, 'color','k', 'linestyle','--');
    line([65 65],ylim, 'color','c', 'linestyle','-');
    line([75 75],ylim, 'color','c', 'linestyle','-');
    %     hold on; plot(posBins(1:end-1), meanSlide.norm(n,:),'k')
    %     plot(posBins(1:end-1), ,'r')
    axis tight
    line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanSlide.low(n,:) meanSlide.high(n,:)])]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    title(num2str(animalIdx))
end
subplot(3,3,9)

plot(2*posBins(1:end-1), mean(meanSlide.high,1),'r')
hold on;
plot(2*posBins(1:end-1), mean(meanSlide.low,1),'b')
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanSlide.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanSlide.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary, variance')
%% Plot stuff, variance
figure(102)
subplot(221)
% errorarea_as(2*posBins(1:end-1), mean(meanSlide.high,1),sem(meanSlide.high),'r')
% hold on;
% errorarea_as(2*posBins(1:end-1), mean(meanSlide.low,1),sem(meanSlide.low),'b')
% errorarea_as(2*posBins(1:end-1), mean(meanSlide.norm,1),sem(meanSlide.norm),'k')
plot(2*posBins(1:end-1), mean(meanSlide.high,1),'r', 'linewidth',1.2)
hold on;
plot(2*posBins(1:end-1), mean(meanSlide.low,1),'b', 'linewidth',1.2)
plot(2*posBins(1:end-1), mean(meanSlide.norm,1),'k', 'linewidth',1.2)

plot(2*posBins(1:end-1), mean(meanSlide.X_high,1),'r--', 'linewidth',1)
plot(2*posBins(1:end-1), mean(meanSlide.X_low,1),'b--', 'linewidth',2)
plot(2*posBins(1:end-1), mean(meanSlide.X_norm,1),'k--', 'linewidth',2)
axis tight
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanSlide.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
% line(xlim, [-0.5 -0.5], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanSlide.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary, variance')
ylabel('Mean abs diff')
xlabel('Position')
%% Difference in variance over time in position
figure
for animalIdx = 1:size(meanSlide.low,1)
    subplot(3,3,animalIdx)
    hold off
    plot(2*posBins(1:end-1), meanSlide.low(animalIdx,:) - meanSlide.high(animalIdx,:),'k')
    %     hold on; plot(posBins(1:end-1), meanSlide.norm(n,:),'k')
    %     plot(posBins(1:end-1), ,'r')
    axis tight
    line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanSlide.low(n,:) meanSlide.high(n,:)])]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    title(num2str(animalIdx))
end
subplot(3,3,9)
errorarea_as(2*posBins(1:end-1), mean(meanSlide.low-meanSlide.high,1), sem(meanSlide.low-meanSlide.high), 'r');
hold on;
errorarea_as(2*posBins(1:end-1), mean(meanSlide.norm-meanSlide.high,1), sem(meanSlide.norm-meanSlide.high), 'k');
% errorarea_as(posBins(1:end-1), mean(meanSlide.low_o-meanSlide.high_o,1), sem(meanSlide.low_o-meanSlide.high_o), 'b');
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanSlide.low(n,:) meanSlide.high(n,:)])]);
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary, diff in variance')

figure(102)
subplot(224)
errorarea_as(2*posBins(1:end-1), mean(meanSlide.low-meanSlide.high,1), sem(meanSlide.low-meanSlide.high), 'r');
% hold on;
% errorarea_as(2*posBins(1:end-1), mean(meanSlide.norm-meanSlide.high,1), sem(meanSlide.norm-meanSlide.high), 'k');
% errorarea_as(posBins(1:end-1), mean(meanSlide.low_o-meanSlide.high_o,1), sem(meanSlide.low_o-meanSlide.high_o), 'b');
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanSlide.low(n,:) meanSlide.high(n,:)])]);
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary, diff in variance')
ylabel('Diff in variance (low-high)')
xlabel('Position')
%% Confidence over time in position
figure
for animalIdx = 1:size(meanconf.low,1)
    subplot(3,3,animalIdx)
    hold off
    plot(2*posBins(1:end-1), meanconf.high(animalIdx,:),'r')
    hold on;
    plot(2*posBins(1:end-1), meanconf.low(animalIdx,:),'b')
    line([20 20],ylim, 'color','k', 'linestyle','--');
    line([40 40],ylim, 'color','k', 'linestyle','--');
    line([60 60],ylim, 'color','k', 'linestyle','--');
    line([80 80],ylim, 'color','k', 'linestyle','--');
    line([65 65],ylim, 'color','c', 'linestyle','-');
    line([75 75],ylim, 'color','c', 'linestyle','-');
    %     hold on; plot(posBins(1:end-1), meanconf.norm(n,:),'k')
    %     plot(posBins(1:end-1), ,'r')
    axis tight
    line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanconf.low(n,:) meanconf.high(n,:)])]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    title(num2str(animalIdx))
end
subplot(3,3,9)

errorarea_as(2*posBins(1:end-1), mean(meanconf.low - meanconf.high,1),sem(meanconf.low - meanconf.high),'k')
% hold on;
% plot(2*posBins(1:end-1), mean(meanconf.low,1),'b')
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanconf.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanconf.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary confidence')

figure(101)
subplot(221)
plot(2*posBins(1:end-1), mean(meanconf.norm, 1),'k')
hold on;
plot(2*posBins(1:end-1), mean(meanconf.low, 1),'b')
plot(2*posBins(1:end-1), mean(meanconf.high, 1),'r')
% plot(2*posBins(1:end-1), mean(meanconf.low,1),'b')
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanconf.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','-'); %axis([1 50 0 max([meanconf.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary confidence')
ylabel('Height of Posterior (Low-High)');


subplot(222)
errorarea_as(2*posBins(1:end-1), mean(meanconf.norm - meanconf.low,1),sem(meanconf.norm - meanconf.low),'b')
hold on;
errorarea_as(2*posBins(1:end-1), mean(meanconf.norm - meanconf.high,1),sem(meanconf.norm - meanconf.high),'r')
% plot(2*posBins(1:end-1), mean(meanconf.low,1),'b')
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanconf.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','-'); %axis([1 50 0 max([meanconf.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary confidence')
ylabel('Height of Posterior (Low-High)');


subplot(223)
errorarea_as(2*posBins(1:end-1), mean(meanconf.low - meanconf.high,1),sem(meanconf.low - meanconf.high),'k')
% hold on;
% plot(2*posBins(1:end-1), mean(meanconf.low,1),'b')
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanconf.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanconf.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary confidence')
ylabel('Height of Posterior (Low-High)');

figure(102)
subplot(223)
hold off;
errorarea_as(2*posBins(1:end-1), -mean(meanSlide.norm - meanSlide.low,1),sem(meanSlide.norm - meanSlide.low),'b')
hold on;
errorarea_as(2*posBins(1:end-1), -mean(meanSlide.norm - meanSlide.high,1),sem(meanSlide.norm - meanSlide.high),'r')
hold off;
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanSlide.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanSlide.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary, low - medium, high - medium')
ylabel('Mean abs diff')
xlabel('Position')
%% Error over time in position
figure
for animalIdx = 1:size(meanerr.low,1)
    subplot(3,3,animalIdx)
    hold off
    plot(2*posBins(1:end-1), meanerr.high(animalIdx,:),'r')
    hold on;
    plot(2*posBins(1:end-1), meanerr.low(animalIdx,:),'b')
    line([20 20],ylim, 'color','k', 'linestyle','--');
    line([40 40],ylim, 'color','k', 'linestyle','--');
    line([60 60],ylim, 'color','k', 'linestyle','--');
    line([80 80],ylim, 'color','k', 'linestyle','--');
    line([65 65],ylim, 'color','c', 'linestyle','-');
    line([75 75],ylim, 'color','c', 'linestyle','-');
    %     hold on; plot(posBins(1:end-1), meanerr.norm(n,:),'k')
    %     plot(posBins(1:end-1), ,'r')
    axis tight
    line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanerr.low(n,:) meanerr.high(n,:)])]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    title(num2str(animalIdx))
end
subplot(3,3,9)

errorarea_as(2*posBins(1:end-1), mean(meanerr.low - meanerr.high,1),sem(meanerr.low - meanerr.high),'k')
% hold on;
% plot(2*posBins(1:end-1), mean(meanerr.low,1),'b')
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanerr.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanerr.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary Error')

figure(102)
subplot(222)
errorarea_as(2*posBins(1:end-1), mean(meanerr.low - meanerr.high,1),sem(meanerr.low - meanerr.high),'k')
% hold on;
% plot(2*posBins(1:end-1), mean(meanerr.low,1),'b')
line([20 20],ylim, 'color','k', 'linestyle','--');
line([40 40],ylim, 'color','k', 'linestyle','--');
line([60 60],ylim, 'color','k', 'linestyle','--');
line([80 80],ylim, 'color','k', 'linestyle','--');
line([65 65],ylim, 'color','c', 'linestyle','-');
line([75 75],ylim, 'color','c', 'linestyle','-');
%     hold on; plot(posBins(1:end-1), meanerr.norm(n,:),'k')
%     plot(posBins(1:end-1), ,'r')
axis tight
line(xlim, [0 0], 'color','k', 'linestyle','--'); %axis([1 50 0 max([meanerr.low(n,:) meanSlide.high(n,:)])]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Summary Error')
ylabel('Mean abs error (Low-High)')
%% Certainty related to licking
