figure('Position',[ 1          41        1920        1084])
X = -50:50;
lowD = 1;
highD = 20;
drange = [(50-highD):51-lowD 51+lowD:(51+highD)];
drange2= [1:51-lowD 51+lowD:100];
smth_win = 0;
gradient_n = 1;
diff_amount = 3;

low_lim = 1;
high_lim = 32;

clear gradients diffs
% for diff_amount = 1:30
clear varOnSpace;
for n = 1:length(Posterior_all)
    maxLim = min([size(Posterior_all(n).Posterior_low,1) ...
        size(Posterior_all(n).Posterior_norm,1) ...
        size(Posterior_all(n).Posterior_high,1)]);
    
    gradients(n).low  = gradient(Posterior_all(n).MAP.low, gradient_n);
    gradients(n).norm = gradient(Posterior_all(n).MAP.norm,gradient_n);
    gradients(n).high = gradient(Posterior_all(n).MAP.high,gradient_n);
    
    [~, MAP.low_orig] = max(Posterior_all(n).Posterior_low_orig,[],2);
    [~, MAP.high_orig] = max(Posterior_all(n).Posterior_high_orig,[],2);
    
    gradients(n).low_o  = gradient(MAP.low_orig, gradient_n);
    gradients(n).high_o = gradient(MAP.high_orig,gradient_n);
    
    diffs(n).low  = Posterior_all(n).MAP.low  - circshift(Posterior_all(n).MAP.low', diff_amount)';
    diffs(n).norm = Posterior_all(n).MAP.norm - circshift(Posterior_all(n).MAP.norm',diff_amount)';
    diffs(n).high = Posterior_all(n).MAP.high - circshift(Posterior_all(n).MAP.high',diff_amount)';
    
    diffs(n).low_o  = MAP.low_orig - circshift(MAP.low_orig, diff_amount);
    diffs(n).high_o = MAP.high_orig - circshift(MAP.high_orig, diff_amount);
    
    dots(n).low  = sum(Posterior_all(n).Posterior_low .* circshift(Posterior_all(n).Posterior_low, diff_amount),2) ...
        ./ 	(sqrt(sum(Posterior_all(n).Posterior_low.^2,2)) .* sqrt(sum(circshift(Posterior_all(n).Posterior_low, diff_amount).^2,2)));
    dots(n).norm = sum(Posterior_all(n).Posterior_norm .* circshift(Posterior_all(n).Posterior_norm,diff_amount),2) ...
        ./ 	(sqrt(sum(Posterior_all(n).Posterior_norm.^2,2)) .* sqrt(sum(circshift(Posterior_all(n).Posterior_norm, diff_amount).^2,2)));
    dots(n).high = sum(Posterior_all(n).Posterior_high .* circshift(Posterior_all(n).Posterior_high,diff_amount),2) ...
        ./ 	(sqrt(sum(Posterior_all(n).Posterior_high.^2,2)) .* sqrt(sum(circshift(Posterior_all(n).Posterior_high, diff_amount).^2,2)));
    
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
    imagesc(Posterior_all(n).Posterior_low(1:maxLim,:)')
    hold on; axis xy
%     plot(Posterior_all(n).X_low(1:maxLim),'k')
    plot(smthInTime(Posterior_all(n).MAP.low(1:maxLim),60,smth_win),'m')
%     set(gca,'CLim',[-0.5 0.5]);
    title(['Animal #:' num2str(n) ', Low Contrast'])
    
    t_low = Posterior_all(n).X_low>low_lim & Posterior_all(n).X_low<high_lim;
    t_high = Posterior_all(n).X_high>low_lim & Posterior_all(n).X_high<high_lim;
    
    subplot(3,5,[4])
    low_hist = hist(gradient(Posterior_all(n).MAP.low(Posterior_all(n).X_low<high_lim),gradient_n),X);
    low_hist = low_hist./sum(low_hist);
%     plot(X,low_hist); 
    plot(Posterior_all(n).X_low(t_low), diffs(n).low(t_low), '.')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    xlabel('Position')
    ylabel('Diff from prev time')
    title('Low Contrast')
    
    subplot(3,5,[6 7 8])
    hold off; 
    imagesc(Posterior_all(n).Posterior_norm(1:maxLim,:)')
    hold on; axis xy
%     plot(Posterior_all(n).X_norm(1:maxLim),'k')
    plot(smthInTime(Posterior_all(n).MAP.norm(1:maxLim),60,smth_win),'m')
%     set(gca,'CLim',[-0.5 0.5]);
    title(['Animal #:' num2str(n) ', Medium Contrast'])
    subplot(3,5,[9 ])%9
    norm_hist = hist(gradient(Posterior_all(n).MAP.norm(Posterior_all(n).X_norm<high_lim),gradient_n),X);
    norm_hist = norm_hist./sum(norm_hist);
%     plot(X,norm_hist); 
    plot(Posterior_all(n).X_norm, diffs(n).norm, '.')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    xlabel('Position')
    ylabel('Diff from prev time')
    title('Normal Contrast')
    
    subplot(3,5,[11 12 13])
    hold off; 
    imagesc(Posterior_all(n).Posterior_high(1:maxLim,:)')
    hold on; axis xy
%     plot(Posterior_all(n).X_high(1:maxLim),'k')
    plot(smthInTime(Posterior_all(n).MAP.high(1:maxLim),60,smth_win),'m')
%     set(gca,'CLim',[-0.5 0.5]);
    title(['Animal #:' num2str(n) ', High Contrast'])
    
    subplot(3,5,[5]) % 15
    high_hist = hist(gradient(Posterior_all(n).MAP.high(Posterior_all(n).X_high<high_lim),gradient_n),X);
    high_hist = high_hist./sum(high_hist);
%     plot(X,high_hist); 
    plot(Posterior_all(n).X_high(t_high), diffs(n).high(t_high), '.')
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
    
    X_low  = Posterior_all(n).X_low;
    X_norm = Posterior_all(n).X_norm;
    X_high = Posterior_all(n).X_high;
    
    low_v(n)  = nanvar(gradients(n).low( ...
         X_low<high_lim & X_low>low_lim));
    norm_v(n)  = nanvar(gradients(n).norm( ...
        X_norm<high_lim & X_norm>low_lim));
    high_v(n)  = nanvar(gradients(n).high( ...
        X_high<high_lim & X_high>low_lim));
    
    low_vo(n)  = nanvar(gradients(n).low_o( ...
        Posterior_all(n).X_low_orig>low_lim  & Posterior_all(n).X_low_orig<high_lim ));
    high_vo(n)  = nanvar(gradients(n).high_o( ...
        Posterior_all(n).X_high_orig>low_lim & Posterior_all(n).X_high_orig<high_lim ));
    
    
    
    low_vd(n)  = nanvar(gradients(n).low( ...
        X_low<high_lim & X_low>low_lim));
    norm_vd(n)  = nanvar(gradients(n).norm( ...
        X_norm<high_lim & X_norm>low_lim));
    high_vd(n)  = nanvar(gradients(n).high( ...
        X_high<high_lim & X_high>low_lim));

    %     low_v(n)  = nanvar(gradient(Posterior_all(n).MAP.low(Posterior_all(n).X_low<high_lim & Posterior_all(n).X_low   >1),gradient_n));
%     norm_v(n) = nanvar(gradient(Posterior_all(n).MAP.norm(Posterior_all(n).X_norm<high_lim & Posterior_all(n).X_norm>1),gradient_n));
%     high_v(n) = nanvar(gradient(Posterior_all(n).MAP.high(Posterior_all(n).X_high<high_lim & Posterior_all(n).X_high>1),gradient_n));
    
%     low_vd(n)  = nanvar(diff(Posterior_all(n).MAP.low(Posterior_all(n).X_low<high_lim)));
%     norm_vd(n) = nanvar(diff(Posterior_all(n).MAP.norm(Posterior_all(n).X_norm<high_lim)));
%     high_vd(n) = nanvar(diff(Posterior_all(n).MAP.high(Posterior_all(n).X_high<high_lim)));
    low_vd(n)  = nanvar(diffs(n).low( ...
        Posterior_all(n).X_low<high_lim & Posterior_all(n).X_low>low_lim));
    norm_vd(n)  = nanvar(diffs(n).norm( ...
        Posterior_all(n).X_norm<high_lim & Posterior_all(n).X_norm>low_lim));
    high_vd(n)  = nanvar(diffs(n).high( ...
        Posterior_all(n).X_high<high_lim & Posterior_all(n).X_high>low_lim));
    
    low_vdo(n)  = nanvar(diffs(n).low_o( ...
        Posterior_all(n).X_low_orig<high_lim & Posterior_all(n).X_low_orig>low_lim));
    high_vdo(n)  = nanvar(diffs(n).high_o( ...
        Posterior_all(n).X_high_orig<high_lim & Posterior_all(n).X_high_orig>low_lim));
    
    posBins = [1 32 38 50];
    for posIdx = 1:length(posBins)-1
        varOnSpace.low(n,posIdx) = nanvar(gradients(n).low( ...
            X_low>=posBins(posIdx) & X_low<posBins(posIdx+1)));
        varOnSpace.norm(n,posIdx) = nanvar(gradients(n).norm( ...
            X_norm>=posBins(posIdx) & X_norm<posBins(posIdx+1)));
        varOnSpace.high(n,posIdx) = nanvar(gradients(n).high( ...
            X_high>=posBins(posIdx) & X_high<posBins(posIdx+1)));
    end
    
    low_dot(n) = mean(dots(n).low);
    norm_dot(n) = mean(dots(n).norm);
    high_dot(n) = mean(dots(n).high);
    axis tight
    set(gca,'XLim',[0 0.1])
    drawnow
    figure('Position',[ 1          41        1920        1084])
%     pause
end
% close
% figure;

% figure;
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

