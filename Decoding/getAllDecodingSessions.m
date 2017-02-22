[Posterior_all_530_103_106] = plot_some_decoding('M140501_BALL', 530, [103:106], 'contrast');
[Posterior_all_531_103_106] = plot_some_decoding('M140501_BALL', 531, [103:106], 'contrast');
[Posterior_all_601_103_106] = plot_some_decoding('M140501_BALL', 601, [103:106], 'contrast');
[Posterior_all_602_102_106] = plot_some_decoding('M140501_BALL', 602, [102:106], 'contrast');
[Posterior_all_604_107_110] = plot_some_decoding('M140502_BALL', 604, [107:110], 'contrast');
[Posterior_all_918_103_105] = plot_some_decoding('M130918_BALL', 1030, [103:105], 'contrast');
[Posterior_all_920_102_103] = plot_some_decoding('M130920_BALL', 1025, [102:103], 'contrast');

Posterior_all(1) = Posterior_all_530_103_106;
Posterior_all(2) = Posterior_all_531_103_106;
Posterior_all(3) = Posterior_all_601_103_106;
Posterior_all(4) = Posterior_all_602_102_106;
Posterior_all(5) = Posterior_all_604_107_110;
Posterior_all(6) = Posterior_all_918_103_105;
Posterior_all(7) = Posterior_all_920_102_103;

%% calculating the mean
meanPost_norm = zeros(size(Posterior_all(1).meanPost.norm));
meanPost_high = meanPost_norm;
meanPost_low = meanPost_norm;

for n = 1:length(Posterior_all)
    meanPost_norm = meanPost_norm + Posterior_all(n).meanPost.norm;
    meanPost_high = meanPost_high + Posterior_all(n).meanPost.high;
    meanPost_low  = meanPost_low + Posterior_all(n).meanPost.low;
end

meanPost_norm = meanPost_norm./length(Posterior_all);
meanPost_high = meanPost_high./length(Posterior_all);
meanPost_low  = meanPost_low./length(Posterior_all);

%% calculating the confidence and accuracy
for idx = 1:length(Posterior_all)
    t = zeros(size(Posterior_all(idx).Posterior_norm));
    t_l = zeros(size(Posterior_all(idx).Posterior_low));
    t_h = zeros(size(Posterior_all(idx).Posterior_high));
    for n = 1:size(t_l,1)
        t_l(n,Posterior_all(idx).X_low(n)) = 1;
    end
    for n = 1:size(t_h,1)
        t_h(n,Posterior_all(idx).X_high(n)) = 1;
    end
    for n = 1:size(t,1)
        t(n,Posterior_all(idx).X_norm(n)) = 1;
    end
    Posterior_all(idx).confidence_dist.low = (max(Posterior_all(idx).Posterior_low') - min(Posterior_all(idx).Posterior_low'));
    Posterior_all(idx).confidence_dist.norm = (max(Posterior_all(idx).Posterior_norm')- min(Posterior_all(idx).Posterior_norm'));
    Posterior_all(idx).confidence_dist.high = (max(Posterior_all(idx).Posterior_high') - min(Posterior_all(idx).Posterior_high'));
    
    Posterior_all(idx).accuracy_dist.low = (Posterior_all(idx).Posterior_low(t_l>0));
    Posterior_all(idx).accuracy_dist.norm = (Posterior_all(idx).Posterior_norm(t>0));
    Posterior_all(idx).accuracy_dist.high = (Posterior_all(idx).Posterior_high(t_h>0));
    
    t_l = zeros(size(Posterior_all(idx).Posterior_low_orig));
    t_h = zeros(size(Posterior_all(idx).Posterior_high_orig));
    for n = 1:size(t_l,1)
        t_l(n,Posterior_all(idx).X_low_orig(n)) = 1;
    end
    for n = 1:size(t_h,1)
        t_h(n,Posterior_all(idx).X_high_orig(n)) = 1;
    end
    Posterior_all(idx).confidence_dist.low_orig = (max(Posterior_all(idx).Posterior_low_orig') - min(Posterior_all(idx).Posterior_low_orig'));
    Posterior_all(idx).confidence_dist.high_orig = (max(Posterior_all(idx).Posterior_high_orig') - min(Posterior_all(idx).Posterior_high_orig'));
    
    Posterior_all(idx).accuracy_dist.low_orig = (Posterior_all(idx).Posterior_low_orig(t_l>0));
    Posterior_all(idx).accuracy_dist.high_orig = (Posterior_all(idx).Posterior_high_orig(t_h>0));
end

%%
figure;
str = '_orig';
for n = 1:7
    subplot(2, 7, n)
    [~,XA] = hist([Posterior_all(n).accuracy_dist.low_orig Posterior_all(n).accuracy_dist.low_orig],10);
    [~,XC] = hist([Posterior_all(n).confidence_dist.low_orig Posterior_all(n).confidence_dist.low_orig],10);
    lowA = hist(Posterior_all(n).accuracy_dist.low_orig,XA);
    highA = hist(Posterior_all(n).accuracy_dist.high_orig,XA);
    plot(XA,lowA./sum(lowA),'b',XA,highA./sum(highA),'r','linewidth',2); axis tight
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none','YTick',[])
    subplot(2, 7, n+7)
    lowC = hist(Posterior_all(n).confidence_dist.low_orig,XC);
    highC = hist(Posterior_all(n).confidence_dist.high_orig,XC);
    plot(XC,lowC./sum(lowC),'b',XC,highC./sum(highC),'r','linewidth',2); axis tight
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none','YTick',[])
end
subplot(2,7,1)
ylabel('Accuracy')
subplot(2,7,8)
ylabel('Confidence')
