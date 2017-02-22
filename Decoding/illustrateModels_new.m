flag_norm = 0;
inputs = -49.9:0.1:50;
rng(2)
%% Lower firing rate
x_low = zeros(750,1000);
x_high = zeros(750,1000);
for n = 1:size(x_low,1)
    centre_low(n) = 500+round(20*randn(1));
    centre_high(n) = 500+round(20*randn(1));
    x_low(n,centre_low(n)) = 1;
    x_high(n,centre_high(n)) = 1;
    
    x_low(n,:) = (smthInTime(x_low(n,:),1000,88+15*randn(1)));
    x_high(n,:) = (smthInTime(x_high(n,:)/3,1000,88+15*randn(1)));
    
end

maxVal_low = mean(max(x_low,[],2));
maxVal_high = mean(max(x_high,[],2));
for n = 1:size(x_low,1)
    if flag_norm
        x_low(n,:) = normalise1var(x_low(n,:));
        x_high(n,:) = normalise1var(x_high(n,:));
    else
        x_low(n,:) = maxVal_low*normalise1var(x_low(n,:))   + maxVal_low*0.1*randn(size(x_low(n,:)));
        x_high(n,:)= maxVal_high*normalise1var(x_high(n,:)) + maxVal_low*0.1*randn(size(x_high(n,:)));
    end
    x_lowDec(n,:)  = circshift(x_low(n,:)',-centre_low(n)+500);
    x_highDec(n,:) = circshift(x_high(n,:)',-centre_high(n)+500);
end
x_low(x_low<0) = 0;
x_high(x_high<0) = 0;
subplot(3,4,1)
hold off
plot(inputs( 1:10:end), x_low(1, 1:10:end)','o','color','b', 'linewidth',1);
hold on;
subplot(3,4,5)
plot(inputs(13:10:end), x_low(2,13:10:end)','o','color','r', 'linewidth',1);
% stem(inputs(27:40:end), x_low(3,27:40:end)','o','color','k', 'linewidth',1);
% stem(inputs(31:40:end), x_low(3,31:40:end)','o','color','g', 'linewidth',1);
ylabel('High certainty trials')
% xlabel('Preferred position - Actual position')

subplot(3,4,2)
% hold off;
hold off
plot(inputs( 1:10:end), x_high(1, 1:10:end)','o','color','b', 'linewidth',1);
subplot(3,4,6)
hold on;
plot(inputs(13:10:end), x_high(2,13:10:end)','o','color','r', 'linewidth',1);
% stem(inputs(27:40:end), x_high(3,27:40:end)','o','color','k', 'linewidth',1);
% stem(inputs(31:40:end), x_high(3,31:40:end)','o','color','g', 'linewidth',1);
ylabel('High certainty trials')
ylabel('Low certainty trials')
hold on;
if flag_norm
    subplot(3,4,9)
    plot(inputs, normalise1var(mean(x_low,1)'),'color',[.5 0 .5], 'linewidth',3.5);
    subplot(3,4,10)
    plot(inputs, normalise1var(mean(x_high,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    subplot(3,4,9)
    plot(inputs, mean(x_low,1)','color',[.5 0 .5], 'linewidth',3.5);
    subplot(3,4,10)
    plot(inputs, mean(x_high,1),'color',[0 .5 .5], 'linewidth',3.5);
end
% ylabel('High uncertainty trials')
xlabel('Preferred position - Actual position')

% subplot(2,3,4)
% hold off;
% plot(inputs, x_highDec(1:3,:)','color','c', 'linewidth',1);
% hold on;
% %ylabel('Average activity')
% ylabel('High uncertainty trials')
% xlabel('Preferred position - Decoded position')
% subplot(2,3,4)
% % hold off;
% plot(inputs, x_lowDec(1:3,:)','color','m', 'linewidth',1);
% hold on;
% if flag_norm
%     plot(inputs, normalise1var(mean(x_highDec,1)),'color',[.5 0 .5], 'linewidth',3.5);
%     plot(inputs, normalise1var(mean(x_lowDec,1)),'color',[0 .5 .5], 'linewidth',3.5);
% else
%     plot(inputs, mean(x_highDec,1),'color',[0 .5 .5], 'linewidth',3.5);
%     plot(inputs, mean(x_lowDec,1),'color',[.5 0 .5], 'linewidth',3.5);
% end
% ylabel('Low uncertainty trials')
% %ylabel('Average activity')
% 
%% High spread
x_low = zeros(750,1000);
x_high = zeros(750,1000);
for n = 1:size(x_low,1)
    centre_low(n) = 500+round(20*randn(1));
    centre_high(n) = 500+round(20*randn(1));
    x_low(n,centre_low(n)) = 1;
    x_high(n,centre_high(n)) = 1;
    
    x_low(n,:) = (smthInTime(x_low(n,:),1000,88+15*randn(1)));
    x_high(n,:) = (smthInTime(x_high(n,:),1000,166+15*randn(1)));
end

maxVal_low = mean(max(x_low,[],2));
maxVal_high = mean(max(x_high,[],2));
for n = 1:size(x_low,1)
    if flag_norm
        x_low(n,:) = normalise1var(x_low(n,:));
        x_high(n,:) = normalise1var(x_high(n,:));
    else
        x_low(n,:) = maxVal_low*normalise1var(x_low(n,:))   + maxVal_low*0.1*randn(size(x_low(n,:)));
        x_high(n,:) = maxVal_high*normalise1var(x_high(n,:))+ maxVal_low*0.1*randn(size(x_high(n,:)));
    end
    x_lowDec(n,:)  = circshift(x_low(n,:)',-centre_low(n)+500);
    x_highDec(n,:) = circshift(x_high(n,:)',-centre_high(n)+500);
end
x_low(x_low<0) = 0;
x_high(x_high<0) = 0;

k = 1;
% subplot(3,4,1+k)
% hold off
% plot(inputs, x_low(1:3,:)','color','m', 'linewidth',1);
% hold on;
% ylabel('Low uncertainty trials')
% xlabel('Preferred position - Actual position')

subplot(3,4,3)
hold off;
plot(inputs( 1:10:end), x_high(1, 1:10:end)','o','color','b', 'linewidth',1);
hold on;
subplot(3,4,7)
plot(inputs(13:10:end), x_high(2,13:10:end)','o','color','r', 'linewidth',1);
% stem(inputs(27:40:end), x_high(3,27:40:end)','o','color','k', 'linewidth',1);
% stem(inputs(31:40:end), x_high(3,31:40:end)','o','color','g', 'linewidth',1);
if flag_norm
    subplot(3,4,11)
    plot(inputs, normalise1var(mean(x_high,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    subplot(3,4,11)
    plot(inputs, mean(x_high,1),'color',[0 .5 .5], 'linewidth',3.5);
end
ylabel('Low certainty trials')
xlabel('Preferred position - Actual position')

% subplot(2,3,4+k)
% hold off;
% plot(inputs, x_highDec(1:3,:)','color','c', 'linewidth',1);
% hold on;
% %ylabel('Average activity')
% ylabel('High uncertainty trials')
% xlabel('Preferred position - Decoded position')
% subplot(2,3,4+k)
% % hold off;
% plot(inputs, x_lowDec(1:3,:)','color','m', 'linewidth',1);
% hold on;
% if flag_norm
%     plot(inputs, normalise1var(mean(x_highDec,1)),'color',[.5 0 .5], 'linewidth',3.5);
%     plot(inputs, normalise1var(mean(x_lowDec,1)),'color',[0 .5 .5], 'linewidth',3.5);
% else
%     plot(inputs, mean(x_highDec,1),'color',[0 .5 .5], 'linewidth',3.5);
%     plot(inputs, mean(x_lowDec,1),'color',[.5 0 .5], 'linewidth',3.5);
% end
% ylabel('Low uncertainty trials')
% %ylabel('Average activity')
%% High variability
x_low = zeros(750,1000);
x_high = zeros(750,1000);
rng(2)
for n = 1:size(x_low,1)
    centre_low(n) = 500+round(20*randn(1));
    centre_high(n) = 500+round(120*randn(1));
    x_low(n,centre_low(n)) = 1;
    x_high(n,centre_high(n)) = 1;
    
    x_low(n,:) = (smthInTime(x_low(n,:),1000,88+15*randn(1)));
    x_high(n,:) = (smthInTime(x_high(n,:),1000,88+15*randn(1)));
end

maxVal_low = mean(max(x_low,[],2));
maxVal_high = mean(max(x_high,[],2));
for n = 1:size(x_low,1)
    if flag_norm
        x_low(n,:) = normalise1var(x_low(n,:));
        x_high(n,:) = normalise1var(x_high(n,:));
    else
        x_low(n,:) = maxVal_low*normalise1var(x_low(n,:))+ maxVal_low*0.1*randn(size(x_low(n,:)));
        x_high(n,:) = maxVal_high*normalise1var(x_high(n,:))+ maxVal_low*0.1*randn(size(x_high(n,:)));
    end
    x_lowDec(n,:)  = circshift(x_low(n,:)',-centre_low(n)+500);
    x_highDec(n,:) = circshift(x_high(n,:)',-centre_high(n)+500);
end
x_low(x_low<0) = 0;
x_high(x_high<0) = 0;

% k = 2;
% subplot(2,3,1+k)
% hold off
% plot(inputs, x_low(1:3,:)','color','m', 'linewidth',1);
% hold on;
% ylabel('Low uncertainty trials')
% xlabel('Preferred position - Actual position')

subplot(3,4,4)
hold off
plot(inputs( 1:10:end), x_high(6, 1:10:end)','o','color','b', 'linewidth',1);
hold on;
subplot(3,4,8)
plot(inputs(13:10:end), x_high(14,13:10:end)','o','color','r', 'linewidth',1);
% stem(inputs(27:40:end), x_high(13,27:40:end)','o','color','k', 'linewidth',1);
% stem(inputs(31:40:end), x_high(3,31:40:end)','o','color','g', 'linewidth',1);
% plot(inputs, x_high(6, :)','color','b', 'linewidth',1);
% hold on;
% plot(inputs, x_high(4, :)','color','r', 'linewidth',1);
% plot(inputs, x_high(14, :)','color','k', 'linewidth',1);
if flag_norm
    subplot(3,4,12)
    plot(inputs, normalise1var(mean(x_high,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    subplot(3,4,12)
    plot(inputs, mean(x_high,1),'color',[0 .5 .5], 'linewidth',3.5);
end
ylabel('Low certainty trials')
xlabel('Preferred position - Actual position')

% subplot(2,3,4+k)
% hold off;
% plot(inputs, x_highDec(1:3,:)','color','c', 'linewidth',1);
% hold on;
% %ylabel('Average activity')
% ylabel('High uncertainty trials')
% xlabel('Preferred position - Decoded position')
% subplot(2,3,4+k)
% % hold off;
% plot(inputs, x_lowDec(1:3,:)','color','m', 'linewidth',1);
% hold on;
% if flag_norm
%     plot(inputs, normalise1var(mean(x_highDec,1)),'color',[.5 0 .5], 'linewidth',3.5);
%     plot(inputs, normalise1var(mean(x_lowDec,1)),'color',[0 .5 .5], 'linewidth',3.5);
% else
%     plot(inputs, mean(x_highDec,1),'color',[0 .5 .5], 'linewidth',3.5);
%     plot(inputs, mean(x_lowDec,1),'color',[.5 0 .5], 'linewidth',3.5);
% end
% ylabel('Low uncertainty trials')
% %ylabel('Average activity')

for n = 1:12
    ax(n) = subplot(3,4,n);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    set(gca, 'YTick', [], 'XTick', [-50 -25 0 25 50]...
        )
    axis tight;
    maxY(n) = max(ylim);%,'XTickLabel',[{'-50'}, {'-25'}, {'0'}, {'25'}, {'50'}])
end
linkaxes(ax,'xy')
subplot(3,4,1)
axis tight;
ylims = ylim;
axis([-50 50 0 max(maxY)]);
title('High certainty')
ylabel('Time 1')
for n = 1:12
    subplot(3,4,n);
    line([0 0], ylim, 'linestyle', '--', 'color','k')
end
subplot(3,4,2)
title('Lower firing rate')
subplot(3,4,3)
title('Independent errors')
subplot(3,4,4)
title('Common errors')
subplot(3,4,5)
ylabel('Time 2')
subplot(3,4,9)
ylabel('Averages')