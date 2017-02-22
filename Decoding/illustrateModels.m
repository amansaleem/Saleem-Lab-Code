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
    
    x_low(n,:) = (smthInTime(x_low(n,:),1000,66+15*randn(1)));
    x_high(n,:) = (smthInTime(x_high(n,:)/3,1000,66+15*randn(1)));
    
end

maxVal_low = mean(max(x_low,[],2));
maxVal_high = mean(max(x_high,[],2));
for n = 1:size(x_low,1)
    if flag_norm
        x_low(n,:) = normalise1var(x_low(n,:));
        x_high(n,:) = normalise1var(x_high(n,:));
    else
        x_low(n,:) = maxVal_low*normalise1var(x_low(n,:));
        x_high(n,:) = maxVal_low*normalise1var(x_high(n,:));
    end
    x_lowDec(n,:)  = circshift(x_low(n,:)',-centre_low(n)+500);
    x_highDec(n,:) = circshift(x_high(n,:)',-centre_high(n)+500);
end

subplot(2,3,1)
hold off
plot(inputs, x_low(1:6,:)','color','m', 'linewidth',1);
hold on;
ylabel('Low uncertainty trials')
% xlabel('Preferred position - Actual position')

subplot(2,3,1)
% hold off;
plot(inputs, x_high(1:6,:)','color','c', 'linewidth',1);
hold on;
if flag_norm
    plot(inputs, normalise1var(mean(x_low,1)'),'color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, normalise1var(mean(x_high,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    plot(inputs, mean(x_low,1)','color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, mean(x_high,1),'color',[0 .5 .5], 'linewidth',3.5);
end
ylabel('High uncertainty trials')
xlabel('Preferred position - Actual position')

subplot(2,3,4)
hold off;
plot(inputs, x_highDec(1:6,:)','color','c', 'linewidth',1);
hold on;
%ylabel('Average activity')
ylabel('High uncertainty trials')
xlabel('Preferred position - Decoded position')
subplot(2,3,4)
% hold off;
plot(inputs, x_lowDec(1:6,:)','color','m', 'linewidth',1);
hold on;
if flag_norm
    plot(inputs, normalise1var(mean(x_highDec,1)),'color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, normalise1var(mean(x_lowDec,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    plot(inputs, mean(x_highDec,1),'color',[0 .5 .5], 'linewidth',3.5);
    plot(inputs, mean(x_lowDec,1),'color',[.5 0 .5], 'linewidth',3.5);
end
ylabel('Low uncertainty trials')
%ylabel('Average activity')

%% High spread
x_low = zeros(750,1000);
x_high = zeros(750,1000);
for n = 1:size(x_low,1)
    centre_low(n) = 500+round(20*randn(1));
    centre_high(n) = 500+round(20*randn(1));
    x_low(n,centre_low(n)) = 1;
    x_high(n,centre_high(n)) = 1;
    
    x_low(n,:) = (smthInTime(x_low(n,:),1000,66+15*randn(1)));
    x_high(n,:) = (smthInTime(x_high(n,:),1000,132+15*randn(1)));
end

maxVal_low = mean(max(x_low,[],2));
maxVal_high = mean(max(x_high,[],2));
for n = 1:size(x_low,1)
    if flag_norm
        x_low(n,:) = normalise1var(x_low(n,:));
        x_high(n,:) = normalise1var(x_high(n,:));
    else
        x_low(n,:) = maxVal_low*normalise1var(x_low(n,:));
        x_high(n,:) = maxVal_high*normalise1var(x_high(n,:));
    end
    x_lowDec(n,:)  = circshift(x_low(n,:)',-centre_low(n)+500);
    x_highDec(n,:) = circshift(x_high(n,:)',-centre_high(n)+500);
end

k = 1;
subplot(2,3,1+k)
hold off
plot(inputs, x_low(1:6,:)','color','m', 'linewidth',1);
hold on;
ylabel('Low uncertainty trials')
% xlabel('Preferred position - Actual position')

subplot(2,3,1+k)
% hold off;
plot(inputs, x_high(1:6,:)','color','c', 'linewidth',1);
hold on;
if flag_norm
    plot(inputs, normalise1var(mean(x_low,1)'),'color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, normalise1var(mean(x_high,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    plot(inputs, mean(x_low,1)','color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, mean(x_high,1),'color',[0 .5 .5], 'linewidth',3.5);
end
ylabel('High uncertainty trials')
xlabel('Preferred position - Actual position')

subplot(2,3,4+k)
hold off;
plot(inputs, x_highDec(1:6,:)','color','c', 'linewidth',1);
hold on;
%ylabel('Average activity')
ylabel('High uncertainty trials')
xlabel('Preferred position - Decoded position')
subplot(2,3,4+k)
% hold off;
plot(inputs, x_lowDec(1:6,:)','color','m', 'linewidth',1);
hold on;
if flag_norm
    plot(inputs, normalise1var(mean(x_highDec,1)),'color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, normalise1var(mean(x_lowDec,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    plot(inputs, mean(x_highDec,1),'color',[0 .5 .5], 'linewidth',3.5);
    plot(inputs, mean(x_lowDec,1),'color',[.5 0 .5], 'linewidth',3.5);
end
ylabel('Low uncertainty trials')
%ylabel('Average activity')
%% High variability
x_low = zeros(750,1000);
x_high = zeros(750,1000);
for n = 1:size(x_low,1)
    centre_low(n) = 500+round(20*randn(1));
    centre_high(n) = 500+round(120*randn(1));
    x_low(n,centre_low(n)) = 1;
    x_high(n,centre_high(n)) = 1;
    
    x_low(n,:) = (smthInTime(x_low(n,:),1000,66+15*randn(1)));
    x_high(n,:) = (smthInTime(x_high(n,:),1000,66+15*randn(1)));
end

maxVal_low = mean(max(x_low,[],2));
maxVal_high = mean(max(x_high,[],2));
for n = 1:size(x_low,1)
    if flag_norm
        x_low(n,:) = normalise1var(x_low(n,:));
        x_high(n,:) = normalise1var(x_high(n,:));
    else
        x_low(n,:) = maxVal_low*normalise1var(x_low(n,:));
        x_high(n,:) = maxVal_high*normalise1var(x_high(n,:));
    end
    x_lowDec(n,:)  = circshift(x_low(n,:)',-centre_low(n)+500);
    x_highDec(n,:) = circshift(x_high(n,:)',-centre_high(n)+500);
end

k = 2;
subplot(2,3,1+k)
hold off
plot(inputs, x_low(1:6,:)','color','m', 'linewidth',1);
hold on;
ylabel('Low uncertainty trials')
% xlabel('Preferred position - Actual position')

subplot(2,3,1+k)
% hold off;
plot(inputs, x_high(1:6,:)','color','c', 'linewidth',1);
hold on;
if flag_norm
    plot(inputs, normalise1var(mean(x_low,1)'),'color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, normalise1var(mean(x_high,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    plot(inputs, mean(x_low,1)','color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, mean(x_high,1),'color',[0 .5 .5], 'linewidth',3.5);
end
ylabel('High uncertainty trials')
xlabel('Preferred position - Actual position')

subplot(2,3,4+k)
hold off;
plot(inputs, x_highDec(1:6,:)','color','c', 'linewidth',1);
hold on;
%ylabel('Average activity')
ylabel('High uncertainty trials')
xlabel('Preferred position - Decoded position')
subplot(2,3,4+k)
% hold off;
plot(inputs, x_lowDec(1:6,:)','color','m', 'linewidth',1);
hold on;
if flag_norm
    plot(inputs, normalise1var(mean(x_highDec,1)),'color',[.5 0 .5], 'linewidth',3.5);
    plot(inputs, normalise1var(mean(x_lowDec,1)),'color',[0 .5 .5], 'linewidth',3.5);
else
    plot(inputs, mean(x_highDec,1),'color',[0 .5 .5], 'linewidth',3.5);
    plot(inputs, mean(x_lowDec,1),'color',[.5 0 .5], 'linewidth',3.5);
end
ylabel('Low uncertainty trials')
%ylabel('Average activity')

for n = 1:6
    ax(n) = subplot(2,3,n);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    set(gca, 'YTick', [], 'XTick', [-50 -25 0 25 50]...
        )
    axis tight;
    maxY(n) = max(ylim);%,'XTickLabel',[{'-50'}, {'-25'}, {'0'}, {'25'}, {'50'}])
end
linkaxes(ax,'xy')
subplot(2,3,1)
axis tight;
ylims = ylim;
axis([-50 50 0 max(maxY)]);