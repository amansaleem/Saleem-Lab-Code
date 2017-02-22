% load Posterior_all_150226_trainedCV_125_noInter.mat
% load Posterior_all_150327_trainedCV_250_noInter
% CA = Posterior_all;
% CA = CA([1:6 8]);
% %
% load Posterior_all_V1_150601_trainedCV_250.mat
% % load Posterior_all_151014_V1_trainedCV_125
% V1 = Posterior_all;

range = 10;
% for delay = 0:-1:-5
delay = 0;
numPoints.cc = zeros(7,3); numPoints.cw = zeros(7,3);
numPoints.wc = zeros(7,3); numPoints.ww = zeros(7,3);

numPoints.ccc = zeros(7,3); numPoints.ccw = zeros(7,3);
numPoints.cwc = zeros(7,3); numPoints.cww = zeros(7,3);

for iseries = 1:7
    C_all.ls = [];    C_all.ln = [];    V_all.ls = [];    V_all.ln = [];
    C_all.ns = [];    C_all.nn = [];    V_all.ns = [];    V_all.nn = [];
    C_all.hs = [];    C_all.hn = [];    V_all.hs = [];    V_all.hn = [];
    C_all.la = [];    C_all.ha = [];    V_all.la = [];    V_all.ha = [];
    C_all.na = [];                      V_all.na = [];
    
    outcome_l = CA(iseries).data.outcome(CA(iseries).t_low);
    outcome_n = CA(iseries).data.outcome(CA(iseries).t_norm);
    outcome_h = CA(iseries).data.outcome(CA(iseries).t_high);
    
    for mid = 1:50
        start = mid-range;
        stop  = mid+range;
        X_ML_l = CA(iseries).X_low_orig;
        V_ML_l = circshift(V1(iseries).MAP.low(1:length(X_ML_l))',delay);
        C_ML_l = CA(iseries).MAP.low(1:length(X_ML_l))';
        O_ML_l = outcome_l(1:length(X_ML_l));
        
        %times when the animal was within a certain range
        temp = X_ML_l==mid;% (X_ML_l>start & X_ML_l<stop);
        % _s: times when either V1 or CA1 was correct
        % _n: times when both V1 and CA1 were wrong
        temp_s = temp & ((C_ML_l>start & C_ML_l<stop) | (V_ML_l>start & V_ML_l<stop));%(V_ML_l>start & V_ML_l<stop) & ~(C_ML_l>start & C_ML_l<stop);%
        temp_n = temp & ~((C_ML_l>start & C_ML_l<stop) | (V_ML_l>start & V_ML_l<stop));%~(V_ML_l>start & V_ML_l<stop) & ((C_ML_l>start & C_ML_l<stop));%
        
        errors_n.com(iseries,1) = sum(temp_n & abs(C_ML_l-V_ML_l)<range)./sum(temp_n);
        errors_n.ind(iseries,1) = sum(temp_n & abs(C_ML_l-V_ML_l)>=range)./sum(temp_n);
        errors_a.com(iseries,1) = sum(temp & abs(C_ML_l-V_ML_l)<range)./sum(temp);
        errors_a.ind(iseries,1) = sum(temp & abs(C_ML_l-V_ML_l)>=range)./sum(temp);
        
        numPoints.cc(iseries,1) = numPoints.cc(iseries,1) + ...
            sum(temp & (V_ML_l>start & V_ML_l<stop) & (C_ML_l>start & C_ML_l<stop));
        numPoints.cw(iseries,1) = numPoints.cw(iseries,1) + ...
            sum(temp & (V_ML_l>start & V_ML_l<stop) & ~(C_ML_l>start & C_ML_l<stop));
        numPoints.wc(iseries,1) = numPoints.wc(iseries,1) + ...
            sum(temp & ~(V_ML_l>start & V_ML_l<stop) & (C_ML_l>start & C_ML_l<stop));
        numPoints.ww(iseries,1) = numPoints.ww(iseries,1) + ...
            sum(temp & ~(V_ML_l>start & V_ML_l<stop) & ~(C_ML_l>start & C_ML_l<stop));
        
        numPoints.ccc(iseries,1) = numPoints.ccc(iseries,1) + ...
            sum(temp & O_ML_l==2 & (V_ML_l>start & V_ML_l<stop) & (C_ML_l>start & C_ML_l<stop));
        numPoints.ccw(iseries,1) = numPoints.ccw(iseries,1) + ...
            sum(temp & O_ML_l==2 & (V_ML_l>start & V_ML_l<stop) & ~(C_ML_l>start & C_ML_l<stop));
        numPoints.cwc(iseries,1) = numPoints.cwc(iseries,1) + ...
            sum(temp & O_ML_l==2 & ~(V_ML_l>start & V_ML_l<stop) & (C_ML_l>start & C_ML_l<stop));
        numPoints.cww(iseries,1) = numPoints.cww(iseries,1) + ...
            sum(temp & O_ML_l==2 & ~(V_ML_l>start & V_ML_l<stop) & ~(C_ML_l>start & C_ML_l<stop));
        
        V_all.la = [V_all.la V_ML_l(temp)'];
        V_all.ls = [V_all.ls V_ML_l(temp_s)'];
        V_all.ln = [V_all.ln V_ML_l(temp_n)'];
        C_all.la = [C_all.la C_ML_l(temp)'];
        C_all.ls = [C_all.ls C_ML_l(temp_s)'];
        C_all.ln = [C_all.ln C_ML_l(temp_n)'];
        
        X_ML_n = CA(iseries).X_norm;
        V_ML_n = circshift(V1(iseries).MAP.norm(1:length(X_ML_n))',delay);
        C_ML_n = CA(iseries).MAP.norm(1:length(X_ML_n))';
        O_ML_n = outcome_n(1:length(X_ML_n));
        
        temp = X_ML_n==mid;%  = (X_ML_n>start & X_ML_n<stop);
        temp_s = temp & ((C_ML_n>start & C_ML_n<stop) | (V_ML_n>start & V_ML_n<stop));%(V_ML_n>start & V_ML_n<stop) & ~(C_ML_n>start & C_ML_n<stop);%
        temp_n = temp & ~((C_ML_n>start & C_ML_n<stop) | (V_ML_n>start & V_ML_n<stop));%~(V_ML_n>start & V_ML_n<stop) & (C_ML_n>start & C_ML_n<stop);%
        
        errors_n.com(iseries,2) = sum(temp_n & abs(C_ML_n-V_ML_n)<range)./sum(temp_n);
        errors_n.ind(iseries,2) = sum(temp_n & abs(C_ML_n-V_ML_n)>=range)./sum(temp_n);
        errors_a.com(iseries,2) = sum(temp & abs(C_ML_n-V_ML_n)<range)./sum(temp);
        errors_a.ind(iseries,2) = sum(temp & abs(C_ML_n-V_ML_n)>=range)./sum(temp);
        
        numPoints.cc(iseries,2) = numPoints.cc(iseries,2) + ...
            sum(temp & (V_ML_n>start & V_ML_n<stop) & (C_ML_n>start & C_ML_n<stop));
        numPoints.cw(iseries,2) = numPoints.cw(iseries,2) + ...
            sum(temp & (V_ML_n>start & V_ML_n<stop) & ~(C_ML_n>start & C_ML_n<stop));
        numPoints.wc(iseries,2) = numPoints.wc(iseries,2) + ...
            sum(temp & ~(V_ML_n>start & V_ML_n<stop) & (C_ML_n>start & C_ML_n<stop));
        numPoints.ww(iseries,2) = numPoints.ww(iseries,2) + ...
            sum(temp & ~(V_ML_n>start & V_ML_n<stop) & ~(C_ML_n>start & C_ML_n<stop));
        
        numPoints.ccc(iseries,2) = numPoints.ccc(iseries,2) + ...
            sum(temp & O_ML_n==2 & (V_ML_n>start & V_ML_n<stop) & (C_ML_n>start & C_ML_n<stop));
        numPoints.ccw(iseries,2) = numPoints.ccw(iseries,2) + ...
            sum(temp & O_ML_n==2 & (V_ML_n>start & V_ML_n<stop) & ~(C_ML_n>start & C_ML_n<stop));
        numPoints.cwc(iseries,2) = numPoints.cwc(iseries,2) + ...
            sum(temp & O_ML_n==2 & ~(V_ML_n>start & V_ML_n<stop) & (C_ML_n>start & C_ML_n<stop));
        numPoints.cww(iseries,2) = numPoints.cww(iseries,2) + ...
            sum(temp & O_ML_n==2 & ~(V_ML_n>start & V_ML_n<stop) & ~(C_ML_n>start & C_ML_n<stop));
        
        V_all.na = [V_all.na V_ML_n(temp)'];
        V_all.ns = [V_all.ns V_ML_n(temp_s)'];
        V_all.nn = [V_all.nn V_ML_n(temp_n)'];
        C_all.na = [C_all.na C_ML_n(temp)'];
        C_all.ns = [C_all.ns C_ML_n(temp_s)'];
        C_all.nn = [C_all.nn C_ML_n(temp_n)'];
        
        X_ML_h = CA(iseries).X_high_orig;
        V_ML_h = circshift(V1(iseries).MAP.high(1:length(X_ML_h))',delay);
        C_ML_h = CA(iseries).MAP.high(1:length(X_ML_h))';
        O_ML_h = outcome_h(1:length(X_ML_h));
        
        temp = X_ML_h==mid;%  = (X_ML_h>start & X_ML_h<stop);
        temp_s = temp & ((C_ML_h>start & C_ML_h<stop) | (V_ML_h>start & V_ML_h<stop));%(V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop);%
        temp_n = temp & ~(V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop);%~((C_ML_h>start & C_ML_h<stop) & (V_ML_h>start & V_ML_h<stop));%
        
        errors_n.com(iseries,3) = sum(temp_n & abs(C_ML_h-V_ML_h)<range)./sum(temp_n);
        errors_n.ind(iseries,3) = sum(temp_n & abs(C_ML_h-V_ML_h)>=range)./sum(temp_n);
        errors_a.com(iseries,3) = sum(temp & abs(C_ML_h-V_ML_h)<range)./sum(temp);
        errors_a.ind(iseries,3) = sum(temp & abs(C_ML_h-V_ML_h)>=range)./sum(temp);
        
        numPoints.cc(iseries,3) = numPoints.cc(iseries,3) + ...
            sum(temp & (V_ML_h>start & V_ML_h<stop) & (C_ML_h>start & C_ML_h<stop));
        numPoints.cw(iseries,3) = numPoints.cw(iseries,3) + ...
            sum(temp & (V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop));
        numPoints.wc(iseries,3) = numPoints.wc(iseries,3) + ...
            sum(temp & ~(V_ML_h>start & V_ML_h<stop) & (C_ML_h>start & C_ML_h<stop));
        numPoints.ww(iseries,3) = numPoints.ww(iseries,3) + ...
            sum(temp & ~(V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop));
        
        numPoints.ccc(iseries,3) = numPoints.ccc(iseries,3) + ...
            sum(temp & O_ML_h==2 & (V_ML_h>start & V_ML_h<stop) & (C_ML_h>start & C_ML_h<stop));
        numPoints.ccw(iseries,3) = numPoints.ccw(iseries,3) + ...
            sum(temp & O_ML_h==2 & (V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop));
        numPoints.cwc(iseries,3) = numPoints.cwc(iseries,3) + ...
            sum(temp & O_ML_h==2 & ~(V_ML_h>start & V_ML_h<stop) & (C_ML_h>start & C_ML_h<stop));
        numPoints.cww(iseries,3) = numPoints.cww(iseries,3) + ...
            sum(temp & O_ML_h==2 & ~(V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop));
        
        V_all.ha = [V_all.ha V_ML_h(temp)'];
        V_all.hs = [V_all.hs V_ML_h(temp_s)'];
        V_all.hn = [V_all.hn V_ML_h(temp_n)'];
        C_all.ha = [C_all.ha C_ML_h(temp)'];
        C_all.hs = [C_all.hs C_ML_h(temp_s)'];
        C_all.hn = [C_all.hn C_ML_h(temp_n)'];
    end
    
    %     for n = 1:6; subplot(2,3,n);
    %         axis equal; axis square; axis([0 50 0 50]); line(xlim,ylim);
    %         hold off;
    %     end
    %     subplot(3,3,2)
    %     hold off
    %     subplot(3,3,3)
    %     hold off
    
    figure(20)
    subplot(3,3,1)
    F = smoothhist2D_corrected([V_all.ls' C_all.ls'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.ls)/(length(V_all.ls)+length(V_all.ln)));
    H.ls(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,4)
    F = smoothhist2D_corrected([V_all.ln' C_all.ln'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.ln)/(length(V_all.ls)+length(V_all.ln)));
    H.ln(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,7)
    F = smoothhist2D_corrected([V_all.la' C_all.la'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.la)/(length(V_all.la)+length(V_all.la)));
    H.la(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,2)
    F = smoothhist2D_corrected([V_all.ns' C_all.ns'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.ns)/(length(V_all.ns)+length(V_all.nn)));
    H.ns(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,5)
    F = smoothhist2D_corrected([V_all.nn' C_all.nn'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.nn)/(length(V_all.ns)+length(V_all.nn)));
    H.nn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,8)
    F = smoothhist2D_corrected([V_all.na' C_all.na'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.na)/(length(V_all.na)+length(V_all.na)));
    H.na(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,3)
    F = smoothhist2D_corrected([V_all.hs' C_all.hs'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.hs)/(length(V_all.hs)+length(V_all.hn)));
    H.hs(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,6)
    F = smoothhist2D_corrected([V_all.hn' C_all.hn'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.hn)/(length(V_all.hs)+length(V_all.hn)));
    H.hn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,9)
    F = smoothhist2D_corrected([V_all.ha' C_all.ha'],4,[50 50],1:50);
%     F = smoothhist2D([V_all.ha' C_all.ha'],4,[50 50],1:50);
    imagesc(1:2:100,1:2:100,F*2500 * length(V_all.ha)/(length(V_all.ha)+length(V_all.ha)));
    H.ha(iseries,:,:) = F; %./ mean(F(:));
    
    for n = 1:9
        subplot(3,3,n);
        axis xy; axis equal; axis square; axis([1+range 50-range 1+range 50-range]); line(xlim,ylim);
        colorbar
    end
    %     drawnow
    %     pause
end
%%
temp = numPoints;
temp.ww = temp.ww ./(numPoints.cc + numPoints.cw + numPoints.wc + numPoints.ww);
temp.cc = temp.cc ./(numPoints.cc + numPoints.cw + numPoints.wc + numPoints.ww);
temp.cw = temp.cw ./(numPoints.cc + numPoints.cw + numPoints.wc + numPoints.ww);
temp.wc = temp.wc ./(numPoints.cc + numPoints.cw + numPoints.wc + numPoints.ww);

temp.cww = temp.cww ./(numPoints.ccc + numPoints.ccw + numPoints.cwc + numPoints.cww);
temp.ccc = temp.ccc ./(numPoints.ccc + numPoints.ccw + numPoints.cwc + numPoints.cww);
temp.ccw = temp.ccw ./(numPoints.ccc + numPoints.ccw + numPoints.cwc + numPoints.cww);
temp.cwc = temp.cwc ./(numPoints.ccc + numPoints.ccw + numPoints.cwc + numPoints.cww);

numPoints = temp;
% for n=1:9; subplot(3,3,n); hold off; end
%%
figure(-delay+1)
clear marginals
for iseries = 1:7
    [marginals.la(iseries,:), X] = get45Marginal(sq(H.la(iseries,:,:)),50);
    [marginals.ls(iseries,:), X] = get45Marginal(sq(H.ls(iseries,:,:)),50);
    [marginals.ln(iseries,:), X] = get45Marginal(sq(H.ln(iseries,:,:)),50);
    
    [marginals.na(iseries,:), X] = get45Marginal(sq(H.na(iseries,:,:)),50);
    [marginals.ns(iseries,:), X] = get45Marginal(sq(H.ns(iseries,:,:)),50);
    [marginals.nn(iseries,:), X] = get45Marginal(sq(H.nn(iseries,:,:)),50);
    
    [marginals.ha(iseries,:), X] = get45Marginal(sq(H.ha(iseries,:,:)),50);
    [marginals.hs(iseries,:), X] = get45Marginal(sq(H.hs(iseries,:,:)),50);
    [marginals.hn(iseries,:), X] = get45Marginal(sq(H.hn(iseries,:,:)),50);
    
%     subplot(3,3,1); imagesc(1:2:100,1:2:100,sq(median(H.ls(iseries,:,:),1)))
%     subplot(3,3,2); imagesc(1:2:100,1:2:100,sq(median(H.ns(iseries,:,:),1)))
%     subplot(3,3,3); imagesc(1:2:100,1:2:100,sq(median(H.hs(iseries,:,:),1)))
%     subplot(3,3,4); imagesc(1:2:100,1:2:100,sq(median(H.ln(iseries,:,:),1)))
%     subplot(3,3,5); imagesc(1:2:100,1:2:100,sq(median(H.nn(iseries,:,:),1)))
%     subplot(3,3,6); imagesc(1:2:100,1:2:100,sq(median(H.hn(iseries,:,:),1)))
%     subplot(3,3,7); imagesc(1:2:100,1:2:100,sq(median(H.la(iseries,:,:),1)))
%     subplot(3,3,8); imagesc(1:2:100,1:2:100,sq(median(H.na(iseries,:,:),1)))
%     subplot(3,3,9); imagesc(1:2:100,1:2:100,sq(median(H.ha(iseries,:,:),1)))
%     for n = 1:9
%         subplot(3,3,n);
%         axis xy; axis equal; axis square; axis([1 50 1 50]); line(xlim,ylim);
%         colorbar;
%     end
%     pause
end
subplot(3,4,1); imagesc(1:2:100,1:2:100,sq(mean(H.la(:,:,:),1)))
title('All (low)')
subplot(3,4,2); imagesc(1:2:100,1:2:100,sq(mean(H.na(:,:,:),1)))
title('All (med)')
subplot(3,4,3); imagesc(1:2:100,1:2:100,sq(mean(H.ha(:,:,:),1)))
title('All (high)')
subplot(3,4,5); imagesc(1:2:100,1:2:100,sq(mean(H.ls(:,:,:),1)))
title('Either correct (low)')
subplot(3,4,6); imagesc(1:2:100,1:2:100,sq(mean(H.ns(:,:,:),1)))
title('Either correct (med)')
subplot(3,4,7); imagesc(1:2:100,1:2:100,sq(mean(H.hs(:,:,:),1)))
title('Either correct (high)')
subplot(3,4,9); imagesc(1:2:100,1:2:100,sq(mean(H.ln(:,:,:),1)))
title('Both wrong (low)')
subplot(3,4,10); imagesc(1:2:100,1:2:100,sq(mean(H.nn(:,:,:),1)))
title('Both wrong (med)')
subplot(3,4,11); imagesc(1:2:100,1:2:100,sq(mean(H.hn(:,:,:),1)))
title('Both wrong (high)')
X = X*2*sqrt(2);
t = X>-90 & X<90;

alpha_val = 1;
subplot(3,4,4)
hold off
meanMar = mean(marginals.la(:,t),1);
semMar = sem(marginals.la(:,t));
errorarea_as(X(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = mean(marginals.na(:,t),1);
semMar = sem(marginals.na(:,t));
errorarea_as(X(t), meanMar,semMar,'k',alpha_val);
meanMar = mean(marginals.ha(:,t),1);
semMar = sem(marginals.ha(:,t));
errorarea_as(X(t), meanMar,semMar,'r',alpha_val);

subplot(3,4,8)
hold off
meanMar = mean(marginals.ls(:,t),1);
semMar = sem(marginals.ls(:,t));
errorarea_as(X(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = mean(marginals.ns(:,t),1);
semMar = sem(marginals.ns(:,t));
errorarea_as(X(t), meanMar,semMar,'k',alpha_val);
meanMar = mean(marginals.hs(:,t),1);
semMar = sem(marginals.hs(:,t));
errorarea_as(X(t), meanMar,semMar,'r',alpha_val);

subplot(3,4,12)
hold off
meanMar = mean(marginals.ln(:,t),1);
semMar = sem(marginals.ln(:,t));
errorarea_as(X(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = mean(marginals.nn(:,t),1);
semMar = sem(marginals.nn(:,t));
errorarea_as(X(t), meanMar,semMar,'k',alpha_val);
meanMar = mean(marginals.hn(:,t),1);
semMar = sem(marginals.hn(:,t));
errorarea_as(X(t), meanMar,semMar,'r',alpha_val);

for n = [1:3 5:7 9:11]
    subplot(3,4,n);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    axis xy; axis equal; axis square; axis([1 100 1 100]); line(xlim,ylim);
    set(gca,'clim',[0.5 2])
    colorbar;
    hold off
end


for n = 4:4:12
    subplot(3,4,n);
    axis tight;
    line(xlim, [0 0], 'color','k','linestyle','--')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
end

%%
figure(5)
subplot(221)
hold off;
errorbar([1 2 3], mean(temp.wc), sem(temp.wc),'k','linewidth',2)
hold on;
errorbar([1 2 3], mean(temp.cw), sem(temp.cw),'k','linewidth',2)
errorbar([1 2 3], mean(temp.cc), sem(temp.cc),'g','linewidth',2)
errorbar([1 2 3], mean(temp.ww), sem(temp.ww),'r','linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis([0.5 3.5 0 0.5])
ylabel('Fraction of time')

subplot(222)
hold off;
errorbar([1 2 3], mean(temp.wc./(temp.wc+temp.ww)),...
    sem(temp.wc./(temp.wc+temp.ww)),'r','linewidth',2)
hold on;
errorbar([1 2 3], mean(temp.cc./(temp.cw+temp.cc)),...
    sem(temp.cc./(temp.cw+temp.cc)),'g','linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis([0.5 3.5 0 0.8])
ylabel('Fraction of time CA1 right')
title('V1 right / wrong')

subplot(223)
hold off;
errorbar([1 2 3], mean(temp.cw./(temp.cw+temp.ww)),...
    sem(temp.cw./(temp.cw+temp.ww)),'r','linewidth',2)
hold on;
errorbar([1 2 3], mean(temp.cc./(temp.wc+temp.cc)),...
    sem(temp.cc./(temp.wc+temp.cc)),'g','linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis([0.5 3.5 0 0.8])
ylabel('Fraction of time V1 right')
title('CA1 right / wrong')

subplot(224)
hold off;
errorbar([1 2 3], mean(errors_a.com), sem(errors_a.com),'k','linewidth',2)
hold on;
errorbar([1 2 3], mean(errors_n.com), sem(errors_n.com),'m','linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis([0.5 3.5 0 0.85])
ylabel('Fraction of common decoded positions')
title('All and both wrong');

%%
%%
figure(6)
subplot(221)
hold off;
errorbar([1 2 3], mean(temp.cwc), sem(temp.cwc),'k','linewidth',2)
hold on;
errorbar([1 2 3], mean(temp.ccw), sem(temp.ccw),'k','linewidth',2)
errorbar([1 2 3], mean(temp.ccc), sem(temp.ccc),'g','linewidth',2)
errorbar([1 2 3], mean(temp.cww), sem(temp.cww),'r','linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis([0.5 3.5 0 0.5])
ylabel('Fraction of time')

subplot(222)
hold off;
errorbar([1 2 3], mean(temp.cwc./(temp.cwc+temp.cww)),...
    sem(temp.cwc./(temp.cwc+temp.cww)),'r','linewidth',2)
hold on;
errorbar([1 2 3], mean(temp.ccc./(temp.ccw+temp.ccc)),...
    sem(temp.ccc./(temp.ccw+temp.ccc)),'g','linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis([0.5 3.5 0 0.8])
ylabel('Fraction of time CA1 right')
title('V1 right / wrong')

subplot(223)
hold off;
errorbar([1 2 3], mean(temp.ccw./(temp.ccw+temp.cww)),...
    sem(temp.ccw./(temp.ccw+temp.cww)),'r','linewidth',2)
hold on;
errorbar([1 2 3], mean(temp.ccc./(temp.cwc+temp.ccc)),...
    sem(temp.ccc./(temp.cwc+temp.ccc)),'g','linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis([0.5 3.5 0 0.8])
ylabel('Fraction of time V1 right')
title('CA1 right / wrong')

subplot(224)
hold off;
errorbar([1 2 3], mean(errors_a.com), sem(errors_a.com),'k','linewidth',2)
hold on;
errorbar([1 2 3], mean(errors_n.com), sem(errors_n.com),'m','linewidth',2)
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis([0.5 3.5 0 0.85])
ylabel('Fraction of common decoded positions')
title('All and both wrong');
