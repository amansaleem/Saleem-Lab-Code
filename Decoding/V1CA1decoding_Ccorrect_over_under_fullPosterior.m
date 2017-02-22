% % % load Posterior_all_150226_trainedCV_125_noInter.mat
% % load Posterior_all_150327_trainedCV_250_noInter
% load Posterior_all_160121_trainedCV_250
% % load Posterior_all_160121_trainedCV_125
% CA = Posterior_all;
% % CA = CA([1:6 8]);
% % % % %
% % % load Posterior_all_V1_150601_trainedCV_250.mat
% % % % % load Posterior_all_151014_V1_trainedCV_125
% % % load Posterior_all_160120_V1_trainedCV_125
% load Posterior_all_160121_V1_trainedCV_250
% % load Posterior_all_160121_V1_trainedCV_125
% V1 = Posterior_all;

alpha_val = 0.2;
range = 10;
range2 = 50;
% for delay = 0:-1:-5
delay = 0;
numPoints.cc = zeros(7,3); numPoints.cw = zeros(7,3);
numPoints.wc = zeros(7,3); numPoints.ww = zeros(7,3);

for iseries = 1:7
    C_all.ls = [];    C_all.ln = [];    V_all.ls = [];    V_all.ln = [];
    C_all.ns = [];    C_all.nn = [];    V_all.ns = [];    V_all.nn = [];
    C_all.hs = [];    C_all.hn = [];    V_all.hs = [];    V_all.hn = [];
    C_all.la = [];    C_all.ha = [];    V_all.la = [];    V_all.ha = [];
    C_all.na = [];                      V_all.na = [];
    X_all.ls = [];    X_all.ln = [];    
    X_all.ns = [];    X_all.nn = [];    
    X_all.hs = [];    X_all.hn = [];    
    X_all.la = [];    X_all.ha = [];    
    X_all.na = [];                      
    X_all.as = [];    X_all.an = [];    X_all.aa = [];
    C_all.as = [];    C_all.an = [];    C_all.aa = [];
    V_all.as = [];    V_all.an = [];    V_all.aa = [];
    
    for mid = 1:50
        start = mid-range;
        stop  = mid+range;
        X_ML_l = V1(iseries).X_low_orig;
        V_ML_l = circshift(V1(iseries).MAP.low(1:length(X_ML_l))',delay);
        C_ML_l = CA(iseries).MAP.low(1:length(X_ML_l))';
        
        %times when the animal was within a certain range and CA1 was
        %correct
        temp = X_ML_l==mid;% & (C_ML_l>start & C_ML_l<stop);% (X_ML_l>start & X_ML_l<stop);
        % _s: CA1 underestimated 
        % _n: CA1 overshot
        temp_s = temp & ((C_ML_l<start)) & C_ML_l>start-range2;% | (V_ML_l>start & V_ML_l<stop));%(V_ML_l>start & V_ML_l<stop) & ~(C_ML_l>start & C_ML_l<stop);%
        temp_n = temp & ((C_ML_l>stop)) & C_ML_l<stop+range2;% | (V_ML_l>start & V_ML_l<stop));%~(V_ML_l>start & V_ML_l<stop) & ((C_ML_l>start & C_ML_l<stop));%
        
        V1p(iseries).la(mid,:) = nanmean(V1(iseries).Posterior_low(temp,:));
        V1p(iseries).ls(mid,:) = nanmean(V1(iseries).Posterior_low(temp_s,:));
        V1p(iseries).ln(mid,:) = nanmean(V1(iseries).Posterior_low(temp_n,:));

        CAp(iseries).la(mid,:) = nanmean(CA(iseries).Posterior_low(temp,:));
        CAp(iseries).ls(mid,:) = nanmean(CA(iseries).Posterior_low(temp_s,:));
        CAp(iseries).ln(mid,:) = nanmean(CA(iseries).Posterior_low(temp_n,:));
        
        % Normal
        X_ML_n = V1(iseries).X_norm;
        V_ML_n = circshift(V1(iseries).MAP.norm(1:length(X_ML_n))',delay);
        C_ML_n = CA(iseries).MAP.norm(1:length(X_ML_n))';
        
        temp = X_ML_n==mid;% & (C_ML_n>start & C_ML_n<stop);%  = (X_ML_n>start & X_ML_n<stop);
        temp_s = temp & ((C_ML_n<start))& C_ML_n>start-range2;% | (V_ML_n>start & V_ML_n<stop));%(V_ML_n>start & V_ML_n<stop) & ~(C_ML_n>start & C_ML_n<stop);%
        temp_n = temp & ((C_ML_n>stop))& C_ML_n<stop+range2;% | (V_ML_n>start & V_ML_n<stop));%~(V_ML_n>start & V_ML_n<stop) & (C_ML_n>start & C_ML_n<stop);%
        
        V1p(iseries).na(mid,:) = nanmean(V1(iseries).Posterior_norm(temp,:));
        V1p(iseries).ns(mid,:) = nanmean(V1(iseries).Posterior_norm(temp_s,:));
        V1p(iseries).nn(mid,:) = nanmean(V1(iseries).Posterior_norm(temp_n,:));

        CAp(iseries).na(mid,:) = nanmean(CA(iseries).Posterior_norm(temp,:));
        CAp(iseries).ns(mid,:) = nanmean(CA(iseries).Posterior_norm(temp_s,:));
        CAp(iseries).nn(mid,:) = nanmean(CA(iseries).Posterior_norm(temp_n,:));
        
        % High
        X_ML_h = V1(iseries).X_high_orig;
        V_ML_h = circshift(V1(iseries).MAP.high(1:length(X_ML_h))',delay);
        C_ML_h = CA(iseries).MAP.high(1:length(X_ML_h))';
        
        temp = X_ML_h==mid;% & (C_ML_h>start & C_ML_h<stop);%  = (X_ML_h>start & X_ML_h<stop);
        temp_s = temp & ((C_ML_h<start)) & C_ML_h>start-range2;% | (V_ML_h>start & V_ML_h<stop));%(V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop);%
        temp_n = temp & ((C_ML_h>stop)) & C_ML_h<stop+range2;% & ~(C_ML_h>start & C_ML_h<stop);%~((C_ML_h>start & C_ML_h<stop) & (V_ML_h>start & V_ML_h<stop));%
        
        V1p(iseries).ha(mid,:) = nanmean(V1(iseries).Posterior_high(temp,:));
        V1p(iseries).hs(mid,:) = nanmean(V1(iseries).Posterior_high(temp_s,:));
        V1p(iseries).hn(mid,:) = nanmean(V1(iseries).Posterior_high(temp_n,:));

        CAp(iseries).ha(mid,:) = nanmean(CA(iseries).Posterior_high(temp,:));
        CAp(iseries).hs(mid,:) = nanmean(CA(iseries).Posterior_high(temp_s,:));
        CAp(iseries).hn(mid,:) = nanmean(CA(iseries).Posterior_high(temp_n,:));
    end
end
%%
clear marginals
% range = range*2;

% for iseries = 1:7
%     [marginals.la(iseries,:), X] = get45Marginal(sq(CAp(iseries).la(:,:))',50);
%     [marginals.ls(iseries,:), X] = get45Marginal(sq(H.ls(iseries,:,:))',50);
%     [marginals.ln(iseries,:), X] = get45Marginal(sq(H.ln(iseries,:,:))',50);
%     
%     [marginals.na(iseries,:), X] = get45Marginal(sq(H.na(iseries,:,:))',50);
%     [marginals.ns(iseries,:), X] = get45Marginal(sq(H.ns(iseries,:,:))',50);
%     [marginals.nn(iseries,:), X] = get45Marginal(sq(H.nn(iseries,:,:))',50);
%     
%     [marginals.ha(iseries,:), X] = get45Marginal(sq(H.ha(iseries,:,:))',50);
%     [marginals.hs(iseries,:), X] = get45Marginal(sq(H.hs(iseries,:,:))',50);
%     [marginals.hn(iseries,:), X] = get45Marginal(sq(H.hn(iseries,:,:))',50);
%     
%     [marginals.xla(iseries,:), X] = get45Marginal(sq(H.xla(iseries,:,:))',50);
%     [marginals.xls(iseries,:), X] = get45Marginal(sq(H.xls(iseries,:,:))',50);
%     [marginals.xln(iseries,:), X] = get45Marginal(sq(H.xln(iseries,:,:))',50);
%     
%     [marginals.xna(iseries,:), X] = get45Marginal(sq(H.xna(iseries,:,:))',50);
%     [marginals.xns(iseries,:), X] = get45Marginal(sq(H.xns(iseries,:,:))',50);
%     [marginals.xnn(iseries,:), X] = get45Marginal(sq(H.xnn(iseries,:,:))',50);
%     
%     [marginals.xha(iseries,:), X] = get45Marginal(sq(H.xha(iseries,:,:))',50);
%     [marginals.xhs(iseries,:), X] = get45Marginal(sq(H.xhs(iseries,:,:))',50);
%     [marginals.xhn(iseries,:), X] = get45Marginal(sq(H.xhn(iseries,:,:))',50);
%     
%     [marginals.aa(iseries,:), X] = get45Marginal(sq(H.aa(iseries,:,:))',50);
%     [marginals.as(iseries,:), X] = get45Marginal(sq(H.as(iseries,:,:))',50);
%     [marginals.an(iseries,:), X] = get45Marginal(sq(H.an(iseries,:,:))',50);
%     
%     [marginals.xaa(iseries,:), X] = get45Marginal(sq(H.xaa(iseries,:,:))',50);
%     [marginals.xas(iseries,:), X] = get45Marginal(sq(H.xas(iseries,:,:))',50);
%     [marginals.xan(iseries,:), X] = get45Marginal(sq(H.xan(iseries,:,:))',50);
%     
% end
%% Plotting V and X
figure(-delay+10)
subplot(3,4,1); imagesc(c1, c2,sq(mean(H.la(:,:,:),1))')
title('All (low)')
subplot(3,4,2); imagesc(c1, c2,sq(mean(H.na(:,:,:),1))')
title('All (med)')
subplot(3,4,3); imagesc(c1, c2,sq(mean(H.ha(:,:,:),1))')
title('All (high)')
subplot(3,4,5); imagesc(c1, c2,sq(mean(H.ls(:,:,:),1))')
title('CA1 undershot (low)')
subplot(3,4,6); imagesc(c1, c2,sq(mean(H.ns(:,:,:),1))')
title('CA1 undershot (med)')
subplot(3,4,7); imagesc(c1, c2,sq(mean(H.hs(:,:,:),1))')
title('CA1 undershot (high)')
subplot(3,4,9); imagesc(c1, c2,sq(mean(H.ln(:,:,:),1))')
title('CA1 overshot (low)')
subplot(3,4,10); imagesc(c1, c2,sq(mean(H.nn(:,:,:),1))')
title('CA1 overshot (med)')
subplot(3,4,11); imagesc(c1, c2,sq(mean(H.hn(:,:,:),1))')
title('CA1 overshot (high)')
Xtemp = X*2*sqrt(2);
t = Xtemp~=-1000; % Xtemp>-90 & Xtemp<90;

subplot(3,4,4)
hold off
meanMar = nanmean(marginals.la(:,t),1);
semMar = nansem(marginals.la(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.na(:,t),1);
semMar = nansem(marginals.na(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.ha(:,t),1);
semMar = nansem(marginals.ha(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

subplot(3,4,8)
hold off
meanMar = nanmean(marginals.ls(:,t),1);
semMar = nansem(marginals.ls(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.ns(:,t),1);
semMar = nansem(marginals.ns(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.hs(:,t),1);
semMar = nansem(marginals.hs(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

subplot(3,4,12)
hold off
meanMar = nanmean(marginals.ln(:,t),1);
semMar = nansem(marginals.ln(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.nn(:,t),1);
semMar = nansem(marginals.nn(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.hn(:,t),1);
semMar = nansem(marginals.hn(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

for n = [1:3 5:7 9:11]
    subplot(3,4,n);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    axis xy; axis equal; axis square; axis([1 50 1 50]); line(xlim,ylim);
%     set(gca,'clim',[0.5 1.5])
    colorbar;
    colormap(jet);
    xlabel('Position');
    ylabel('V1 decoded position');
    hold off
end


for n = 4:4:12
    subplot(3,4,n);
    axis tight;
    line(xlim, [0 0], 'color','k','linestyle','--')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line([0 0], ylim, 'color','k','linestyle','--');
    line(xlim, [1 1], 'color','k','linestyle','--');
end
%% Plotting X and C
figure(-delay+11)
subplot(3,4,1); imagesc(c1, c2,sq(mean(H.xla(:,:,:),1))')
title('All (low)')
subplot(3,4,2); imagesc(c1, c2,sq(mean(H.xna(:,:,:),1))')
title('All (med)')
subplot(3,4,3); imagesc(c1, c2,sq(mean(H.xha(:,:,:),1))')
title('All (high)')
subplot(3,4,5); imagesc(c1, c2,sq(mean(H.xls(:,:,:),1))')
title('CA1 undershot (low)')
subplot(3,4,6); imagesc(c1, c2,sq(mean(H.xns(:,:,:),1))')
title('CA1 undershot (med)')
subplot(3,4,7); imagesc(c1, c2,sq(mean(H.xhs(:,:,:),1))')
title('CA1 undershot (high)')
subplot(3,4,9); imagesc(c1, c2,sq(mean(H.xln(:,:,:),1))')
title('CA1 overshot (low)')
subplot(3,4,10); imagesc(c1, c2,sq(mean(H.xnn(:,:,:),1))')
title('CA1 overshot (med)')
subplot(3,4,11); imagesc(c1, c2,sq(mean(H.xhn(:,:,:),1))')
title('CA1 overshot (high)')
Xtemp = X*2*sqrt(2);
t = Xtemp>-1000;%Xtemp>-90 & Xtemp<90;

subplot(3,4,4)
hold off
meanMar = nanmean(marginals.xla(:,t),1);
semMar = nansem(marginals.xla(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.na(:,t),1);
semMar = nansem(marginals.na(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.ha(:,t),1);
semMar = nansem(marginals.ha(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

subplot(3,4,8)
hold off
meanMar = nanmean(marginals.xls(:,t),1);
semMar = nansem(marginals.xls(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.xns(:,t),1);
semMar = nansem(marginals.xns(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.xhs(:,t),1);
semMar = nansem(marginals.xhs(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

subplot(3,4,12)
hold off
meanMar = nanmean(marginals.xln(:,t),1);
semMar = nansem(marginals.xln(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.xnn(:,t),1);
semMar = nansem(marginals.xnn(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.xhn(:,t),1);
semMar = nansem(marginals.xhn(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

for n = [1:3 5:7 9:11]
    subplot(3,4,n);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    axis xy; axis equal; axis square; axis([1 50 1 50]); line(xlim,ylim);
%     set(gca,'clim',[0.5 1.5])
    colorbar;
    colormap(jet);
    xlabel('Position');
    ylabel('CA1 decoded position');
    hold off
end


for n = 4:4:12
    subplot(3,4,n);
    axis tight;
    line(xlim, [0 0], 'color','k','linestyle','--')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line([0 0], ylim, 'color','k','linestyle','--');
    line(xlim, [1 1], 'color','k','linestyle','--');
end
%% Plotting V and X for all contrasts
figure(-delay+13)
subplot(3,2,1); imagesc(c1, c2,sq(mean(H.aa(:,:,:),1))')
title('All (low)')
subplot(3,2,3); imagesc(c1, c2,sq(mean(H.as(:,:,:),1))')
title('CA1 undershot (low)')
subplot(3,2,5); imagesc(c1, c2,sq(mean(H.an(:,:,:),1))')
title('CA1 overshot (low)')
Xtemp = X*2*sqrt(2);
t = Xtemp~=-1000; % Xtemp>-90 & Xtemp<90;

subplot(3,2,2)
hold off
meanMar = nanmean(marginals.la(:,t),1);
semMar = nansem(marginals.la(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.na(:,t),1);
semMar = nansem(marginals.na(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.ha(:,t),1);
semMar = nansem(marginals.ha(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

subplot(3,2,4)
hold off
meanMar = nanmean(marginals.ls(:,t),1);
semMar = nansem(marginals.ls(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.ns(:,t),1);
semMar = nansem(marginals.ns(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.hs(:,t),1);
semMar = nansem(marginals.hs(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

subplot(3,2,6)
hold off
meanMar = nanmean(marginals.ln(:,t),1);
semMar = nansem(marginals.ln(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'b',alpha_val); hold on;
meanMar = nanmean(marginals.nn(:,t),1);
semMar = nansem(marginals.nn(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'k',alpha_val);
meanMar = nanmean(marginals.hn(:,t),1);
semMar = nansem(marginals.hn(:,t));
errorarea_as(Xtemp(t), meanMar,semMar,'r',alpha_val);

for n = [1:2:6]
    subplot(3,2,n);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    axis xy; axis equal; axis square; axis([1 50 1 50]); line(xlim,ylim);
%     set(gca,'clim',[0.5 1.5])
    colorbar;
    colormap(jet);
    xlabel('Position');
    ylabel('V1 decoded position');
    hold off
end


for n = 2:2:6
    subplot(3,2,n);
    axis tight;
    line(xlim, [0 0], 'color','k','linestyle','--')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line([0 0], ylim, 'color','k','linestyle','--');
    line(xlim, [1 1], 'color','k','linestyle','--');
end

%%
subplot(331)
imagesc(V1p(5).la');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)
subplot(332)
imagesc(V1p(5).na');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)
subplot(333)
imagesc(V1p(5).ha');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)
subplot(334)
imagesc(V1p(5).ls');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)
subplot(335)
imagesc(V1p(5).ns');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)
subplot(336)
imagesc(V1p(5).hs');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)
subplot(337)
imagesc(V1p(5).ln');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)
subplot(338)
imagesc(V1p(5).nn');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)
subplot(339)
imagesc(V1p(5).hn');
axis equal; axis xy; RedWhiteBlue; axis tight; set(gca, 'Clim',[-0.5 0.5]);line(xlim, ylim)