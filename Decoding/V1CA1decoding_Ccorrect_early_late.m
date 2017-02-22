% % load Posterior_all_150226_trainedCV_125_noInter.mat
% load Posterior_all_150327_trainedCV_250_noInter
% CA = Posterior_all;
% CA = CA([1:6 8]);
% %
% load Posterior_all_V1_150601_trainedCV_250.mat
% % load Posterior_all_151014_V1_trainedCV_125
% V1 = Posterior_all;

alpha_val = 0.2;
range = 15;
range2 = 50;
win = 3;
% for delay = 0:-1:-5
delay = 0;
numPoints.cc = zeros(7,3); numPoints.cw = zeros(7,3);
numPoints.wc = zeros(7,3); numPoints.ww = zeros(7,3);

for iseries = [1:7]
    C_all.ls = [];    C_all.ln = [];    V_all.ls = [];    V_all.ln = [];
    C_all.ns = [];    C_all.nn = [];    V_all.ns = [];    V_all.nn = [];
    C_all.hs = [];    C_all.hn = [];    V_all.hs = [];    V_all.hn = [];
    C_all.la = [];    C_all.ha = [];    V_all.la = [];    V_all.ha = [];
    C_all.na = [];                      V_all.na = [];
    X_all.ls = [];    X_all.ln = [];    V_all.ls = [];    V_all.ln = [];
    X_all.ns = [];    X_all.nn = [];    V_all.ns = [];    V_all.nn = [];
    X_all.hs = [];    X_all.hn = [];    V_all.hs = [];    V_all.hn = [];
    X_all.la = [];    X_all.ha = [];    V_all.la = [];    V_all.ha = [];
    X_all.na = [];                      V_all.na = [];
    outcome_l = CA(iseries).data.outcome(CA(iseries).t_low);
    outcome_n = CA(iseries).data.outcome(CA(iseries).t_norm);
    outcome_h = CA(iseries).data.outcome(CA(iseries).t_high);
    
    for mid = 1:50
        X_ML_l = CA(iseries).X_low_orig;
        V_ML_l = circshift(V1(iseries).MAP.low(1:length(X_ML_l))',delay);
        C_ML_l = CA(iseries).MAP.low(1:length(X_ML_l))';
        O_ML_l = outcome_l(1:length(X_ML_l));
        
        %times when the animal was within a certain range and CA1 was
        %correct
        temp = X_ML_l==mid & O_ML_l==2;% & (C_ML_l>start & C_ML_l<stop);% (X_ML_l>start & X_ML_l<stop);
        % _s: CA1 underestimated 
        % _n: Animal Late
        temp_s = X_ML_l==mid & O_ML_l==1;%& ((C_ML_l<start)) & C_ML_l>start-range2;% | (V_ML_l>start & V_ML_l<stop));%(V_ML_l>start & V_ML_l<stop) & ~(C_ML_l>start & C_ML_l<stop);%
        temp_n = X_ML_l==mid & O_ML_l==3;%& ((C_ML_l>stop)) & C_ML_l<stop+range2;% | (V_ML_l>start & V_ML_l<stop));%~(V_ML_l>start & V_ML_l<stop) & ((C_ML_l>start & C_ML_l<stop));%
        
        V_all.la = [V_all.la V_ML_l(temp)'];
        V_all.ls = [V_all.ls V_ML_l(temp_s)'];
        V_all.ln = [V_all.ln V_ML_l(temp_n)'];
        
        C_all.la = [C_all.la C_ML_l(temp)'];
        C_all.ls = [C_all.ls C_ML_l(temp_s)'];
        C_all.ln = [C_all.ln C_ML_l(temp_n)'];
        
        X_all.la = [X_all.la X_ML_l(temp)'];
        X_all.ls = [X_all.ls X_ML_l(temp_s)'];
        X_all.ln = [X_all.ln X_ML_l(temp_n)'];
        
        % Normal
        X_ML_n = CA(iseries).X_norm;
        V_ML_n = circshift(V1(iseries).MAP.norm(1:length(X_ML_n))',delay);
        C_ML_n = CA(iseries).MAP.norm(1:length(X_ML_n))';
        O_ML_n = outcome_n(1:length(X_ML_n));
        
        temp = X_ML_n==mid & O_ML_n==2;%;% & (C_ML_n>start & C_ML_n<stop);%  = (X_ML_n>start & X_ML_n<stop);
        temp_s = X_ML_n==mid & O_ML_n==1;% & ((C_ML_n<start))& C_ML_n>start-range2;% | (V_ML_n>start & V_ML_n<stop));%(V_ML_n>start & V_ML_n<stop) & ~(C_ML_n>start & C_ML_n<stop);%
        temp_n = X_ML_n==mid & O_ML_n==3;% & ((C_ML_n>stop))& C_ML_n<stop+range2;% | (V_ML_n>start & V_ML_n<stop));%~(V_ML_n>start & V_ML_n<stop) & (C_ML_n>start & C_ML_n<stop);%
        
        V_all.na = [V_all.na V_ML_n(temp)'];
        V_all.ns = [V_all.ns V_ML_n(temp_s)'];
        V_all.nn = [V_all.nn V_ML_n(temp_n)'];
        
        C_all.na = [C_all.na C_ML_n(temp)'];
        C_all.ns = [C_all.ns C_ML_n(temp_s)'];
        C_all.nn = [C_all.nn C_ML_n(temp_n)'];
        
        X_all.na = [X_all.na X_ML_n(temp)'];
        X_all.ns = [X_all.ns X_ML_n(temp_s)'];
        X_all.nn = [X_all.nn X_ML_n(temp_n)'];
        
        % High
        X_ML_h = CA(iseries).X_high_orig;
        V_ML_h = circshift(V1(iseries).MAP.high(1:length(X_ML_h))',delay);
        C_ML_h = CA(iseries).MAP.high(1:length(X_ML_h))';
        O_ML_h = outcome_h(1:length(X_ML_h));
        
        temp = X_ML_h==mid & O_ML_h==2;%;% & (C_ML_h>start & C_ML_h<stop);%  = (X_ML_h>start & X_ML_h<stop);
        temp_s = X_ML_h==mid & O_ML_h==1;% & ((C_ML_h<start)) & C_ML_h>start-range2;% | (V_ML_h>start & V_ML_h<stop));%(V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop);%
        temp_n = X_ML_h==mid & O_ML_h==3;% & ((C_ML_h>stop)) & C_ML_h<stop+range2;% & ~(C_ML_h>start & C_ML_h<stop);%~((C_ML_h>start & C_ML_h<stop) & (V_ML_h>start & V_ML_h<stop));%
        
        V_all.ha = [V_all.ha V_ML_h(temp)'];
        V_all.hs = [V_all.hs V_ML_h(temp_s)'];
        V_all.hn = [V_all.hn V_ML_h(temp_n)'];
        
        C_all.ha = [C_all.ha C_ML_h(temp)'];
        C_all.hs = [C_all.hs C_ML_h(temp_s)'];
        C_all.hn = [C_all.hn C_ML_h(temp_n)'];
        
        X_all.ha = [X_all.ha X_ML_h(temp)'];
        X_all.hs = [X_all.hs X_ML_h(temp_s)'];
        X_all.hn = [X_all.hn X_ML_h(temp_n)'];
    end
    
    figure(20)
    subplot(3,3,1)
    [F, c1, c2] = smoothhist2D_corrected([V_all.ls' X_all.ls'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ls)/(length(V_all.ls)+length(V_all.ln)));
    H.ls(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,4)
    [F, c1, c2] = smoothhist2D_corrected([V_all.ln' X_all.ln'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ln)/(length(V_all.ls)+length(V_all.ln)));
    H.ln(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,7)
    [F, c1, c2] = smoothhist2D_corrected([V_all.la' X_all.la'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.la)/(length(V_all.la)+length(V_all.la)));
    H.la(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,2)
    [F, c1, c2] = smoothhist2D_corrected([V_all.ns' X_all.ns'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ns)/(length(V_all.ns)+length(V_all.nn)));
    H.ns(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,5)
    [F, c1, c2] = smoothhist2D_corrected([V_all.nn' X_all.nn'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.nn)/(length(V_all.ns)+length(V_all.nn)));
    H.nn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,8)
    [F, c1, c2] = smoothhist2D_corrected([V_all.na' X_all.na'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.na)/(length(V_all.na)+length(V_all.na)));
    H.na(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,3)
    [F, c1, c2] = smoothhist2D_corrected([V_all.hs' X_all.hs'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.hs)/(length(V_all.hs)+length(V_all.hn)));
    H.hs(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,6)
    [F, c1, c2] = smoothhist2D_corrected([V_all.hn' X_all.hn'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.hn)/(length(V_all.hs)+length(V_all.hn)));
    H.hn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,9)
    [F, c1, c2] = smoothhist2D_corrected([V_all.ha' X_all.ha'],win,[50 50], 1:50);
%     F = smoothhist2D_corrected([V_all.ha' C_all.ha'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ha)/(length(V_all.ha)+length(V_all.ha)));
    H.ha(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,1)
    [F, c1, c2] = smoothhist2D_corrected([C_all.ls' X_all.ls'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.ls)/(length(X_all.ls)+length(X_all.ln)));
    H.xls(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,4)
    [F, c1, c2] = smoothhist2D_corrected([C_all.ln' X_all.ln'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.ln)/(length(X_all.ls)+length(X_all.ln)));
    H.xln(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,7)
    [F, c1, c2] = smoothhist2D_corrected([C_all.la' X_all.la'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.la)/(length(X_all.la)+length(X_all.la)));
    H.xla(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,2)
    [F, c1, c2] = smoothhist2D_corrected([C_all.ns' X_all.ns'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.ns)/(length(X_all.ns)+length(X_all.nn)));
    H.xns(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,5)
    [F, c1, c2] = smoothhist2D_corrected([C_all.nn' X_all.nn'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.nn)/(length(X_all.ns)+length(X_all.nn)));
    H.xnn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,8)
    [F, c1, c2] = smoothhist2D_corrected([C_all.na' X_all.na'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.na)/(length(X_all.na)+length(X_all.na)));
    H.xna(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,3)
    [F, c1, c2] = smoothhist2D_corrected([C_all.hs' X_all.hs'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.hs)/(length(X_all.hs)+length(X_all.hn)));
    H.xhs(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,6)
    [F, c1, c2] = smoothhist2D_corrected([C_all.hn' X_all.hn'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.hn)/(length(X_all.hs)+length(X_all.hn)));
    H.xhn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,9)
    [F, c1, c2] = smoothhist2D_corrected([C_all.ha' X_all.ha'],win,[50 50], 1:50);
%     F = smoothhist2D_corrected([X_all.ha' C_all.ha'],win,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.ha)/(length(X_all.ha)+length(X_all.ha)));
    H.xha(iseries,:,:) = F; %./ mean(F(:));
    for n = 1:9
        subplot(3,3,n);
        axis xy; axis equal; axis square; axis([1+range 50-range 1+range 50-range]); line(xlim,ylim);
        colorbar
    end
    %     drawnow
    %     pause
%     H.xhn(iseries,50-1*range:end,:) = nan;
%     H.xln(iseries,50-1*range:end,:) = nan;
%     H.xnn(iseries,50-range:end,:) = nan;
%     
%     H.hn(iseries,50-1*range:end,:) = nan;
%     H.ln(iseries,50-1*range:end,:) = nan;
%     H.nn(iseries,50-1*range:end,:) = nan;
%     
%     H.xhs(iseries,1:1*range,:) = nan;
%     H.xls(iseries,1:1*range,:) = nan;
%     H.xns(iseries,1:1*range,:) = nan;
%     
%     H.hs(iseries,1:1*range,:) = nan;
%     H.ls(iseries,1:1*range,:) = nan;
%     H.ns(iseries,1:1*range,:) = nan;
%     
    H.xhs(iseries,50-(20/15)*range:end,:) = nan;
    H.xls(iseries,50-(20/15)*range:end,:) = nan;
    H.xns(iseries,50-(20/15)*range:end,:) = nan;
    
    H.hs(iseries,50-(20/15)*range:end,:) = nan;
    H.ls(iseries,50-(20/15)*range:end,:) = nan;
    H.ns(iseries,50-(20/15)*range:end,:) = nan;
%     
%     H.xhn(iseries,1:1*range,:) = nan;
%     H.xln(iseries,1:1*range,:) = nan;
%     H.xnn(iseries,1:1*range,:) = nan;
%     
%     H.hn(iseries,1:1*range,:) = nan;
%     H.ln(iseries,1:1*range,:) = nan;
%     H.nn(iseries,1:1*range,:) = nan;
end
%%
clear marginals
% range = range*2;

for iseries = 1:7
    [marginals.la(iseries,:), X] = get45Marginal(sq(H.la(iseries,:,:))',50);
    [marginals.ls(iseries,:), X] = get45Marginal(sq(H.ls(iseries,:,:))',50);
    [marginals.ln(iseries,:), X] = get45Marginal(sq(H.ln(iseries,:,:))',50);
    
    [marginals.na(iseries,:), X] = get45Marginal(sq(H.na(iseries,:,:))',50);
    [marginals.ns(iseries,:), X] = get45Marginal(sq(H.ns(iseries,:,:))',50);
    [marginals.nn(iseries,:), X] = get45Marginal(sq(H.nn(iseries,:,:))',50);
    
    [marginals.ha(iseries,:), X] = get45Marginal(sq(H.ha(iseries,:,:))',50);
    [marginals.hs(iseries,:), X] = get45Marginal(sq(H.hs(iseries,:,:))',50);
    [marginals.hn(iseries,:), X] = get45Marginal(sq(H.hn(iseries,:,:))',50);
    
    [marginals.xla(iseries,:), X] = get45Marginal(sq(H.xla(iseries,:,:))',50);
    [marginals.xls(iseries,:), X] = get45Marginal(sq(H.xls(iseries,:,:))',50);
    [marginals.xln(iseries,:), X] = get45Marginal(sq(H.xln(iseries,:,:))',50);
    
    [marginals.xna(iseries,:), X] = get45Marginal(sq(H.xna(iseries,:,:))',50);
    [marginals.xns(iseries,:), X] = get45Marginal(sq(H.xns(iseries,:,:))',50);
    [marginals.xnn(iseries,:), X] = get45Marginal(sq(H.xnn(iseries,:,:))',50);
    
    [marginals.xha(iseries,:), X] = get45Marginal(sq(H.xha(iseries,:,:))',50);
    [marginals.xhs(iseries,:), X] = get45Marginal(sq(H.xhs(iseries,:,:))',50);
    [marginals.xhn(iseries,:), X] = get45Marginal(sq(H.xhn(iseries,:,:))',50);
    
end
%% Plotting V and X
figure(-delay+10)
subplot(3,4,1); imagesc(c1, c2,sq(nanmean(H.la(:,:,:),1))')
title('Animal Correct(low)')
subplot(3,4,2); imagesc(c1, c2,sq(nanmean(H.na(:,:,:),1))')
title('Animal Correct(med)')
subplot(3,4,3); imagesc(c1, c2,sq(nanmean(H.ha(:,:,:),1))')
title('Animal Correct(high)')
subplot(3,4,5); imagesc(c1, c2,sq(nanmean(H.ls(:,:,:),1))')
title('Animal Early (low)')
subplot(3,4,6); imagesc(c1, c2,sq(nanmean(H.ns(:,:,:),1))')
title('Animal Early (med)')
subplot(3,4,7); imagesc(c1, c2,sq(nanmean(H.hs(:,:,:),1))')
title('Animal Early (high)')
subplot(3,4,9); imagesc(c1, c2,sq(nanmean(H.ln(:,:,:),1))')
title('Animal Late (low)')
subplot(3,4,10); imagesc(c1, c2,sq(nanmean(H.nn(:,:,:),1))')
title('Animal Late (med)')
subplot(3,4,11); imagesc(c1, c2,sq(nanmean(H.hn(:,:,:),1))')
title('Animal Late (high)')
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
    set(gca,'clim',[0.5 3.0])
    colorbar;
    colormap(jet);
    xlabel('Position');
    ylabel('V1 decoded position');
    hold off
end


for n = 4:4:12
    subplot(3,4,n);
    axis tight;
    set(gca, 'xlim',[-30 30])
    line(xlim, [0 0], 'color','k','linestyle','--')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line([0 0], ylim, 'color','k','linestyle','--');
    line(xlim, [1 1], 'color','k','linestyle','--');
end
%% Plotting X and C
figure(-delay+11)
subplot(3,4,1); imagesc(c1, c2,sq(nanmean(H.xla(:,:,:),1))')
title('Animal Correct(low)')
subplot(3,4,2); imagesc(c1, c2,sq(nanmean(H.xna(:,:,:),1))')
title('Animal Correct(med)')
subplot(3,4,3); imagesc(c1, c2,sq(nanmean(H.xha(:,:,:),1))')
title('Animal Correct(high)')
subplot(3,4,5); imagesc(c1, c2,sq(nanmean(H.xls(:,:,:),1))')
title('Animal Early (low)')
subplot(3,4,6); imagesc(c1, c2,sq(nanmean(H.xns(:,:,:),1))')
title('Animal Early (med)')
subplot(3,4,7); imagesc(c1, c2,sq(nanmean(H.xhs(:,:,:),1))')
title('Animal Early (high)')
subplot(3,4,9); imagesc(c1, c2,sq(nanmean(H.xln(:,:,:),1))')
title('Animal Late (low)')
subplot(3,4,10); imagesc(c1, c2,sq(nanmean(H.xnn(:,:,:),1))')
title('Animal Late (med)')
subplot(3,4,11); imagesc(c1, c2,sq(nanmean(H.xhn(:,:,:),1))')
title('Animal Late (high)')
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
    set(gca,'clim',[0.5 3])
    colorbar;
    colormap(jet);
    xlabel('Position');
    ylabel('CA1 decoded position');
    hold off
end


for n = 4:4:12
    subplot(3,4,n);
    axis tight;
    set(gca, 'xlim',[-30 30])
    line(xlim, [0 0], 'color','k','linestyle','--')
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line([0 0], ylim, 'color','k','linestyle','--');
    line(xlim, [1 1], 'color','k','linestyle','--');
end