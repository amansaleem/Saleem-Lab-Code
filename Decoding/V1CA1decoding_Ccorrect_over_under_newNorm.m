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

shuffle=0; 
% 0-no shuffle; 
% 1:shuffle within bin; 
% 2:full shuffle, 
% 3:Completely V1 to be independent
% 4:Completely V1 & CA1 to be independent

switch shuffle
    case 0
        display('Actual (no shuffle )');
    case 1
        display('Shuffle within position bin');
    case 2
        display('Fully random');
    case 3
        display('Fully independent V1');
    case 4
        display('Fully independent V1 and CA1');
end
for iseries = 1:8
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
        
        if shuffle>=3
            % Just X + noise;
            V_ML_l = X_ML_l + 10*randn(size(X_ML_l));
            V_ML_l(V_ML_l>50) = V_ML_l(V_ML_l>50) - 50;
            V_ML_l(V_ML_l<1) = 50 + V_ML_l(V_ML_l<1);
            if shuffle==4
                C_ML_l = X_ML_l + 10*randn(size(X_ML_l));
                C_ML_l(C_ML_l>50) = C_ML_l(C_ML_l>50) - 50;
                C_ML_l(C_ML_l<1) = 50 + C_ML_l(C_ML_l<1);
            end
        end
        
        %times when the animal was within a certain range and CA1 was
        %correct
        temp = X_ML_l==mid;% & (C_ML_l>start & C_ML_l<stop);% (X_ML_l>start & X_ML_l<stop);
        % _s: CA1 underestimated 
        % _n: CA1 overshot
        temp_s = temp & ((C_ML_l<start)) & C_ML_l>start-range2;% | (V_ML_l>start & V_ML_l<stop));%(V_ML_l>start & V_ML_l<stop) & ~(C_ML_l>start & C_ML_l<stop);%
        temp_n = temp & ((C_ML_l>stop)) & C_ML_l<stop+range2;% | (V_ML_l>start & V_ML_l<stop));%~(V_ML_l>start & V_ML_l<stop) & ((C_ML_l>start & C_ML_l<stop));%
        
        if shuffle==2
            % Complete Shuffle of Low Contrast, no relation between X and V_ML
            Rperm = randperm(length(V_ML_l));
            V_ML_l = V_ML_l(Rperm);
            C_ML_l = C_ML_l(Rperm);
        elseif shuffle==1
            % Shuffle only current X of Low Contrast, no relation between X and V_ML
            Rperm = randperm(sum(temp));
            PermSubset = V_ML_l(temp);
            V_ML_l(temp) = PermSubset(Rperm);
            PermSubset = C_ML_l(temp);
            C_ML_l(temp) = PermSubset(Rperm);
        end
        V_all.la = [V_all.la V_ML_l(temp)'];
        V_all.ls = [V_all.ls V_ML_l(temp_s)'];
        V_all.ln = [V_all.ln V_ML_l(temp_n)'];

        C_all.la = [C_all.la C_ML_l(temp)'];
        C_all.ls = [C_all.ls C_ML_l(temp_s)'];
        C_all.ln = [C_all.ln C_ML_l(temp_n)'];
        
        X_all.la = [X_all.la X_ML_l(temp)'];
        X_all.ls = [X_all.ls X_ML_l(temp_s)'];
        X_all.ln = [X_all.ln X_ML_l(temp_n)'];
        
        V_all.aa = [V_all.aa V_ML_l(temp)'];
        V_all.as = [V_all.as V_ML_l(temp_s)'];
        V_all.an = [V_all.an V_ML_l(temp_n)'];
        
        C_all.aa = [C_all.aa C_ML_l(temp)'];
        C_all.as = [C_all.as C_ML_l(temp_s)'];
        C_all.an = [C_all.an C_ML_l(temp_n)'];
        
        X_all.aa = [X_all.aa X_ML_l(temp)'];
        X_all.as = [X_all.as X_ML_l(temp_s)'];
        X_all.an = [X_all.an X_ML_l(temp_n)'];
        
        % Normal
        X_ML_n = V1(iseries).X_norm;
        V_ML_n = circshift(V1(iseries).MAP.norm(1:length(X_ML_n))',delay);
        C_ML_n = CA(iseries).MAP.norm(1:length(X_ML_n))';
        
        if shuffle>=3
            % Just X + noise;
            V_ML_n = X_ML_n + 7*randn(size(X_ML_n));
            V_ML_n(V_ML_n>50) = V_ML_n(V_ML_n>50) - 50;
            V_ML_n(V_ML_n<1) = 50 + V_ML_n(V_ML_n<1);
            if shuffle==4
                C_ML_n = X_ML_n + 10*randn(size(X_ML_n));
                C_ML_n(C_ML_n>50) = C_ML_n(C_ML_n>50) - 50;
                C_ML_n(C_ML_n<1) = 50 + C_ML_n(C_ML_n<1);
            end
        end
        
        temp   = X_ML_n==mid;% & (C_ML_n>start & C_ML_n<stop);%  = (X_ML_n>start & X_ML_n<stop);
        temp_s = temp & ((C_ML_n<start))& C_ML_n>start-range2;% | (V_ML_n>start & V_ML_n<stop));%(V_ML_n>start & V_ML_n<stop) & ~(C_ML_n>start & C_ML_n<stop);%
        temp_n = temp & ((C_ML_n>stop))& C_ML_n<stop+range2;% | (V_ML_n>start & V_ML_n<stop));%~(V_ML_n>start & V_ML_n<stop) & (C_ML_n>start & C_ML_n<stop);%
        
        if shuffle==2
            % Complete Shuffle of Low Contrast, no relation between X and V_ML
            Rperm = randperm(length(V_ML_n));
            V_ML_n = V_ML_n(Rperm);
            C_ML_n = C_ML_n(Rperm);
        elseif shuffle==1
            % Shuffle only current X of Low Contrast, no relation between X and V_ML
            Rperm = randperm(sum(temp));
            PermSubset = V_ML_n(temp);
            V_ML_n(temp) = PermSubset(Rperm);
            PermSubset = C_ML_n(temp);
            C_ML_n(temp) = PermSubset(Rperm);
        end
        
        V_all.na = [V_all.na V_ML_n(temp)'];
        V_all.ns = [V_all.ns V_ML_n(temp_s)'];
        V_all.nn = [V_all.nn V_ML_n(temp_n)'];
        
        C_all.na = [C_all.na C_ML_n(temp)'];
        C_all.ns = [C_all.ns C_ML_n(temp_s)'];
        C_all.nn = [C_all.nn C_ML_n(temp_n)'];
        
        X_all.na = [X_all.na X_ML_n(temp)'];
        X_all.ns = [X_all.ns X_ML_n(temp_s)'];
        X_all.nn = [X_all.nn X_ML_n(temp_n)'];
        
        V_all.aa = [V_all.aa V_ML_n(temp)'];
        V_all.as = [V_all.as V_ML_n(temp_s)'];
        V_all.an = [V_all.an V_ML_n(temp_n)'];
        
        C_all.aa = [C_all.aa C_ML_n(temp)'];
        C_all.as = [C_all.as C_ML_n(temp_s)'];
        C_all.an = [C_all.an C_ML_n(temp_n)'];
        
        X_all.aa = [X_all.aa X_ML_n(temp)'];
        X_all.as = [X_all.as X_ML_n(temp_s)'];
        X_all.an = [X_all.an X_ML_n(temp_n)'];
        
        % High
        X_ML_h = V1(iseries).X_high_orig;
        V_ML_h = circshift(V1(iseries).MAP.high(1:length(X_ML_h))',delay);
        C_ML_h = CA(iseries).MAP.high(1:length(X_ML_h))';
        
        if shuffle>=3
            % Just X + noise;
            V_ML_h = X_ML_h + 4*randn(size(X_ML_h));
            V_ML_h(V_ML_h>50) = V_ML_h(V_ML_h>50) - 50;
            V_ML_h(V_ML_h<1) = 50 + V_ML_h(V_ML_h<1);
            if shuffle==4
                C_ML_h = X_ML_h + 10*randn(size(X_ML_h));
                C_ML_h(C_ML_h>50) = C_ML_h(C_ML_h>50) - 50;
                C_ML_h(C_ML_h<1) = 50 + C_ML_h(C_ML_h<1);
            end
        end
        temp = X_ML_h==mid;% & (C_ML_h>start & C_ML_h<stop);%  = (X_ML_h>start & X_ML_h<stop);
        temp_s = temp & ((C_ML_h<start)) & C_ML_h>start-range2;% | (V_ML_h>start & V_ML_h<stop));%(V_ML_h>start & V_ML_h<stop) & ~(C_ML_h>start & C_ML_h<stop);%
        temp_n = temp & ((C_ML_h>stop)) & C_ML_h<stop+range2;% & ~(C_ML_h>start & C_ML_h<stop);%~((C_ML_h>start & C_ML_h<stop) & (V_ML_h>start & V_ML_h<stop));%
        
        if shuffle==2
            % Complete Shuffle of Low Contrast, no relation between X and V_ML
            Rperm = randperm(length(V_ML_h));
            V_ML_h = V_ML_h(Rperm);
        elseif shuffle==1
            % Shuffle only current X of Low Contrast, no relation between X and V_ML
            Rperm = randperm(sum(temp));
            PermSubset = V_ML_h(temp);
            V_ML_h(temp) = PermSubset(Rperm);
        end
        
        V_all.ha = [V_all.ha V_ML_h(temp)'];
        V_all.hs = [V_all.hs V_ML_h(temp_s)'];
        V_all.hn = [V_all.hn V_ML_h(temp_n)'];
        
        C_all.ha = [C_all.ha C_ML_h(temp)'];
        C_all.hs = [C_all.hs C_ML_h(temp_s)'];
        C_all.hn = [C_all.hn C_ML_h(temp_n)'];
        
        X_all.ha = [X_all.ha X_ML_h(temp)'];
        X_all.hs = [X_all.hs X_ML_h(temp_s)'];
        X_all.hn = [X_all.hn X_ML_h(temp_n)'];
        
        
        V_all.aa = [V_all.aa V_ML_h(temp)'];
        V_all.as = [V_all.as V_ML_h(temp_s)'];
        V_all.an = [V_all.an V_ML_h(temp_n)'];
        
        C_all.aa = [C_all.aa C_ML_h(temp)'];
        C_all.as = [C_all.as C_ML_h(temp_s)'];
        C_all.an = [C_all.an C_ML_h(temp_n)'];
        
        X_all.aa = [X_all.aa X_ML_h(temp)'];
        X_all.as = [X_all.as X_ML_h(temp_s)'];
        X_all.an = [X_all.an X_ML_h(temp_n)'];
    end
    
    figure(20)
    subplot(3,3,1)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.ls' X_all.ls'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ls)/(length(V_all.ls)+length(V_all.ln)));
    H.ls(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,4)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.ln' X_all.ln'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ln)/(length(V_all.ls)+length(V_all.ln)));
    H.ln(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,7)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.la' X_all.la'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.la)/(length(V_all.la)+length(V_all.la)));
    H.la(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,2)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.ns' X_all.ns'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ns)/(length(V_all.ns)+length(V_all.nn)));
    H.ns(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,5)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.nn' X_all.nn'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.nn)/(length(V_all.ns)+length(V_all.nn)));
    H.nn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,8)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.na' X_all.na'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.na)/(length(V_all.na)+length(V_all.na)));
    H.na(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,3)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.hs' X_all.hs'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.hs)/(length(V_all.hs)+length(V_all.hn)));
    H.hs(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,6)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.hn' X_all.hn'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.hn)/(length(V_all.hs)+length(V_all.hn)));
    H.hn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,9)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.ha' X_all.ha'],4,[50 50], 1:50);
%     F = smoothhist2D_SVDcorrected([V_all.ha' C_all.ha'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ha)/(length(V_all.ha)+length(V_all.ha)));
    H.ha(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,1)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.ls' X_all.ls'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.ls)/(length(X_all.ls)+length(X_all.ln)));
    H.xls(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,4)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.ln' X_all.ln'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.ln)/(length(X_all.ls)+length(X_all.ln)));
    H.xln(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,7)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.la' X_all.la'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.la)/(length(X_all.la)+length(X_all.la)));
    H.xla(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,2)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.ns' X_all.ns'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.ns)/(length(X_all.ns)+length(X_all.nn)));
    H.xns(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,5)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.nn' X_all.nn'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.nn)/(length(X_all.ns)+length(X_all.nn)));
    H.xnn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,8)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.na' X_all.na'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.na)/(length(X_all.na)+length(X_all.na)));
    H.xna(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,3)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.hs' X_all.hs'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.hs)/(length(X_all.hs)+length(X_all.hn)));
    H.xhs(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,6)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.hn' X_all.hn'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.hn)/(length(X_all.hs)+length(X_all.hn)));
    H.xhn(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,9)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.ha' X_all.ha'],4,[50 50], 1:50);
%     F = smoothhist2D_SVDcorrected([X_all.ha' C_all.ha'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(X_all.ha)/(length(X_all.ha)+length(X_all.ha)));
    H.xha(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,1)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.as' X_all.as'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ls)/(length(V_all.ls)+length(V_all.ln)));
    H.as(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,4)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.an' X_all.an'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ln)/(length(V_all.ls)+length(V_all.ln)));
    H.an(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,7)
    [F, c1, c2] = smoothhist2D_SVDcorrected([V_all.aa' X_all.aa'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.la)/(length(V_all.la)+length(V_all.la)));
    H.aa(iseries,:,:) = F; %./ mean(F(:));
    
    subplot(3,3,1)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.as' X_all.as'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ls)/(length(V_all.ls)+length(V_all.ln)));
    H.xas(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,4)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.an' X_all.an'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.ln)/(length(V_all.ls)+length(V_all.ln)));
    H.xan(iseries,:,:) = F; %./ mean(F(:));
    subplot(3,3,7)
    [F, c1, c2] = smoothhist2D_SVDcorrected([C_all.aa' X_all.aa'],4,[50 50], 1:50);
    imagesc(c1, c2,F*2500 * length(V_all.la)/(length(V_all.la)+length(V_all.la)));
    H.xaa(iseries,:,:) = F; %./ mean(F(:));
    
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
%     H.xhs(iseries,50-1*range:end,:) = nan;
%     H.xls(iseries,50-1*range:end,:) = nan;
%     H.xns(iseries,50-1*range:end,:) = nan;
%     
%     H.hs(iseries,50-1*range:end,:) = nan;
%     H.ls(iseries,50-1*range:end,:) = nan;
%     H.ns(iseries,50-1*range:end,:) = nan;
%     
%     H.xhn(iseries,1:1*range,:) = nan;
%     H.xln(iseries,1:1*range,:) = nan;
%     H.xnn(iseries,1:1*range,:) = nan;
%     
%     H.hn(iseries,1:1*range,:) = nan;
%     H.ln(iseries,1:1*range,:) = nan;
%     H.nn(iseries,1:1*range,:) = nan;
%     
%     H.as(iseries,1:1*range,:) = nan;
%     H.as(iseries,50-1*range:end,:) = nan;
%     H.xas(iseries,1:1*range,:) = nan;
%     H.xas(iseries,50-1*range:end,:) = nan;
%     
%     H.an(iseries,1:1*range,:) = nan;
%     H.an(iseries,50-1*range:end,:) = nan;
%     H.xan(iseries,1:1*range,:) = nan;
%     H.xan(iseries,50-1*range:end,:) = nan;
%     
%     H.aa(iseries,1:1*range,:) = nan;
%     H.aa(iseries,50-1*range:end,:) = nan;
%     H.xaa(iseries,1:1*range,:) = nan;
%     H.xaa(iseries,50-1*range:end,:) = nan;
%     
end
%%
clear marginals
% range = range*2;

for iseries = 1:8
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
    
    [marginals.aa(iseries,:), X] = get45Marginal(sq(H.aa(iseries,:,:))',50);
    [marginals.as(iseries,:), X] = get45Marginal(sq(H.as(iseries,:,:))',50);
    [marginals.an(iseries,:), X] = get45Marginal(sq(H.an(iseries,:,:))',50);
    
    [marginals.xaa(iseries,:), X] = get45Marginal(sq(H.xaa(iseries,:,:))',50);
    [marginals.xas(iseries,:), X] = get45Marginal(sq(H.xas(iseries,:,:))',50);
    [marginals.xan(iseries,:), X] = get45Marginal(sq(H.xan(iseries,:,:))',50);
    
end
%% Plotting V and X
figure(-delay+110+shuffle)
subplot(3,4,1); imagesc(c1, c2,sq(mean(H.la,1))')
title('All (low)')
subplot(3,4,2); imagesc(c1, c2,sq(mean(H.na,1))')
title('All (med)')
subplot(3,4,3); imagesc(c1, c2,sq(mean(H.ha,1))')
title('All (high)')
subplot(3,4,5); imagesc(c1, c2,sq(mean(H.ls,1))')
title('CA1 undershot (low)')
subplot(3,4,6); imagesc(c1, c2,sq(mean(H.ns,1))')
title('CA1 undershot (med)')
subplot(3,4,7); imagesc(c1, c2,sq(mean(H.hs,1))')
title('CA1 undershot (high)')
subplot(3,4,9); imagesc(c1, c2,sq(mean(H.ln,1))')
title('CA1 overshot (low)')
subplot(3,4,10); imagesc(c1, c2,sq(mean(H.nn,1))')
title('CA1 overshot (med)')
subplot(3,4,11); imagesc(c1, c2,sq(mean(H.hn,1))')
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
figure(-delay+213)
subplot(3,2,1); imagesc(c1, c2,sq(mean(H.aa,1))')
title('All (low)')
subplot(3,2,3); imagesc(c1, c2,sq(mean(H.as./H.aa,1))')
title('CA1 undershot (low)')
subplot(3,2,5); imagesc(c1, c2,sq(mean(H.an./H.aa,1))')
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