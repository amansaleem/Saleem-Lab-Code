if ~exist('CA')
%     load Posterior_all_150226_trainedCorrect_gaussSmth_250_noInter
    load Posterior_all_160121_trainedCV_250
    CA = Posterior_all;
    
    load Posterior_all_160121_V1_trainedCV_250
    V1 = Posterior_all;
    clear Post*
end

only_plots = 0;

if ~only_plots
    delay = 0;
    Noise_norm_l   = zeros(5,8,50,50,50);
    Noise_norm_n   = zeros(5,8,50,50,50);
    Noise_norm_h   = zeros(5,8,50,50,50);
    
    Signal_norm_l = zeros(5,8,50,50,50);
    Signal_norm_n = zeros(5,8,50,50,50);
    Signal_norm_h = zeros(5,8,50,50,50);
    
    All_norm_l = zeros(5,8,50,50,50);
    All_norm_n = zeros(5,8,50,50,50);
    All_norm_h = zeros(5,8,50,50,50);
    
    for shuffle = 0:4;
        switch shuffle
            case 0
                display('Actual (no shuffle )');% 0-no shuffle;
            case 1
                display('Shuffle within position bin');% 1:shuffle within bin;
            case 2
                display('Fully random shuffle'); % 2:full shuffle,
            case 3
                display('Synthetic data: 0 noise correlations');% 3:V1 to be independent
            case 4
                display('Synthetic data: with noise correlations');% 4:Completely V1 & CA1 to be independent
        end
        figure;
        for iseries = 1:8
            X_ML_l = V1(iseries).X_low_orig;
            V_ML_l = V1(iseries).MAP.low(1:length(X_ML_l))';
            S_ML_l = V1(iseries).data.smthBallSpd(V1(iseries).t_low);
            C_ML_l = CA(iseries).MAP.low(1:length(X_ML_l))';
            %         V_ML_l = S_ML_l(1:length(X_ML_l));
            
            X_ML_h = V1(iseries).X_high_orig;
            V_ML_h = V1(iseries).MAP.high(1:length(X_ML_h))';
            S_ML_h = V1(iseries).data.smthBallSpd(V1(iseries).t_high);
            C_ML_h = CA(iseries).MAP.high(1:length(X_ML_h))';
            %         V_ML_h = S_ML_h(1:length(X_ML_h));
            
            X_ML_n = V1(iseries).X_norm;
            V_ML_n = V1(iseries).MAP.norm(1:length(X_ML_n))';
            S_ML_n = V1(iseries).data.smthBallSpd(V1(iseries).t_norm);
            C_ML_n = CA(iseries).MAP.norm(1:length(X_ML_n))';
            %         V_ML_n = S_ML_n(1:length(X_ML_n));
            
            if shuffle==3
                % Just X + noise;
                V_ML_l = X_ML_l + 12*randn(size(X_ML_l));
                V_ML_l(V_ML_l>50) = V_ML_l(V_ML_l>50) - 50;
                V_ML_l(V_ML_l<1) = 50 + V_ML_l(V_ML_l<1);
                V_ML_l(V_ML_l>50) = V_ML_l(V_ML_l>50) - 50;
                V_ML_l(V_ML_l<1) = 50 + V_ML_l(V_ML_l<1);
                
                V_ML_n = X_ML_n + 8*randn(size(X_ML_n));
                V_ML_n(V_ML_n>50) = V_ML_n(V_ML_n>50) - 50;
                V_ML_n(V_ML_n<1) = 50 + V_ML_n(V_ML_n<1);
                V_ML_n(V_ML_n>50) = V_ML_n(V_ML_n>50) - 50;
                V_ML_n(V_ML_n<1) = 50 + V_ML_n(V_ML_n<1);
                
                V_ML_h = X_ML_h + 5*randn(size(X_ML_h));
                V_ML_h(V_ML_h>50) = V_ML_h(V_ML_h>50) - 50;
                V_ML_h(V_ML_h<1) = 50 + V_ML_h(V_ML_h<1);
                V_ML_h(V_ML_h>50) = V_ML_h(V_ML_h>50) - 50;
                V_ML_h(V_ML_h<1) = 50 + V_ML_h(V_ML_h<1);
                
                C_ML_h = X_ML_h + 2*randn(size(X_ML_h));
                C_ML_h(C_ML_h>50) = C_ML_h(C_ML_h>50) - 50;
                C_ML_h(C_ML_h<1) = 50 + C_ML_h(C_ML_h<1);
                C_ML_h(C_ML_h>50) = C_ML_h(C_ML_h>50) - 50;
                C_ML_h(C_ML_h<1) = 50 + C_ML_h(C_ML_h<1);
                
                C_ML_n = X_ML_n + 4*randn(size(X_ML_n));
                C_ML_n(C_ML_n>50) = C_ML_n(C_ML_n>50) - 50;
                C_ML_n(C_ML_n<1) = 50 + C_ML_n(C_ML_n<1);
                C_ML_n(C_ML_n>50) = C_ML_n(C_ML_n>50) - 50;
                C_ML_n(C_ML_n<1) = 50 + C_ML_n(C_ML_n<1);
                
                C_ML_l = X_ML_l + 7*randn(size(X_ML_l));
                C_ML_l(C_ML_l>50) = C_ML_l(C_ML_l>50) - 50;
                C_ML_l(C_ML_l<1) = 50 + C_ML_l(C_ML_l<1);
                C_ML_l(C_ML_l>50) = C_ML_l(C_ML_l>50) - 50;
                C_ML_l(C_ML_l<1) = 50 + C_ML_l(C_ML_l<1);
            end
            if shuffle==4
                % X + noise + common_noise;
                Cnoise_l = 1.5*randn(size(X_ML_l));
                Cnoise_n = 1.5*randn(size(X_ML_n));
                Cnoise_h = 1.5*randn(size(X_ML_h));
                
                V_ML_l = X_ML_l + 12*randn(size(X_ML_l)) + Cnoise_l;
                V_ML_l(V_ML_l>50) = V_ML_l(V_ML_l>50) - 50;
                V_ML_l(V_ML_l<1) = 50 + V_ML_l(V_ML_l<1);
                V_ML_l(V_ML_l>50) = V_ML_l(V_ML_l>50) - 50;
                V_ML_l(V_ML_l<1) = 50 + V_ML_l(V_ML_l<1);
                
                V_ML_n = X_ML_n + 8*randn(size(X_ML_n)) + Cnoise_n;
                V_ML_n(V_ML_n>50) = V_ML_n(V_ML_n>50) - 50;
                V_ML_n(V_ML_n<1) = 50 + V_ML_n(V_ML_n<1);
                V_ML_n(V_ML_n>50) = V_ML_n(V_ML_n>50) - 50;
                V_ML_n(V_ML_n<1) = 50 + V_ML_n(V_ML_n<1);
                
                V_ML_h = X_ML_h + 5*randn(size(X_ML_h)) + Cnoise_h;
                V_ML_h(V_ML_h>50) = V_ML_h(V_ML_h>50) - 50;
                V_ML_h(V_ML_h<1) = 50 + V_ML_h(V_ML_h<1);
                V_ML_h(V_ML_h>50) = V_ML_h(V_ML_h>50) - 50;
                V_ML_h(V_ML_h<1) = 50 + V_ML_h(V_ML_h<1);
                
                C_ML_h = X_ML_h + 2*randn(size(X_ML_h)) + Cnoise_h;
                C_ML_h(C_ML_h>50) = C_ML_h(C_ML_h>50) - 50;
                C_ML_h(C_ML_h<1) = 50 + C_ML_h(C_ML_h<1);
                C_ML_h(C_ML_h>50) = C_ML_h(C_ML_h>50) - 50;
                C_ML_h(C_ML_h<1) = 50 + C_ML_h(C_ML_h<1);
                
                C_ML_n = X_ML_n + 4*randn(size(X_ML_n)) + Cnoise_n;
                C_ML_n(C_ML_n>50) = C_ML_n(C_ML_n>50) - 50;
                C_ML_n(C_ML_n<1) = 50 + C_ML_n(C_ML_n<1);
                C_ML_n(C_ML_n>50) = C_ML_n(C_ML_n>50) - 50;
                C_ML_n(C_ML_n<1) = 50 + C_ML_n(C_ML_n<1);
                
                C_ML_l = X_ML_l + 7*randn(size(X_ML_l)) + Cnoise_l;
                C_ML_l(C_ML_l>50) = C_ML_l(C_ML_l>50) - 50;
                C_ML_l(C_ML_l<1) = 50 + C_ML_l(C_ML_l<1);
                C_ML_l(C_ML_l>50) = C_ML_l(C_ML_l>50) - 50;
                C_ML_l(C_ML_l<1) = 50 + C_ML_l(C_ML_l<1);
            end
            if shuffle==2
                % Complete Shuffle of Contrast, no relation between X and V_ML,
                % C_ML
                Rperm = randperm(length(V_ML_l));
                V_ML_l = V_ML_l(Rperm);
                Rperm = randperm(length(V_ML_l));
                C_ML_l = C_ML_l(Rperm);
                
                Rperm = randperm(length(V_ML_n));
                V_ML_n = V_ML_n(Rperm);
                Rperm = randperm(length(V_ML_n));
                C_ML_n = C_ML_n(Rperm);
                
                Rperm = randperm(length(V_ML_h));
                V_ML_h = V_ML_h(Rperm);
                Rperm = randperm(length(V_ML_h));
                C_ML_h = C_ML_h(Rperm);
            end
            %%
            for mid = 1:50
                %times when the animal was at a certain position
                temp = X_ML_l==mid;%
                
                if shuffle==1
                    % Shuffle only current X of Low Contrast, no relation between X and V_ML
                    Rperm = randperm(sum(temp));
                    PermSubset = V_ML_l(temp);
                    V_ML_l(temp) = PermSubset(Rperm);
                    %                 Rperm = randperm(sum(temp));
                    %                 PermSubset = C_ML_l(temp);
                    %                 C_ML_l(temp) = PermSubset(Rperm);
                end
                [Noise_norm_l(shuffle+1, iseries,mid,:,:),~,~,...
                    Signal_norm_l(shuffle+1, iseries,mid,:,:),...
                    All_norm_l(shuffle+1, iseries,mid,:,:)]...
                    = smoothhist2D_SVDcorrected(...
                    [C_ML_l(temp) V_ML_l(temp)],4,[50 50], 1:50);
                
                % Normal
                temp   = X_ML_n==mid;% & (C_ML_n>start & C_ML_n<stop);%  = (X_ML_n>start & X_ML_n<stop);
                if shuffle==1
                    % Shuffle only current X of Low Contrast, no relation between X and V_ML
                    Rperm = randperm(sum(temp));
                    PermSubset = V_ML_n(temp);
                    V_ML_n(temp) = PermSubset(Rperm);
                    Rperm = randperm(sum(temp));
                    PermSubset = C_ML_n(temp);
                    C_ML_n(temp) = PermSubset(Rperm);
                end
                [Noise_norm_n(shuffle+1, iseries,mid,:,:),~,~,...
                    Signal_norm_n(shuffle+1, iseries,mid,:,:),...
                    All_norm_n(shuffle+1, iseries,mid,:,:)]...
                    = smoothhist2D_SVDcorrected(...
                    [C_ML_n(temp) V_ML_n(temp)],4,[50 50], 1:50);
                % High
                temp = X_ML_h==mid;% & (C_ML_h>start & C_ML_h<stop);%  = (X_ML_h>start & X_ML_h<stop);
                if shuffle==1
                    % Shuffle only current X of Low Contrast, no relation between X and V_ML
                    % Shuffle only current X of Low Contrast, no relation between X and V_ML
                    Rperm = randperm(sum(temp));
                    PermSubset = V_ML_h(temp);
                    V_ML_h(temp) = PermSubset(Rperm);
                    Rperm = randperm(sum(temp));
                    PermSubset = C_ML_h(temp);
                    C_ML_h(temp) = PermSubset(Rperm);
                end
                [Noise_norm_h(shuffle+1, iseries,mid,:,:),~,~,...
                    Signal_norm_h(shuffle+1, iseries,mid,:,:),...
                    All_norm_h(shuffle+1, iseries,mid,:,:)]...
                    = smoothhist2D_SVDcorrected(...
                    [C_ML_h(temp) V_ML_h(temp)],4,[50 50], 1:50);
            end
        end
        drawnow
        close
    end
end
%%
clear marginals
% range = range*2;
for shuffle=0:4
    for iseries = 1:8
        [marginals.l(shuffle+1,iseries,:), X2] = get45Marginal(sq(sum(Noise_norm_l(shuffle+1, iseries,:,:,:),3))',50);
        [marginals.n(shuffle+1,iseries,:), X2] = get45Marginal(sq(sum(Noise_norm_n(shuffle+1, iseries,:,:,:),3))',50);
        [marginals.h(shuffle+1,iseries,:), X2] = get45Marginal(sq(sum(Noise_norm_h(shuffle+1, iseries,:,:,:),3))',50);
    end
end
X = X2*2*sqrt(2);
figure(17)
for n = 2:5
    plot(X, mean(sq(marginals.l(n,:,:))), 'color',[.5 .5 .5])
    hold on;
    plot(X, mean(sq(marginals.n(n,:,:))), 'color',[.5 .5 .5])
    plot(X, mean(sq(marginals.h(n,:,:))), 'color',[.5 .5 .5])
end
errorarea_as(X, mean(sq(marginals.l(1,:,:))),sem(sq(marginals.l(1,:,:))),'b')
errorarea_as(X, mean(sq(marginals.n(1,:,:))),sem(sq(marginals.n(1,:,:))),'k')
errorarea_as(X, mean(sq(marginals.h(1,:,:))),sem(sq(marginals.h(1,:,:))),'r')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
set(gca,'XLim',[-75 75])
for n = 1:50
    p_l(n) = signrank(sq(marginals.l(1,:,n)));
    p_n(n) = signrank(sq(marginals.n(1,:,n)));
    p_h(n) = signrank(sq(marginals.h(1,:,n)));
end
%%
for shuffle = 0:4
    figure(16)
    subplot(5,3,3*(shuffle)+1)
    imagesc(sq(mean(mean(Noise_norm_l(shuffle+1, :,:,:,:),3),2))');
    axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
    title(['Noise']);
    subplot(5,3,3*(shuffle)+2)
    imagesc(sq(mean(mean(Noise_norm_n(shuffle+1, :,:,:,:),3),2))');
    axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
    title(['Norm, Shuffle = ' num2str(shuffle)]);
    subplot(5,3,3*(shuffle)+3)
    imagesc(sq(mean(mean(Noise_norm_h(shuffle+1, :,:,:,:),3),2))');
    axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
    
    switch shuffle
        case 0
            title('Actual (no shuffle )');
        case 1
            title('Shuffle within position bin V1');
        case 2
            title('Fully random');
        case 3
            title('Independent V1');
        case 4
            title('Fully independent V1 and CA1');
    end
end
for n = 1:15
    subplot(5,3,n)
    set(gca, 'clim', [-5 5]/1000)
end
%%
clear marginals
% range = range*2;
for shuffle=0:4
    for iseries = 1:8
        [marginals.l(shuffle+1,iseries,:), X2] = get45Marginal(sq(sum(Signal_norm_l(shuffle+1, iseries,:,:,:),3))',50);
        [marginals.n(shuffle+1,iseries,:), X2] = get45Marginal(sq(sum(Signal_norm_n(shuffle+1, iseries,:,:,:),3))',50);
        [marginals.h(shuffle+1,iseries,:), X2] = get45Marginal(sq(sum(Signal_norm_h(shuffle+1, iseries,:,:,:),3))',50);
    end
end
X = X2*2*sqrt(2);
figure(27)
for n = 2:5
    plot(X, mean(sq(marginals.l(n,:,:))), 'color',[.5 .5 .5])
    hold on;
    plot(X, mean(sq(marginals.n(n,:,:))), 'color',[.5 .5 .5])
    plot(X, mean(sq(marginals.h(n,:,:))), 'color',[.5 .5 .5])
end
errorarea_as(X, mean(sq(marginals.l(1,:,:))),sem(sq(marginals.l(1,:,:))),'b')
errorarea_as(X, mean(sq(marginals.n(1,:,:))),sem(sq(marginals.n(1,:,:))),'k')
errorarea_as(X, mean(sq(marginals.h(1,:,:))),sem(sq(marginals.h(1,:,:))),'r')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
set(gca,'XLim',[-75 75])
for n = 1:50
    p_l(n) = signrank(sq(marginals.l(1,:,n)));
    p_n(n) = signrank(sq(marginals.n(1,:,n)));
    p_h(n) = signrank(sq(marginals.h(1,:,n)));
end
%%
for shuffle = 0:4
    figure(12)
    subplot(5,3,3*(shuffle)+1)
    imagesc(sq(mean(mean(Signal_norm_l(shuffle+1, :,:,:,:),3),2))');
    axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
    title(['Signal']);
    subplot(5,3,3*(shuffle)+2)
    imagesc(sq(mean(mean(Signal_norm_n(shuffle+1, :,:,:,:),3),2))');
    axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
    title(['Norm, Shuffle = ' num2str(shuffle)]);
    subplot(5,3,3*(shuffle)+3)
    imagesc(sq(mean(mean(Signal_norm_h(shuffle+1, :,:,:,:),3),2))');
    axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
    
    switch shuffle
        case 0
            title('Actual (no shuffle )');
        case 1
            title('Shuffle within position bin V1');
        case 2
            title('Fully random');
        case 3
            title('Synthetic: 0 noise correlations');
        case 4
            title('Synthetic: with noise correlations');
    end
end
for n = 1:15
    subplot(5,3,n)
    set(gca, 'clim', [0 0.04])
end
% %%
% for iseries = 1:8
%     figure(iseries)
%     for shuffle = 0:4
%         subplot(5,3,3*(shuffle)+1)
%         imagesc(sq(mean(mean(Noise_norm_l(shuffle+1, iseries,:,:,:),3),2))');
%         axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
%         subplot(5,3,3*(shuffle)+2)
%         imagesc(sq(mean(mean(Noise_norm_n(shuffle+1, iseries,:,:,:),3),2))');
%         axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
%         title(['Norm, Shuffle = ' num2str(shuffle)]);
%         subplot(5,3,3*(shuffle)+3)
%         imagesc(sq(mean(mean(Noise_norm_h(shuffle+1, iseries,:,:,:),3),2))');
%         axis equal; axis xy; axis tight; line(xlim, ylim); colormap(jet); colorbar
%
%         switch shuffle
%             case 0
%                 title('Actual (no shuffle )');
%             case 1
%                 title('Shuffle within position bin V1');
%             case 2
%                 title('Fully random');
%             case 3
%                 title('Independent V1');
%             case 4
%                 title('Fully independent V1 and CA1');
%         end
%     end
%     for n = 1:15
%         subplot(5,3,n)
%         set(gca, 'clim', [-5 5]/10000)
%     end
% end
%%
% figure(30)
shuffle = 4;
sLim = max(max(mean(mean(sum(Signal_norm_l(1:3, :,n,:,:),3),2),1)))*2;
nLim = max(max(mean(mean(sum(Noise_norm_l(1:3, :,n,:,:),3),2),1)));

for n = 1:50;
    subplot(331)
    imagesc(sq(mean(mean(Signal_norm_l(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    set(gca, 'clim', [0 sLim])
    plot(n,n,'ko','MarkerSize',15)
    hold off
    title('Signal correlation, Low Contrast')
    subplot(332)
    imagesc(sq(mean(mean(Noise_norm_l(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    plot(n,n,'ko','MarkerSize',15)
    hold off
    set(gca, 'clim', [-1 1]*nLim)
    title('Noise correlation, Low Contrast')
    subplot(333)
    imagesc(sq(mean(mean(All_norm_l(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    set(gca, 'clim', [0 sLim])
    plot(n,n,'ko','MarkerSize',15)
    hold off
    title('All correlation, Low Contrast')
    
    subplot(334)
    imagesc(sq(mean(mean(Signal_norm_n(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    plot(n,n,'ko','MarkerSize',15)
    set(gca, 'clim', [0 sLim])
    hold off
    title('Signal correlation, Medium Contrast')
    subplot(335)
    imagesc(sq(mean(mean(Noise_norm_n(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    plot(n,n,'ko','MarkerSize',15)
    hold off
    set(gca, 'clim', [-1 1]*nLim)
    title('Noise correlation, Medium Contrast')
    subplot(336)
    imagesc(sq(mean(mean(All_norm_n(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    plot(n,n,'ko','MarkerSize',15)
    hold off
    set(gca, 'clim', [0 1]*sLim)
%     title('Noise correlation, Medium Contrast')
    
    subplot(337)
    imagesc(sq(mean(mean(Signal_norm_h(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    set(gca, 'clim', [0 sLim])
    plot(n,n,'ko','MarkerSize',15)
    hold off
    title('Signal correlation, High Contrast')
    subplot(338)
    imagesc(sq(mean(mean(Noise_norm_h(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    plot(n,n,'ko','MarkerSize',15)
    set(gca, 'clim', [-1 1]*nLim)
    hold off
    title('Noise correlation, High Contrast')
    subplot(339)
    imagesc(sq(mean(mean(All_norm_h(shuffle+1, :,n,:,:),3),2))'); axis xy; axis square
    hold on;line(xlim,ylim,'color','k')
    plot(n,n,'ko','MarkerSize',15)
    set(gca, 'clim', [0 1]*sLim)
    hold off
    title('All correlation, High Contrast')
    
    pause(0.05)
%     pause
    for n = 1:6
        subplot(3,2,6)
        colorbar off;
    end
end

%
subplot(321)
imagesc(sq(mean(mean(Signal_norm_l(shuffle+1, :,:,:,:),3),2))'); axis xy; axis square
hold on;line(xlim,ylim,'color','k')
set(gca, 'clim', [0 sLim*15])
hold off
title('Signal correlation, Low Contrast')
subplot(322)
imagesc(sq(mean(mean(Noise_norm_l(shuffle+1, :,:,:,:),3),2))'); axis xy; axis square
hold on;line(xlim,ylim,'color','k')
hold off
set(gca, 'clim', [-1 1]*nLim*5)
title('Noise correlation, Low Contrast')

subplot(323)
imagesc(sq(mean(mean(Signal_norm_n(shuffle+1, :,:,:,:),3),2))'); axis xy; axis square
hold on;line(xlim,ylim,'color','k')
set(gca, 'clim', [0 sLim*15])
hold off
title('Signal correlation, Medium Contrast')
subplot(324)
imagesc(sq(mean(mean(Noise_norm_n(shuffle+1, :,:,:,:),3),2))'); axis xy; axis square
hold on;line(xlim,ylim,'color','k')
hold off
set(gca, 'clim', [-1 1]*nLim*5)
title('Noise correlation, Medium Contrast')

subplot(325)
imagesc(sq(mean(mean(Signal_norm_h(shuffle+1, :,:,:,:),3),2))'); axis xy; axis square
hold on;line(xlim,ylim,'color','k')
set(gca, 'clim', [0 sLim*15])
hold off
title('Signal correlation, High Contrast')
subplot(326)
imagesc(sq(mean(mean(Noise_norm_h(shuffle+1, :,:,:,:),3),2))'); axis xy; axis square
hold on;line(xlim,ylim,'color','k')
set(gca, 'clim', [-1 1]*nLim*5)
hold off
title('Noise correlation, High Contrast')

for n = 1:6
    subplot(3,2,6)
    colorbar off;
end

%%
Centred_signal_l = zeros(5,8,50,100,100);%
Centred_signal_n = zeros(5,8,50,100,100);%
Centred_signal_h = zeros(5,8,50,100,100);%

Centred_noise_l = zeros(5,8,50,100,100);%
Centred_noise_n = zeros(5,8,50,100,100);%
Centred_noise_h = zeros(5,8,50,100,100);%

Centred_all_l = zeros(5,8,50,100,100);%
Centred_all_n = zeros(5,8,50,100,100);%
Centred_all_h = zeros(5,8,50,100,100);%

Centred_signal_l(:,:,:,51:100,51:100) = Signal_norm_l;%
Centred_signal_n(:,:,:,51:100,51:100) = Signal_norm_n;%
Centred_signal_h(:,:,:,51:100,51:100) = Signal_norm_h;%

Centred_all_l(:,:,:,51:100,51:100) = All_norm_l;%
Centred_all_n(:,:,:,51:100,51:100) = All_norm_n;%
Centred_all_h(:,:,:,51:100,51:100) = All_norm_h;%

Centred_noise_l(:,:,:,51:100,51:100) = Noise_norm_l;%
Centred_noise_n(:,:,:,51:100,51:100) = Noise_norm_n;%
Centred_noise_h(:,:,:,51:100,51:100) = Noise_norm_h;%

for ishuffle = 1:5
    for iseries = 1:8
        for bin = 1:50
            Centred_signal_l(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_signal_l(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
            Centred_signal_n(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_signal_n(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
            Centred_signal_h(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_signal_h(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
            
            Centred_noise_l(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_noise_l(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
            Centred_noise_n(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_noise_n(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
            Centred_noise_h(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_noise_h(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
            
            Centred_all_l(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_all_l(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
            Centred_all_n(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_all_n(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
            Centred_all_h(ishuffle,iseries,bin,:,:) = circshift(sq(Centred_all_h(ishuffle,iseries,bin,:,:)),[-bin+1 -bin+1]);
        end
    end
end
%%
% High contrast
h = figure(29);
set(h, 'position', [1716 49 837 1308]);

subplot(4,2,1)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_h(1, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_s1 = get(gca,'clim');
title('Signal: Data')

subplot(4,2,2)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_h(1, :,:,:,:),3),2))');
clim_n1 = get(gca,'clim');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Noise: Data')

subplot(4,2,3)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_h(2, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_s2 = get(gca,'clim');
title('Signal: Shuffled data')

subplot(4,2,4)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_h(2, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_n2 = get(gca,'clim');
title('Noise: Shuffled data')

subplot(4,2,5)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_h(4, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Signal: Synth no noise')
clim_s3 = get(gca,'clim');

subplot(4,2,6)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_h(4, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Noise: Synth no noise')
clim_n3 = get(gca,'clim');

subplot(4,2,7)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_h(5, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Signal: Synth with noise')
clim_s4 = get(gca,'clim');

subplot(4,2,8)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_h(5, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_n4 = get(gca,'clim');
title('Noise: Synth with noise')

s12_lims = [min([clim_s1 clim_s2]) max([clim_s1 clim_s2])];
subplot(4,2,1)
set(gca, 'clim', s12_lims);
subplot(4,2,3)
set(gca, 'clim', s12_lims);

s34_lims = [min([clim_s3 clim_s4]) max([clim_s3 clim_s4])];
subplot(4,2,5)
set(gca, 'clim', s34_lims);
subplot(4,2,7)
set(gca, 'clim', s34_lims);


n12_lims = [min([clim_n1 clim_n2]) max([clim_n1 clim_n2])];
n1234_lims = [min([clim_n1 clim_n2 clim_n3 clim_n4]) max([clim_n1 clim_n2 clim_n3 clim_n4])];
subplot(4,2,2)
set(gca, 'clim', n12_lims);
subplot(4,2,4)
set(gca, 'clim', n12_lims);

n34_lims = [min([clim_n3 clim_n4]) max([clim_n3 clim_n4])];
subplot(4,2,6)
set(gca, 'clim', n12_lims);
subplot(4,2,8)
set(gca, 'clim', n12_lims);

for n = 1:8
    subplot(4,2,n)
    axis xy; axis square; colorbar; 
    line(xlim, [0 0], 'color','k'); 
    line([0 0], ylim, 'color','k'); 
    axis([-50 50 -50 50]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
end

% medium contrast
h = figure(30);
set(h, 'position', [862 49 837 1308]);

subplot(4,2,1)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_n(1, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_s1 = get(gca,'clim');
title('Signal: Data')

subplot(4,2,2)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_n(1, :,:,:,:),3),2))');
clim_n1 = get(gca,'clim');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Noise: Data')

subplot(4,2,3)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_n(2, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_s2 = get(gca,'clim');
title('Signal: Shuffled data')

subplot(4,2,4)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_n(2, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_n2 = get(gca,'clim');
title('Noise: Shuffled data')

subplot(4,2,5)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_n(4, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Signal: Synth no noise')
clim_s3 = get(gca,'clim');

subplot(4,2,6)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_n(4, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Noise: Synth no noise')
clim_n3 = get(gca,'clim');

subplot(4,2,7)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_n(5, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Signal: Synth with noise')
clim_s4 = get(gca,'clim');

subplot(4,2,8)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_n(5, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_n4 = get(gca,'clim');
title('Noise: Synth with noise')

s12_lims = [min([clim_s1 clim_s2]) max([clim_s1 clim_s2])];
subplot(4,2,1)
set(gca, 'clim', s12_lims);
subplot(4,2,3)
set(gca, 'clim', s12_lims);

s34_lims = [min([clim_s3 clim_s4]) max([clim_s3 clim_s4])];
subplot(4,2,5)
set(gca, 'clim', s34_lims);
subplot(4,2,7)
set(gca, 'clim', s34_lims);


n12_lims = [min([clim_n1 clim_n2]) max([clim_n1 clim_n2])];
n1234_lims = [min([clim_n1 clim_n2 clim_n3 clim_n4]) max([clim_n1 clim_n2 clim_n3 clim_n4])];
subplot(4,2,2)
set(gca,'clim',[-0.002 0.003]); % set(gca, 'clim', n1234_lims);
subplot(4,2,4)
set(gca,'clim',[-0.002 0.003]); % set(gca, 'clim', n1234_lims);

n34_lims = [min([clim_n3 clim_n4]) max([clim_n3 clim_n4])];
subplot(4,2,6)
set(gca,'clim',[-0.002 0.003]); % set(gca, 'clim', n1234_lims);
subplot(4,2,8)
set(gca,'clim',[-0.002 0.003]); % set(gca, 'clim', n1234_lims);

for n = 1:8
    subplot(4,2,n)
    axis xy; axis square; colorbar; 
    line(xlim, [0 0], 'color','k'); 
    line([0 0], ylim, 'color','k'); 
    axis([-50 50 -50 50]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
end
% Low contrast
h = figure(31);
set(h, 'position', [9 49 837 1308]);

subplot(4,2,1)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_l(1, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_s1 = get(gca,'clim');
title('Signal: Data')

subplot(4,2,2)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_l(1, :,:,:,:),3),2))');
clim_n1 = get(gca,'clim');
axis xy; axis square; colorbar;
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Noise: Data')

subplot(4,2,3)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_l(2, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_s2 = get(gca,'clim');
title('Signal: Shuffled data')

subplot(4,2,4)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_l(2, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_n2 = get(gca,'clim');
title('Noise: Shuffled data')

subplot(4,2,5)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_l(4, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Signal: Synth no noise')
clim_s3 = get(gca,'clim');

subplot(4,2,6)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_l(4, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Noise: Synth no noise')
clim_n3 = get(gca,'clim');

subplot(4,2,7)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_signal_l(5, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
title('Signal: Synth with noise')
clim_s4 = get(gca,'clim');set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');

subplot(4,2,8)
imagesc(-100:2:100,-100:2:100,sq(mean(mean(Centred_noise_l(5, :,:,:,:),3),2))');
axis xy; axis square; colorbar; 
line(xlim, [0 0], 'color','k'); line([0 0], ylim, 'color','k'); axis([-50 50 -50 50]);
clim_n4 = get(gca,'clim');set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
title('Noise: Synth with noise')

s12_lims = [min([clim_s1 clim_s2]) max([clim_s1 clim_s2])];
subplot(4,2,1)
set(gca, 'clim', s12_lims);set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
subplot(4,2,3)
set(gca, 'clim', s12_lims);set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');

s34_lims = [min([clim_s3 clim_s4]) max([clim_s3 clim_s4])];
subplot(4,2,5)
set(gca, 'clim', s34_lims);set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
subplot(4,2,7)
set(gca, 'clim', s34_lims);


n12_lims = [min([clim_n1 clim_n2]) max([clim_n1 clim_n2])];
n1234_lims = [min([clim_n1 clim_n2 clim_n3 clim_n4]) max([clim_n1 clim_n2 clim_n3 clim_n4])];
subplot(4,2,2)
set(gca,'clim',[-0.002 0.003]);%set(gca,'clim',[-0.002 0.003]); % set(gca, 'clim', n1234_lims);
subplot(4,2,4)
set(gca,'clim',[-0.002 0.003]);%set(gca,'clim',[-0.002 0.003]); % set(gca, 'clim', n1234_lims);

n34_lims = [min([clim_n3 clim_n4]) max([clim_n3 clim_n4])];
subplot(4,2,6)
set(gca,'clim',[-0.002 0.003]); % set(gca, 'clim', n1234_lims);
subplot(4,2,8)
set(gca,'clim',[-0.002 0.003]); % set(gca, 'clim', n1234_lims);
for n = 1:8
    subplot(4,2,n)
    axis xy; axis square; colorbar; 
    line(xlim, [0 0], 'color','k'); 
    line([0 0], ylim, 'color','k'); 
    axis([-50 50 -50 50]);
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
end