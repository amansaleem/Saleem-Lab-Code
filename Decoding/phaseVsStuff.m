function Posterior_all = phaseVsStuff(Posterior_all, stop_point, start_point, spread)
if nargin<2
    stop_point = 65;
end
if nargin<3
    start_point = 0;
end
if nargin<4
    spread = 35;
end

% Posterior_all.thetaBins = Posterior_all.thetaBins*(180/pi)+ 180;
Posterior_all.meanPost_diff_norm = zeros(length(Posterior_all.thetaBins),99);
Posterior_all.meanPost_diff_low  = zeros(length(Posterior_all.thetaBins),99);
Posterior_all.meanPost_diff_high = zeros(length(Posterior_all.thetaBins),99);

%% Getting the mean error for each phase
for idx = 1:length(Posterior_all.thetaBins)
    Posterior_all.meanErrorforPhase_norm(idx) = nanmean((Posterior_all.error.norm(Posterior_all.thetaPhase(Posterior_all.t_norm)==idx & ...
        Posterior_all.data.outcome(Posterior_all.t_norm) == 2 & ... 
        (Posterior_all.data.traj(Posterior_all.t_norm)>=start_point & Posterior_all.data.traj(Posterior_all.t_norm)<=stop_point))));
    Posterior_all.meanErrorforPhase_low(idx)  = nanmean((Posterior_all.error.low(Posterior_all.thetaPhase(Posterior_all.t_low)==idx & ...
        Posterior_all.data.outcome(Posterior_all.t_low) == 2 & ...
        (Posterior_all.data.traj(Posterior_all.t_low)>=start_point & Posterior_all.data.traj(Posterior_all.t_low)<=stop_point))));
    Posterior_all.meanErrorforPhase_high(idx) = nanmean((Posterior_all.error.high(Posterior_all.thetaPhase(Posterior_all.t_high)==idx & ...
        Posterior_all.data.outcome(Posterior_all.t_high) == 2 & ...
        (Posterior_all.data.traj(Posterior_all.t_high)>=start_point & Posterior_all.data.traj(Posterior_all.t_high)<=stop_point))));
    
    Posterior_all.meanConfForPhase_norm(idx) = nanmean((Posterior_all.confidence.norm(Posterior_all.thetaPhase(Posterior_all.t_norm)==idx & ...
        (Posterior_all.data.traj(Posterior_all.t_norm)>=start_point & Posterior_all.data.traj(Posterior_all.t_norm)<=stop_point))));
    Posterior_all.meanConfForPhase_low(idx)  = nanmean((Posterior_all.confidence.low(Posterior_all.thetaPhase(Posterior_all.t_low)==idx & ...
        (Posterior_all.data.traj(Posterior_all.t_low)>=start_point & Posterior_all.data.traj(Posterior_all.t_low)<=stop_point))));
    Posterior_all.meanConfForPhase_high(idx) = nanmean((Posterior_all.confidence.high(Posterior_all.thetaPhase(Posterior_all.t_high)==idx & ...
        (Posterior_all.data.traj(Posterior_all.t_high)>=start_point & Posterior_all.data.traj(Posterior_all.t_high)<=stop_point))));
    
    Posterior_all.meanPost_diff_norm(idx,:) = nanmean(Posterior_all.Posterior_diff_norm(Posterior_all.thetaPhase(Posterior_all.t_norm)==idx & ...
        Posterior_all.data.outcome(Posterior_all.t_norm) == 2 & ...
        (Posterior_all.data.traj(Posterior_all.t_norm)>=start_point & Posterior_all.data.traj(Posterior_all.t_norm)<=stop_point),:),1);
    Posterior_all.meanPost_diff_low(idx,:)  = nanmean(Posterior_all.Posterior_diff_low(Posterior_all.thetaPhase(Posterior_all.t_low)==idx & ...
        Posterior_all.data.outcome(Posterior_all.t_low) == 2 & ...
        (Posterior_all.data.traj(Posterior_all.t_low)>=start_point & Posterior_all.data.traj(Posterior_all.t_low)<=stop_point),:),1);
    Posterior_all.meanPost_diff_high(idx,:) = nanmean(Posterior_all.Posterior_diff_high(Posterior_all.thetaPhase(Posterior_all.t_high)==idx & ...
        Posterior_all.data.outcome(Posterior_all.t_high) == 2 & ...
        (Posterior_all.data.traj(Posterior_all.t_high)>=start_point & Posterior_all.data.traj(Posterior_all.t_high)<=stop_point),:),1);
end

%% Plotting stuff

% Plotting the mean error vs. phase
figure('Position',[267 126 1200 800]);
subplot(2,2,1)
plot( Posterior_all.meanErrorforPhase_norm,Posterior_all.thetaBins,'.-k')
hold on;
plot( Posterior_all.meanErrorforPhase_norm,360+Posterior_all.thetaBins,'.-k')

plot( Posterior_all.meanErrorforPhase_low,Posterior_all.thetaBins,'.-b')
plot( Posterior_all.meanErrorforPhase_low,360+Posterior_all.thetaBins,'.-b')

plot( Posterior_all.meanErrorforPhase_high,Posterior_all.thetaBins,'.-r')
plot( Posterior_all.meanErrorforPhase_high,360+Posterior_all.thetaBins,'.-r')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');

% subplot(2,2,3)
% imagesc(Posterior_all.thetaBins, (-10:10), Posterior_all.meanPost_diff_norm(:,41:60)')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% % set(gca,'Clim',[-0.25 0.25])
% RedWhiteBlue; axis xy
% grid on
% 
% subplot(2,2,2)
% imagesc(Posterior_all.thetaBins, (-10:10), Posterior_all.meanPost_diff_low(:,41:60)')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% % set(gca,'Clim',[-0.25 0.25])
% RedWhiteBlue; axis xy;
% grid on
% 
% subplot(2,2,4)
% imagesc(Posterior_all.thetaBins, (-10:10), Posterior_all.meanPost_diff_high(:,41:60)')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% % set(gca,'Clim',[-0.25 0.25])
% RedWhiteBlue; axis xy;
% grid on

%% Plotting the difference in posterior vs. theta phase
subplot(2,2,3)
imagesc((-spread:spread),[Posterior_all.thetaBins 360+ Posterior_all.thetaBins],...
    [Posterior_all.meanPost_diff_norm(:,(49-spread):(49+spread))' Posterior_all.meanPost_diff_norm(:,(49-spread):(49+spread))']')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% set(gca,'Clim',[-0.25 0.25])
RedWhiteBlue; axis xy; colorbar
title('Normal Contrast')
grid on

subplot(2,2,2)
imagesc((-spread:spread),[Posterior_all.thetaBins 360+ Posterior_all.thetaBins],...
    [Posterior_all.meanPost_diff_low(:,(49-spread):(49+spread))' Posterior_all.meanPost_diff_low(:,(49-spread):(49+spread))']')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% set(gca,'Clim',[-0.25 0.25])
RedWhiteBlue; axis xy; colorbar
title('Low Contrast')
grid on

subplot(2,2,4)
imagesc((-spread:spread),[Posterior_all.thetaBins 360+ Posterior_all.thetaBins],...
    [Posterior_all.meanPost_diff_high(:,(49-spread):(49+spread))' Posterior_all.meanPost_diff_high(:,(49-spread):(49+spread))']')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% set(gca,'Clim',[-0.25 0.25])
RedWhiteBlue; axis xy; colorbar
title('High Contrast')
grid on

for m=2:4; subplot(2,2,m); tmp = get(gca, 'CLim'); if tmp(1)~=tmp(2); set(gca, 'CLim', [nanmax([0 (tmp(2)-0.4)]) tmp(2)]); end; end
%%
% figure
% subplot(122)
% plot((-49:49),nanmean(Posterior_all.meanPost_diff_high),'r')
% hold on;
% plot((-49:49),nanmean(Posterior_all.meanPost_diff_low),'b')
% plot((-49:49),nanmean(Posterior_all.meanPost_diff_norm),'k')
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% % set(gca,'YLim',[-0.5 0.5]);
% subplot(121)
% imagesc(Posterior_all.meanPost.norm); 
% axis xy; 
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% set(gca, 'CLim',[-0.5 0.5]);
% RedWhiteBlue;
% axis equal; 
% line([1 50], [1 50], 'linestyle','--', 'color', 'k');
% axis off;
end


