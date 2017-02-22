function meanPostTrend = Posterior_predictions(Posterior_all, start, stop, subIdx)
%% Defining conditions

if nargin<4
    subIdx = [1 1];
    figure('Position',[323         678        1206         420]);
end
if nargin<3
    start = 25;
    stop  = 30;
end
es = Posterior_all.data;

numBins = Posterior_all.decoder.numBins;
baseline = 1./numBins;

base = es.traj~=0 & es.contrast~=0 & ~isnan(es.traj) & es.gain==1 & es.roomLength==1 & es.smthBallSpd>5 & es.trajspeed>=0;
t_low  = es.contrast==0.18 & base;
t      = es.contrast==0.6 & base;
t_high = es.contrast==0.72 & base;


%% To get the posterior in a region and licks

P_high = Posterior_all.Posterior_high;
P_norm = Posterior_all.Posterior_norm;
P_low  = Posterior_all.Posterior_low;

P_high = P_high.*((ones(size(P_high,2),1)*(max(P_high') - min(P_high')))');
P_norm = P_norm.*((ones(size(P_norm,2),1)*(max(P_norm') - min(P_norm')))');
P_low  = P_low .*((ones(size(P_low ,2),1)*(max(P_low')  - min(P_low' )))');

% meanPostTrend.P_high = P_high;
% meanPostTrend.P_low  = P_low ;
% meanPostTrend.P_norm = P_norm;
% meanPostTrend.es     = es;

trialEnds = [find(diff(es.trialID)>=1)-15];
temp = zeros(size(es.traj));
temp(trialEnds) = 1;
trialEnds = temp;
output.outcome.high= (es.outcome(trialEnds & es.contrast>0.6 & es.contrast>0));
output.outcome.low = (es.outcome(trialEnds & es.contrast<0.6 & es.contrast>0));
output.outcome.norm= (es.outcome(trialEnds & es.contrast==0.6 & es.gain==1 & es.roomLength==1));

output.outcome.complete_high= (es.trajPercent(trialEnds & es.contrast>0.6 & es.contrast>0)>80);
output.outcome.complete_low = (es.trajPercent(trialEnds & es.contrast<0.6 & es.contrast>0)>80);
output.outcome.complete_norm= (es.trajPercent(trialEnds & es.contrast==0.6 & es.gain==1 & es.roomLength==1)>80);

output.outcome.correctTrials_high = es.trialID(trialEnds & t_high   & es.outcome==2);
output.outcome.correctTrials_low  = es.trialID(trialEnds & t_low    & es.outcome==2);
output.outcome.correctTrials_norm = es.trialID(trialEnds & t        & es.outcome==2);

output.outcome.wrongTrials_high = es.trialID(trialEnds & t_high     & (es.outcome==0 & es.trajPercent<80));
output.outcome.wrongTrials_low  = es.trialID(trialEnds & t_low      & (es.outcome==0 & es.trajPercent<80));
output.outcome.wrongTrials_norm = es.trialID(trialEnds & t          & (es.outcome==0 & es.trajPercent<80));

output.outcome.missTrials_high = es.trialID(trialEnds & t_high  & es.outcome==0 & es.trajPercent>80);
output.outcome.missTrials_low  = es.trialID(trialEnds & t_low   & es.outcome==0 & es.trajPercent>80);
output.outcome.missTrials_norm = es.trialID(trialEnds & t       & es.outcome==0 & es.trajPercent>80);

%%
meanPostTrend.low.correct = nan*ones(length(output.outcome.correctTrials_low),100);
meanPostTrend.low.wrong = nan*ones(length(output.outcome.wrongTrials_low),100);
meanPostTrend.low.miss = nan*ones(length(output.outcome.missTrials_low),100);

meanPostTrend.norm.correct = nan*ones(length(output.outcome.correctTrials_norm),100);
meanPostTrend.norm.wrong = nan*ones(length(output.outcome.wrongTrials_norm),100);
meanPostTrend.norm.miss = nan*ones(length(output.outcome.missTrials_norm),100);

meanPostTrend.high.correct = nan*ones(length(output.outcome.correctTrials_high ),100);
meanPostTrend.high.wrong = nan*ones(length(output.outcome.wrongTrials_high),100);
meanPostTrend.high.miss = nan*ones(length(output.outcome.missTrials_high),100);

temp_low = zeros(size(t));
meanPostTrend.low.correctIdx = zeros(size(t));
meanPostTrend.low.wrongIdx = zeros(size(t));
meanPostTrend.low.missIdx = zeros(size(t));

meanPostTrend.licks.traj = [];
meanPostTrend.licks.meanPost = [];
meanPostTrend.licks.trialType = [];

subplot(subIdx(1),3,3*(subIdx(2)-1)+1)
for tIdx = 1:length(output.outcome.correctTrials_low)'
    meanPostTrend.low.correct(tIdx,:) = ...
        interp1q(es.traj(t_low & es.trialID==output.outcome.correctTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.correctTrials_low(tIdx),start:stop),2), (1:100)');
    
    trajs = es.traj(t_low & es.trialID==output.outcome.correctTrials_low(tIdx) & es.lick);
    meanPosts = nanmean(P_low(es.trialID(t_low & es.lick)==output.outcome.correctTrials_low(tIdx),start:stop),2);
    trialType = ones(sum(t_low & es.trialID==output.outcome.correctTrials_low(tIdx) & es.lick),1)*1;
    meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
    meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
    meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
    
    plot(es.traj(t_low & es.trialID==output.outcome.correctTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.correctTrials_low(tIdx),start:stop),2),'color',[.5 .5 .5]);
    hold on;
end
for tIdx = 1:length(output.outcome.missTrials_low)'
    meanPostTrend.low.miss(tIdx,:) = ...
        interp1q(es.traj(t_low & es.trialID==output.outcome.missTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.missTrials_low(tIdx),start:stop),2), (1:100)');
    trajs = es.traj(t_low & es.trialID==output.outcome.missTrials_low(tIdx) & es.lick);
    
    meanPosts = nanmean(P_low(es.trialID(t_low & es.lick)==output.outcome.missTrials_low(tIdx),start:stop),2);
    trialType = ones(sum(t_low & es.trialID==output.outcome.missTrials_low(tIdx) & es.lick),1)*0;
    meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
    meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
    meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
    
    plot(es.traj(t_low & es.trialID==output.outcome.missTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.missTrials_low(tIdx),start:stop),2),'color',[0 0 1]);
end
for tIdx = 1:length(output.outcome.wrongTrials_low)'
    meanPostTrend.low.wrong(tIdx,:) = ...
        interp1q(es.traj(t_low & es.trialID==output.outcome.wrongTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.wrongTrials_low(tIdx),start:stop),2), (1:100)');
    trajs = es.traj(t_low & es.trialID==output.outcome.wrongTrials_low(tIdx) & es.lick);
    
    meanPosts = nanmean(P_low(es.trialID(t_low & es.lick)==output.outcome.wrongTrials_low(tIdx),start:stop),2);
    trialType = ones(sum(t_low & es.trialID==output.outcome.wrongTrials_low(tIdx) & es.lick),1)*2;
    meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
    meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
    meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
    
    plot(es.traj(t_low & es.trialID==output.outcome.wrongTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.wrongTrials_low(tIdx),start:stop),2),'color',[1 0 0]);
    %     2*Posterior_all.X_low(es.trialID(t_low)==tIdx),
end
ylims_low = ylim;
meanPostTrend.all.correct = meanPostTrend.low.correct;
meanPostTrend.all.miss = meanPostTrend.low.miss;
meanPostTrend.all.wrong = meanPostTrend.low.wrong;

try
    errorbar(1:100, nanmean(meanPostTrend.low.correct,1), nansem(meanPostTrend.low.correct),'k')
catch;end
try
    errorbar(1:100, nanmean(meanPostTrend.low.miss,1), nansem(meanPostTrend.low.miss),'b')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.low.wrong,1), nansem(meanPostTrend.low.wrong),'r')
catch; end
title('Low Contrast')

subplot(subIdx(1),3,3*(subIdx(2)-1)+2)
for tIdx = 1:length(output.outcome.correctTrials_norm)'
    try
        meanPostTrend.low.correct(tIdx,:) = ...
            interp1q(es.traj(t & es.trialID==output.outcome.correctTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.correctTrials_norm(tIdx),start:stop),2), (1:100)');
        
        trajs = es.traj(t & es.trialID==output.outcome.correctTrials_norm(tIdx) & es.lick);
        meanPosts = nanmean(P_norm(es.trialID(t & es.lick)==output.outcome.correctTrials_norm(tIdx),start:stop),2);
        trialType = ones(sum(t & es.trialID==output.outcome.correctTrials_norm(tIdx) & es.lick),1)*1;
        meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
        meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
        meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
        
        plot(es.traj(t & es.trialID==output.outcome.correctTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.correctTrials_norm(tIdx),start:stop),2),'color',[.5 .5 .5]);
        hold on;
    catch; end
end
for tIdx = 1:length(output.outcome.missTrials_norm)'
    try
        meanPostTrend.low.miss(tIdx,:) = ...
            interp1q(es.traj(t & es.trialID==output.outcome.missTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.missTrials_norm(tIdx),start:stop),2), (1:100)');
        trajs = es.traj(t & es.trialID==output.outcome.missTrials_norm(tIdx) & es.lick);
        
        meanPosts = nanmean(P_norm(es.trialID(t & es.lick)==output.outcome.missTrials_norm(tIdx),start:stop),2);
        trialType = ones(sum(t & es.trialID==output.outcome.missTrials_norm(tIdx) & es.lick),1)*0;
        meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
        meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
        meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
        
        plot(es.traj(t & es.trialID==output.outcome.missTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.missTrials_norm(tIdx),start:stop),2),'color',[0 0 1]);
    catch; end
end
for tIdx = 1:length(output.outcome.wrongTrials_norm)'
    try
        meanPostTrend.low.wrong(tIdx,:) = ...
            interp1q(es.traj(t & es.trialID==output.outcome.wrongTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.wrongTrials_norm(tIdx),start:stop),2), (1:100)');
        trajs = es.traj(t & es.trialID==output.outcome.wrongTrials_norm(tIdx) & es.lick);
        
        meanPosts = nanmean(P_norm(es.trialID(t & es.lick)==output.outcome.wrongTrials_norm(tIdx),start:stop),2);
        trialType = ones(sum(t & es.trialID==output.outcome.wrongTrials_norm(tIdx) & es.lick),1)*2;
        meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
        meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
        meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
        
        plot(es.traj(t & es.trialID==output.outcome.wrongTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.wrongTrials_norm(tIdx),start:stop),2),'color',[1 0 0]);
        %     2*Posterior_all.X_norm(es.trialID(t)==tIdx),
    catch; end
end

meanPostTrend.all.correct = [meanPostTrend.all.correct' meanPostTrend.norm.correct']';
meanPostTrend.all.miss =    [meanPostTrend.all.miss' meanPostTrend.norm.miss']';
meanPostTrend.all.wrong =   [meanPostTrend.all.wrong' meanPostTrend.norm.wrong']';

ylims_norm = ylim;
title('Baseline')
try
    errorbar(1:100, nanmean(meanPostTrend.norm.correct,1), nansem(meanPostTrend.norm.correct),'k')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.norm.miss,1), nansem(meanPostTrend.norm.miss),'b')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.norm.wrong,1), nansem(meanPostTrend.norm.wrong),'r')
catch; end

subplot(subIdx(1),3,3*(subIdx(2)-1)+3)
for tIdx = 1:length(output.outcome.correctTrials_high)'
    meanPostTrend.low.correct(tIdx,:) = ...
        interp1q(es.traj(t_high & es.trialID==output.outcome.correctTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.correctTrials_high(tIdx),start:stop),2), (1:100)');
    
    trajs = es.traj(t_high & es.trialID==output.outcome.correctTrials_high(tIdx) & es.lick);
    meanPosts = nanmean(P_high(es.trialID(t_high & es.lick)==output.outcome.correctTrials_high(tIdx),start:stop),2);
    trialType = ones(sum(t_high & es.trialID==output.outcome.correctTrials_high(tIdx) & es.lick),1)*1;
    meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
    meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
    meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
    
    plot(es.traj(t_high & es.trialID==output.outcome.correctTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.correctTrials_high(tIdx),start:stop),2),'color',[.5 .5 .5]);
    hold on;
end
for tIdx = 1:length(output.outcome.missTrials_high)'
    meanPostTrend.low.miss(tIdx,:) = ...
        interp1q(es.traj(t_high & es.trialID==output.outcome.missTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.missTrials_high(tIdx),start:stop),2), (1:100)');
    trajs = es.traj(t_high & es.trialID==output.outcome.missTrials_high(tIdx) & es.lick);
    
    meanPosts = nanmean(P_high(es.trialID(t_high & es.lick)==output.outcome.missTrials_high(tIdx),start:stop),2);
    trialType = ones(sum(t_high & es.trialID==output.outcome.missTrials_high(tIdx) & es.lick),1)*0;
    meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
    meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
    meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
    
    plot(es.traj(t_high & es.trialID==output.outcome.missTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.missTrials_high(tIdx),start:stop),2),'color',[0 0 1]);
end
for tIdx = 1:length(output.outcome.wrongTrials_high)'
    meanPostTrend.low.wrong(tIdx,:) = ...
        interp1q(es.traj(t_high & es.trialID==output.outcome.wrongTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.wrongTrials_high(tIdx),start:stop),2), (1:100)');
    trajs = es.traj(t_high & es.trialID==output.outcome.wrongTrials_high(tIdx) & es.lick);
    
    meanPosts = nanmean(P_high(es.trialID(t_high & es.lick)==output.outcome.wrongTrials_high(tIdx),start:stop),2);
    trialType = ones(sum(t_high & es.trialID==output.outcome.wrongTrials_high(tIdx) & es.lick),1)*2;
    meanPostTrend.licks.traj        = [meanPostTrend.licks.traj trajs'];
    meanPostTrend.licks.meanPost    = [meanPostTrend.licks.meanPost meanPosts'];
    meanPostTrend.licks.trialType   = [meanPostTrend.licks.trialType trialType'];
    
    plot(es.traj(t_high & es.trialID==output.outcome.wrongTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.wrongTrials_high(tIdx),start:stop),2),'color',[1 0 0]);
    %     2*Posterior_all.X_high(es.trialID(t_high)==tIdx),
end
meanPostTrend.all.correct = [meanPostTrend.all.correct' meanPostTrend.high.correct']';
meanPostTrend.all.miss =    [meanPostTrend.all.miss' meanPostTrend.high.miss']';
meanPostTrend.all.wrong =   [meanPostTrend.all.wrong' meanPostTrend.high.wrong']';

errorbar(1:100, nanmean(meanPostTrend.high.correct,1), nansem(meanPostTrend.high.correct),'k')
try
    errorbar(1:100, nanmean(meanPostTrend.high.miss,1), nansem(meanPostTrend.high.miss),'b')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.high.wrong,1), nansem(meanPostTrend.high.wrong),'r')
catch; end

title('High contrast')
ylims_high = ylim;
ylims = [min([ylims_low ylims_norm ylims_high]) max([ylims_low ylims_norm ylims_high])];
for pIdx = 3:-1:1
    subplot(subIdx(1),3,3*(subIdx(2)-1)+pIdx)
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none','YLim',ylims);
    xlabel('Position (cm)')
    line([70-5 70-5],ylim,'color','c');
    line([70+5 70+5],ylim,'color','c');
end
ylabel(['Posterior in ' num2str(2*start) ' to ' num2str(2*stop)]);

norm_factor = nanmean(meanPostTrend.all.correct(:,70),1);
meanPostTrend.all.correct = meanPostTrend.all.correct./norm_factor;
meanPostTrend.all.miss = meanPostTrend.all.miss./norm_factor;
meanPostTrend.all.wrong = meanPostTrend.all.wrong./norm_factor;
meanPostTrend.norm_factor = norm_factor;
figure;
try
    errorbar(1:100, nanmean(meanPostTrend.all.correct,1), nansem(meanPostTrend.all.correct),'k');
    hold on;
    errorbar(1:100, nanmean(meanPostTrend.all.miss,1), nansem(meanPostTrend.all.miss),'b')
    errorbar(1:100, nanmean(meanPostTrend.all.wrong,1), nansem(meanPostTrend.all.wrong),'r')
catch; end;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlabel('Position (cm)')
ylabel(['Posterior in ' num2str(2*start) ' to ' num2str(2*stop)]);
line([70-5 70-5],ylim,'color','c');
line([70+5 70+5],ylim,'color','c');

    function output = getsumPost(X)
        Y = (2.^X)*baseline;
        output = log2(nansum(Y,2)./baseline);
    end
end