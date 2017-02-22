function meanPostTrend = Posterior_predictions_safe(Posterior_all, start, stop, subIdx)
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
% t_low  = es.contrast==0.18 & base;
% t      = es.contrast==0.6 & base;
% t_high = es.contrast==0.72 & base;
t_low = Posterior_all.t_low;
t = Posterior_all.t_norm;
t_high = Posterior_all.t_high;

es_new.traj = es.traj;
es_new.lick = es.lick;
es_new.t_low = t_low;
es_new.t_norm = t;
es_new.t_high = t_high;
es_new.outcome = nan*ones(size(es_new.traj));
es_new.posterior = nan*ones(size(es_new.traj));
es_new.X_ML = nan*ones(size(es_new.traj));

diffs.miss = [];
diffs.correct = [];
diffs.wrong = [];

%% To get the posterior in a region and licks

P_norm = Posterior_all.Posterior_norm;

% if nanmin(P_norm)>-.5
%     P_high = log2(Posterior_all.Posterior_high*50);
%     P_norm = log2(Posterior_all.Posterior_norm*50);
%     P_low  = log2(Posterior_all.Posterior_low*50);
% else
    P_high = Posterior_all.Posterior_high;
    P_norm = Posterior_all.Posterior_norm;
    P_low  = Posterior_all.Posterior_low;
% end
if sum(t)>size(P_norm,1)
    temp = find(t);
    t(temp(size(P_norm,1)+1:end)) = 0;
end
es_new.posterior(t_low) = nanmean(P_low(:,start:stop),2);
es_new.posterior(t) = nanmean(P_norm(:,start:stop),2);
es_new.posterior(t_high) = nanmean(P_high(:,start:stop),2);

[~,es_new.X_ML(t_low)]  = max(Posterior_all.Posterior_low');
[~,es_new.X_ML(t)]      = max(Posterior_all.Posterior_norm');
[~,es_new.X_ML(t_high)] = max(Posterior_all.Posterior_high');

% P_high = P_high.*((ones(size(P_high,2),1)*(max(P_high') - min(P_high')))');
% P_norm = P_norm.*((ones(size(P_norm,2),1)*(max(P_norm') - min(P_norm')))');
% P_low  = P_low .*((ones(size(P_low ,2),1)*(max(P_low')  - min(P_low' )))');

% meanPostTrend.P_high = P_high;
% meanPostTrend.P_low  = P_low ;
% meanPostTrend.P_norm = P_norm;
% meanPostTrend.es     = es;

%% redefining the outcome
trialEnds = [find(diff(es.trialID)>=1)-15];
temp = zeros(size(es.traj));
temp(trialEnds) = 1;
trialEnds = temp;
output.outcome.high= Posterior_all.outcome.highContrast;% (es.outcome(trialEnds & es.contrast>0.6 & es.contrast>0));
output.outcome.low = Posterior_all.outcome.lowContrast;%(es.outcome(trialEnds & es.contrast<0.6 & es.contrast>0));
output.outcome.norm= Posterior_all.outcome.normContrast;%(es.outcome(trialEnds & es.contrast==0.6 & es.gain==1 & es.roomLength==1));

output.outcome.complete_high= intersect(Posterior_all.outcome.highContrast,Posterior_all.outcome.complete);%(es.trajPercent(trialEnds & es.contrast>0.6 & es.contrast>0)>80);
output.outcome.complete_low = intersect(Posterior_all.outcome.lowContrast,Posterior_all.outcome.complete);%(es.trajPercent(trialEnds & es.contrast<0.6 & es.contrast>0)>80);
output.outcome.complete_norm= intersect(Posterior_all.outcome.normContrast,Posterior_all.outcome.complete);%(es.trajPercent(trialEnds & es.contrast==0.6 & es.gain==1 & es.roomLength==1)>80);

output.outcome.correctTrials_high = intersect(Posterior_all.outcome.highContrast,Posterior_all.outcome.correctTrials); %es.trialID(trialEnds & t_high   & es.outcome==2);
output.outcome.correctTrials_low  = intersect(Posterior_all.outcome.lowContrast,Posterior_all.outcome.correctTrials); %es.trialID(trialEnds & t_low    & es.outcome==2);
output.outcome.correctTrials_norm = intersect(Posterior_all.outcome.normContrast,Posterior_all.outcome.correctTrials); % es.trialID(trialEnds & t        & es.outcome==2);

output.outcome.wrongTrials_high = intersect(Posterior_all.outcome.highContrast,Posterior_all.outcome.earlyTrials);%es.trialID(trialEnds & t_high     & (es.outcome==0 & es.trajPercent<80));
output.outcome.wrongTrials_low  = intersect(Posterior_all.outcome.lowContrast,Posterior_all.outcome.earlyTrials);%es.trialID(trialEnds & t_low      & (es.outcome==0 & es.trajPercent<80));
output.outcome.wrongTrials_norm = intersect(Posterior_all.outcome.normContrast,Posterior_all.outcome.earlyTrials);%es.trialID(trialEnds & t          & (es.outcome==0 & es.trajPercent<80));

output.outcome.lateTrials_high = intersect(Posterior_all.outcome.highContrast,Posterior_all.outcome.lateTrials);%es.trialID(trialEnds & t_high  & es.outcome==0 & es.trajPercent>80);
output.outcome.lateTrials_low  = intersect(Posterior_all.outcome.lowContrast,Posterior_all.outcome.lateTrials);%es.trialID(trialEnds & t_low   & es.outcome==0 & es.trajPercent>80);
output.outcome.lateTrials_norm = intersect(Posterior_all.outcome.normContrast,Posterior_all.outcome.lateTrials);%es.trialID(trialEnds & t       & es.outcome==0 & es.trajPercent>80);

output.outcome.missTrials_high = intersect(Posterior_all.outcome.highContrast,Posterior_all.outcome.missTrials);%es.trialID(trialEnds & t_high  & es.outcome==0 & es.trajPercent>80);
output.outcome.missTrials_low  = intersect(Posterior_all.outcome.lowContrast,Posterior_all.outcome.missTrials);%es.trialID(trialEnds & t_low   & es.outcome==0 & es.trajPercent>80);
output.outcome.missTrials_norm = intersect(Posterior_all.outcome.normContrast,Posterior_all.outcome.missTrials);%es.trialID(trialEnds & t       & es.outcome==0 & es.trajPercent>80);

for tIdx = [output.outcome.correctTrials_high' output.outcome.correctTrials_norm' output.outcome.correctTrials_low']
    es_new.outcome(es.trialID==tIdx) = 1;
end
for tIdx = [output.outcome.wrongTrials_high' output.outcome.wrongTrials_norm' output.outcome.wrongTrials_low']
    es_new.outcome(es.trialID==tIdx)   = 2;
end
try
for tIdx = [output.outcome.lateTrials_high' output.outcome.lateTrials_norm' output.outcome.lateTrials_low']
    es_new.outcome(es.trialID==tIdx)    = 3;
end
catch; end
for tIdx = [output.outcome.missTrials_high' output.outcome.missTrials_norm' output.outcome.missTrials_low']
    es_new.outcome(es.trialID==tIdx)    = 4;
end
%%
meanPostTrend.low.correct = nan*ones(length(output.outcome.correctTrials_low),100);
meanPostTrend.low.wrong = nan*ones(length(output.outcome.wrongTrials_low),100);
meanPostTrend.low.miss = nan*ones(length(output.outcome.missTrials_low),100);
meanPostTrend.low.late = nan*ones(length(output.outcome.lateTrials_low),100);

meanPostTrend.norm.correct = nan*ones(length(output.outcome.correctTrials_norm),100);
meanPostTrend.norm.wrong = nan*ones(length(output.outcome.wrongTrials_norm),100);
meanPostTrend.norm.miss = nan*ones(length(output.outcome.missTrials_norm),100);
meanPostTrend.norm.late = nan*ones(length(output.outcome.lateTrials_norm),100);

meanPostTrend.high.correct = nan*ones(length(output.outcome.correctTrials_high ),100);
meanPostTrend.high.wrong = nan*ones(length(output.outcome.wrongTrials_high),100);
meanPostTrend.high.miss = nan*ones(length(output.outcome.missTrials_high),100);
meanPostTrend.high.late = nan*ones(length(output.outcome.lateTrials_high),100);

temp_low = zeros(size(t));
meanPostTrend.low.correctIdx = zeros(size(t));
meanPostTrend.low.wrongIdx = zeros(size(t));
meanPostTrend.low.missIdx = zeros(size(t));
meanPostTrend.low.lateIdx = zeros(size(t));

meanPostTrend.licks.traj = [];
meanPostTrend.licks.meanPost = [];
meanPostTrend.licks.trailType = [];

subplot(subIdx(1),3,3*(subIdx(2)-1)+1)
for tIdx = 1:length(output.outcome.correctTrials_low)'
    meanPostTrend.low.correct(tIdx,:) = ...
        interp1q(es.traj(t_low & es.trialID==output.outcome.correctTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.correctTrials_low(tIdx),start:stop),2), (1:100)');
    trajs = es.traj(t_low & es.trialID==output.outcome.correctTrials_low(tIdx) & es.lick);
    meanPosts = es.traj(t_low & es.trialID==output.outcome.correctTrials_low(tIdx) & es.lick);
   
%     diffs.correct = [diffs.correct diff(Posterior_all.MAP.low(es.trialID(t_low)==output.outcome.correctTrials_low(tIdx) & (Posterior_all.X_low<start)))];
    plot(es.traj(t_low & es.trialID==output.outcome.correctTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.correctTrials_low(tIdx),start:stop),2),'color',[.5 .5 .5]);
    hold on;
end
for tIdx = 1:length(output.outcome.missTrials_low)'
%     diffs.miss = [diffs.miss diff(Posterior_all.MAP.low(es.trialID(t_low)==output.outcome.missTrials_low(tIdx) & (Posterior_all.X_low<start)))];
    meanPostTrend.low.miss(tIdx,:) = ...
        interp1q(es.traj(t_low & es.trialID==output.outcome.missTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.missTrials_low(tIdx),start:stop),2), (1:100)');
    plot(es.traj(t_low & es.trialID==output.outcome.missTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.missTrials_low(tIdx),start:stop),2),'color',[0 0 1]);
end
for tIdx = 1:length(output.outcome.lateTrials_low)'
%     diffs.miss = [diffs.miss diff(Posterior_all.MAP.low(es.trialID(t_low)==output.outcome.missTrials_low(tIdx) & (Posterior_all.X_low<start)))];
    meanPostTrend.low.late(tIdx,:) = ...
        interp1q(es.traj(t_low & es.trialID==output.outcome.lateTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.lateTrials_low(tIdx),start:stop),2), (1:100)');
    plot(es.traj(t_low & es.trialID==output.outcome.lateTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.lateTrials_low(tIdx),start:stop),2),'color',[0 0 1]);
end
for tIdx = 1:length(output.outcome.wrongTrials_low)'
%     diffs.wrong = [diffs.wrong diff(Posterior_all.MAP.low(es.trialID(t_low)==output.outcome.wrongTrials_low(tIdx) & (Posterior_all.X_low<start)))];
    meanPostTrend.low.wrong(tIdx,:) = ...
        interp1q(es.traj(t_low & es.trialID==output.outcome.wrongTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.wrongTrials_low(tIdx),start:stop),2), (1:100)');
    plot(es.traj(t_low & es.trialID==output.outcome.wrongTrials_low(tIdx)), nanmean(P_low(es.trialID(t_low)==output.outcome.wrongTrials_low(tIdx),start:stop),2),'color',[1 0 0]);
    %     2*Posterior_all.X_low(es.trialID(t_low)==tIdx),
end
ylims_low = ylim;
meanPostTrend.all.correct = meanPostTrend.low.correct;
meanPostTrend.all.miss = meanPostTrend.low.miss;
meanPostTrend.all.late = meanPostTrend.low.late;
meanPostTrend.all.wrong = meanPostTrend.low.wrong;

try
    errorbar(1:100, nanmean(meanPostTrend.low.correct,1), nansem(meanPostTrend.low.correct),'k')
catch;end
try
    errorbar(1:100, nanmean(meanPostTrend.low.miss,1), nansem(meanPostTrend.low.miss),'b')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.low.late,1), nansem(meanPostTrend.low.late),'c')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.low.wrong,1), nansem(meanPostTrend.low.wrong),'r')
catch; end
title('Low Contrast')

subplot(subIdx(1),3,3*(subIdx(2)-1)+2)
for tIdx = 1:length(output.outcome.correctTrials_norm)'
    try
        meanPostTrend.norm.correct(tIdx,:) = ...
            interp1q(es.traj(t & es.trialID==output.outcome.correctTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.correctTrials_norm(tIdx),start:stop),2), (1:100)');
        plot(es.traj(t & es.trialID==output.outcome.correctTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.correctTrials_norm(tIdx),start:stop),2),'color',[.5 .5 .5]);
%         diffs.correct = [diffs.correct diff(Posterior_all.MAP.norm(es.trialID(t)==output.outcome.correctTrials_low(tIdx) & (Posterior_all.X_norm<start)))];
        hold on;
    catch
        meanPostTrend.norm.correct(tIdx,:) = 0;
    end
end
for tIdx = 1:length(output.outcome.missTrials_norm)'
    try
        meanPostTrend.norm.miss(tIdx,:) = ...
            interp1q(es.traj(t & es.trialID==output.outcome.missTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.missTrials_norm(tIdx),start:stop),2), (1:100)');
%         diffs.miss = [diffs.miss diff(Posterior_all.MAP.norm(es.trialID(t)==output.outcome.missTrials_low(tIdx) & (Posterior_all.X_norm<start)))];
    plot(es.traj(t & es.trialID==output.outcome.missTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.missTrials_norm(tIdx),start:stop),2),'color',[0 0 1]);
    catch
        meanPostTrend.norm.miss(tIdx,:) = 0;
    end
end

for tIdx = 1:length(output.outcome.lateTrials_norm)'
    try
        meanPostTrend.norm.late(tIdx,:) = ...
            interp1q(es.traj(t & es.trialID==output.outcome.lateTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.lateTrials_norm(tIdx),start:stop),2), (1:100)');
%         diffs.late = [diffs.late diff(Posterior_all.MAP.norm(es.trialID(t)==output.outcome.lateTrials_low(tIdx) & (Posterior_all.X_norm<start)))];
    plot(es.traj(t & es.trialID==output.outcome.lateTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.lateTrials_norm(tIdx),start:stop),2),'color',[0 0 1]);
    catch
        meanPostTrend.norm.late(tIdx,:) = 0;
    end
end
for tIdx = 1:length(output.outcome.wrongTrials_norm)'
    try
        meanPostTrend.norm.wrong(tIdx,:) = ...
            interp1q(es.traj(t & es.trialID==output.outcome.wrongTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.wrongTrials_norm(tIdx),start:stop),2), (1:100)');
%         diffs.wrong = [diffs.wrong diff(Posterior_all.MAP.norm(es.trialID(t)==output.outcome.wrongTrials_low(tIdx) & (Posterior_all.X_norm<start)))];
        plot(es.traj(t & es.trialID==output.outcome.wrongTrials_norm(tIdx)), nanmean(P_norm(es.trialID(t)==output.outcome.wrongTrials_norm(tIdx),start:stop),2),'color',[1 0 0]);
    catch
        meanPostTrend.norm.wrong(tIdx,:) = 0;
    end
    %     2*Posterior_all.X_norm(es.trialID(t)==tIdx),
end

meanPostTrend.all.correct = [meanPostTrend.all.correct' meanPostTrend.norm.correct']';
meanPostTrend.all.miss =    [meanPostTrend.all.miss' meanPostTrend.norm.miss']';
meanPostTrend.all.late =    [meanPostTrend.all.late' meanPostTrend.norm.late']';
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
    errorbar(1:100, nanmean(meanPostTrend.norm.late,1), nansem(meanPostTrend.norm.late),'c')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.norm.wrong,1), nansem(meanPostTrend.norm.wrong),'r')
catch; end

subplot(subIdx(1),3,3*(subIdx(2)-1)+3)
for tIdx = 1:length(output.outcome.correctTrials_high)'
    meanPostTrend.high.correct(tIdx,:) = ...
        interp1q(es.traj(t_high & es.trialID==output.outcome.correctTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.correctTrials_high(tIdx),start:stop),2), (1:100)');
    plot(es.traj(t_high & es.trialID==output.outcome.correctTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.correctTrials_high(tIdx),start:stop),2),'color',[.5 .5 .5]);
%     diffs.correct = [diffs.correct diff(Posterior_all.MAP.high(es.trialID(t_high)==output.outcome.correctTrials_high(tIdx) & (Posterior_all.X_high<start)))];
    hold on;
end
for tIdx = 1:length(output.outcome.missTrials_high)'
    meanPostTrend.high.miss(tIdx,:) = ...
        interp1q(es.traj(t_high & es.trialID==output.outcome.missTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.missTrials_high(tIdx),start:stop),2), (1:100)');
    plot(es.traj(t_high & es.trialID==output.outcome.missTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.missTrials_high(tIdx),start:stop),2),'color',[0 0 1]);
%     diffs.miss = [diffs.miss diff(Posterior_all.MAP.high(es.trialID(t_high)==output.outcome.missTrials_high(tIdx) & (Posterior_all.X_high<start)))];
        
end
for tIdx = 1:length(output.outcome.lateTrials_high)'
    meanPostTrend.high.late(tIdx,:) = ...
        interp1q(es.traj(t_high & es.trialID==output.outcome.lateTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.lateTrials_high(tIdx),start:stop),2), (1:100)');
    plot(es.traj(t_high & es.trialID==output.outcome.lateTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.lateTrials_high(tIdx),start:stop),2),'color',[0 0 1]);
%     diffs.late = [diffs.late diff(Posterior_all.MAP.high(es.trialID(t_high)==output.outcome.lateTrials_high(tIdx) & (Posterior_all.X_high<start)))];
        
end
for tIdx = 1:length(output.outcome.wrongTrials_high)'
    meanPostTrend.high.wrong(tIdx,:) = ...
        interp1q(es.traj(t_high & es.trialID==output.outcome.wrongTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.wrongTrials_high(tIdx),start:stop),2), (1:100)');
    plot(es.traj(t_high & es.trialID==output.outcome.wrongTrials_high(tIdx)), nanmean(P_high(es.trialID(t_high)==output.outcome.wrongTrials_high(tIdx),start:stop),2),'color',[1 0 0]);
%     diffs.wrong = [diffs.wrong diff(Posterior_all.MAP.high(es.trialID(t_high)==output.outcome.wrongTrials_high(tIdx) & (Posterior_all.X_high<start)))];
        %     2*Posterior_all.X_high(es.trialID(t_high)==tIdx),
end
meanPostTrend.all.correct = [meanPostTrend.all.correct' meanPostTrend.high.correct']';
meanPostTrend.all.miss =    [meanPostTrend.all.miss' meanPostTrend.high.miss']';
meanPostTrend.all.late =    [meanPostTrend.all.late' meanPostTrend.high.late']';
meanPostTrend.all.wrong =   [meanPostTrend.all.wrong' meanPostTrend.high.wrong']';
try
errorbar(1:100, nanmean(meanPostTrend.high.correct,1), nansem(meanPostTrend.high.correct),'k')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.high.miss,1), nansem(meanPostTrend.high.miss),'b')
catch; end
try
    errorbar(1:100, nanmean(meanPostTrend.high.late,1), nansem(meanPostTrend.high.late),'c')
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
meanPostTrend.all.late = meanPostTrend.all.late./norm_factor;
es_new.norm_factor = norm_factor;
meanPostTrend.es = es_new;

figure;
try
    errorbar(1:100, nanmean(meanPostTrend.all.correct,1), nansem(meanPostTrend.all.correct),'k');
    hold on;
    errorbar(1:100, nanmean(meanPostTrend.all.miss,1), nansem(meanPostTrend.all.miss),'b')
    errorbar(1:100, nanmean(meanPostTrend.all.late,1), nansem(meanPostTrend.all.late),'c')
    errorbar(1:100, nanmean(meanPostTrend.all.wrong,1), nansem(meanPostTrend.all.wrong),'r')
catch; end;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
xlabel('Position (cm)')
ylabel(['Posterior in ' num2str(2*start) ' to ' num2str(2*stop)]);
line([70-5 70-5],ylim,'color','c');
line([70+5 70+5],ylim,'color','c');

%% plot licks on decoded position vs. actual position
figure;
subplot(121)
patch(2*[start stop stop start],2*[start start stop stop],'c');
hold on;

axis equal; axis square;
set(gca,'XLim',[0 100],'YLim',[0 100])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, ylim, 'linestyle','--', 'color',[.5 .5 .5]);
line(xlim, 2*[start start ],'color','c');
line(xlim, 2*[stop stop ],'color','c');
line(2*[start start],ylim, 'color','c');
line(2*[stop stop],ylim, 'color','c');
correct_times = zeros(size(es.traj));
for temp = [output.outcome.correctTrials_high' output.outcome.correctTrials_norm' output.outcome.correctTrials_low']
    correct_times(es.trialID==temp) = 1;
end
meanPostTrend.licks.correct.traj = es.traj(correct_times & es.lick & es.traj~=0);
meanPostTrend.licks.correct.X_ML = 2*es_new.X_ML(correct_times & es.lick & es.traj~=0);
plot(es.traj(correct_times & es.lick)+randn(1)*0.1, 2*es_new.X_ML(correct_times & es.lick)+randn(1)*0.1,'k.'); 
xlabel('Actual position (cm)');
ylabel('MLE decoded position (cm)');

subplot(122)
patch(2*[start stop stop start],2*[start start stop stop],'c');
hold on;
axis equal; axis square;
set(gca,'XLim',[0 100],'YLim',[0 100])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, ylim, 'linestyle','--', 'color',[.5 .5 .5]);
line(xlim, 2*[start start ],'color','c');
line(xlim, 2*[stop stop ],'color','c');
line(2*[start start],ylim, 'color','c');
line(2*[stop stop],ylim, 'color','c');
wrong_times = zeros(size(es.traj));
for temp = [output.outcome.wrongTrials_high' output.outcome.wrongTrials_norm' output.outcome.wrongTrials_low']
    wrong_times(es.trialID==temp) = 1;
end
meanPostTrend.licks.wrong.traj = es.traj(wrong_times & es.lick & es.traj~=0);
meanPostTrend.licks.wrong.X_ML = 2*es_new.X_ML(wrong_times & es.lick & es.traj~=0);

late_times = zeros(size(es.traj));
for temp = [output.outcome.lateTrials_high' output.outcome.lateTrials_norm' output.outcome.lateTrials_low']
    if ~isempty(temp)
    late_times(es.trialID==temp) = 1;
    end
end
meanPostTrend.licks.late.traj = es.traj(late_times & es.lick & es.traj~=0);
meanPostTrend.licks.late.X_ML = 2*es_new.X_ML(late_times & es.lick & es.traj~=0);
plot(es.traj(wrong_times & es.lick)+randn(1)*0.1, 2*es_new.X_ML(wrong_times & es.lick)+randn(1)*0.1,'r.'); 
plot(es.traj(late_times & es.lick)+randn(1)*0.1, 2*es_new.X_ML(late_times & es.lick)+randn(1)*0.1,'g.'); 
xlabel('Actual position (cm)');
ylabel('MLE decoded position (cm)');

% meanPostTrend.diffs = diffs;
% 
% X = 0:5:25;
% dm = hist(abs(diffs.miss),X);
% dc = hist(abs(diffs.correct),X);
% dw = hist(abs(diffs.wrong),X);
% dm = dm./sum(dm);
% dc = dc./sum(dc);
% dw = dw./sum(dw);
% figure;
% plot(X, dm, 'k', X, dc, 'c',X, dw,'r')

%%
    function output = getsumPost(X)
        Y = (2.^X)*baseline;
        output = log2(nansum(Y,2)./baseline);
    end
end
