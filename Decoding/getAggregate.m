function [corr_time] = getAggregate(Posterior_all, lags, time_bin)

if nargin<2
    lags = 10;
end
if nargin<3
    time_bin = 1;
end
    
for animalIdx = 1:length(Posterior_all)
    
    bump_stats = getActivityBumpStats(Posterior_all(animalIdx));
%     close
    X = nan*ones(size(Posterior_all(animalIdx).t_norm));
    X2 = nan*ones(length(Posterior_all(animalIdx).t_norm), size(Posterior_all(animalIdx).Posterior_norm,2));
    X3 = nan*ones(length(Posterior_all(animalIdx).t_norm), size(bump_stats.bumps.low,1));
    
    Aggregate(animalIdx).X         = X;
    Aggregate(animalIdx).MAP       = X;
    Aggregate(animalIdx).error     = X;
    Aggregate(animalIdx).Posterior = X2;
    Aggregate(animalIdx).nonNormPosterior = X2;
    Aggregate(animalIdx).Bumps     = X3;
    Aggregate(animalIdx).es     = Posterior_all(animalIdx).data;
    
    t_low = Posterior_all(animalIdx).t_low;%   & ~(Posterior_all(animalIdx).data.outcome>3);
    t_norm = Posterior_all(animalIdx).t_norm;% & ~(Posterior_all(animalIdx).data.outcome>3);
    t_high = Posterior_all(animalIdx).t_high;% & ~(Posterior_all(animalIdx).data.outcome>3);
    t_gray = Posterior_all(animalIdx).data.contrast==0 & Posterior_all(animalIdx).data.traj>0;
    
    t_low_noMiss  = Posterior_all(animalIdx).t_low  & ~(Posterior_all(animalIdx).data.outcome>3);
    t_norm_noMiss = Posterior_all(animalIdx).t_norm & ~(Posterior_all(animalIdx).data.outcome>3);
    t_high_noMiss = Posterior_all(animalIdx).t_high & ~(Posterior_all(animalIdx).data.outcome>3);
    
%     Aggregate(animalIdx).t_low_correct  = Posterior_all(animalIdx).t_low   & (Posterior_all(animalIdx).data.outcome==2);
%     Aggregate(animalIdx).t_low_wrong    = Posterior_all(animalIdx).t_low   & (Posterior_all(animalIdx).data.outcome==1 | Posterior_all(animalIdx).data.outcome==3);
%     Aggregate(animalIdx).t_norm_correct = Posterior_all(animalIdx).t_norm  & (Posterior_all(animalIdx).data.outcome==2);
%     Aggregate(animalIdx).t_norm_wrong   = Posterior_all(animalIdx).t_norm  & (Posterior_all(animalIdx).data.outcome==1 | Posterior_all(animalIdx).data.outcome==3);
%     Aggregate(animalIdx).t_high_correct = Posterior_all(animalIdx).t_high  & (Posterior_all(animalIdx).data.outcome==2);
%     Aggregate(animalIdx).t_high_wrong   = Posterior_all(animalIdx).t_high  & (Posterior_all(animalIdx).data.outcome==1 | Posterior_all(animalIdx).data.outcome==3);
    
    Aggregate(animalIdx).t_correct = (t_low | t_norm | t_high) & (Posterior_all(animalIdx).data.outcome==2);
    Aggregate(animalIdx).t_wrong   = (t_low | t_norm | t_high) & (Posterior_all(animalIdx).data.outcome==1 | Posterior_all(animalIdx).data.outcome==3);
    Aggregate(animalIdx).t_correct_l = (t_low ) & (Posterior_all(animalIdx).data.outcome==2);
    Aggregate(animalIdx).t_wrong_l   = (t_low ) & (Posterior_all(animalIdx).data.outcome==1 | Posterior_all(animalIdx).data.outcome==3);
    Aggregate(animalIdx).t_correct_n = (t_norm) & (Posterior_all(animalIdx).data.outcome==2);
    Aggregate(animalIdx).t_wrong_n   = (t_norm) & (Posterior_all(animalIdx).data.outcome==1 | Posterior_all(animalIdx).data.outcome==3);
    Aggregate(animalIdx).t_correct_h = (t_high) & (Posterior_all(animalIdx).data.outcome==2);
    Aggregate(animalIdx).t_wrong_h   = (t_high) & (Posterior_all(animalIdx).data.outcome==1 | Posterior_all(animalIdx).data.outcome==3);
        
    Aggregate(animalIdx).t_low  = t_low;  %Posterior_all(animalIdx).t_low & ~(Posterior_all(animalIdx).data.outcome>3);
    Aggregate(animalIdx).t_norm = t_norm; %Posterior_all(animalIdx).t_norm & ~(Posterior_all(animalIdx).data.outcome>3);
    Aggregate(animalIdx).t_high = t_high; %Posterior_all(animalIdx).t_high & ~(Posterior_all(animalIdx).data.outcome>3);
    Aggregate(animalIdx).t_low_noMiss  = t_low_noMiss;  %Posterior_all(animalIdx).t_low & ~(Posterior_all(animalIdx).data.outcome>3);
    Aggregate(animalIdx).t_norm_noMiss = t_norm_noMiss; %Posterior_all(animalIdx).t_norm & ~(Posterior_all(animalIdx).data.outcome>3);
    Aggregate(animalIdx).t_high_noMiss = t_high_noMiss; %Posterior_all(animalIdx).t_high & ~(Posterior_all(animalIdx).data.outcome>3);
    Aggregate(animalIdx).t_gray = t_gray;
    
    Aggregate(animalIdx).bump_bins = bump_stats.bins;
    
    % Normal condition
    Aggregate(animalIdx).X(t_norm)         = Posterior_all(animalIdx).X_norm;
    Aggregate(animalIdx).MAP(t_norm)       = Posterior_all(animalIdx).MAP.norm;
    Aggregate(animalIdx).error(t_norm)       = Posterior_all(animalIdx).error.norm;
    Aggregate(animalIdx).Posterior(t_norm,:) = Posterior_all(animalIdx).Posterior_norm;
    Aggregate(animalIdx).Bumps(t_norm,:)   = bump_stats.bumps.norm';
    
    % Low condition.
    Aggregate(animalIdx).X(t_low)         = Posterior_all(animalIdx).X_low;
    Aggregate(animalIdx).MAP(t_low)       = Posterior_all(animalIdx).MAP.low;
    Aggregate(animalIdx).error(t_low)       = Posterior_all(animalIdx).error.low;
    Aggregate(animalIdx).Posterior(t_low,:) = Posterior_all(animalIdx).Posterior_low;
    Aggregate(animalIdx).Bumps(t_low,:)   = bump_stats.bumps.low';
    
    % High condition.
    Aggregate(animalIdx).X(t_high)         = Posterior_all(animalIdx).X_high;
    Aggregate(animalIdx).MAP(t_high)       = Posterior_all(animalIdx).MAP.high;
    Aggregate(animalIdx).error(t_high)       = Posterior_all(animalIdx).error.high;
    Aggregate(animalIdx).Posterior(t_high,:) = Posterior_all(animalIdx).Posterior_high;
    Aggregate(animalIdx).Bumps(t_high,:)   = bump_stats.bumps.high';
    
    try
        Aggregate(animalIdx).nonNormPosterior(t_norm,:) = Posterior_all(animalIdx).nonNormPosterior_norm;
        Aggregate(animalIdx).nonNormPosterior(t_low,:) = Posterior_all(animalIdx).nonNormPosterior_low;
        Aggregate(animalIdx).nonNormPosterior(t_high,:) = Posterior_all(animalIdx).nonNormPosterior_high;
    catch; end
end

%% plot posterior and actual position
for n = 1:length(Posterior_all)
    figure(n+20)
    
    subplot(3,3,[1 2])
    t = Aggregate(n).t_low_noMiss;
    plot(Aggregate(n).X(t),'k');
    hold on;plot(Aggregate(n).MAP(t),'r', 'linewidth',1.1);
    licks = Posterior_all(n).data.lick(t);
    X = Aggregate(n).X(t); MAP = Aggregate(n).MAP(t);
    plot(find(licks>0),MAP(licks>0), 'b+')
    subplot(3,3,3)
    plot(Aggregate(n).X(t), Aggregate(n).MAP(t),'k');
    axis equal;    axis([0 50 0 50]);
    line(xlim, ylim, 'color', 'k');
    
    subplot(3,3,[4 5])
    t = Aggregate(n).t_norm_noMiss;
    plot(Aggregate(n).X(t),'k');
    hold on; plot(Aggregate(n).MAP(t),'r', 'linewidth',1.1);
    licks = Posterior_all(n).data.lick(t);
    X = Aggregate(n).X(t); MAP = Aggregate(n).MAP(t);
    plot(find(licks>0),MAP(licks>0), 'b+')
    subplot(3,3,6)
    plot(Aggregate(n).X(t), Aggregate(n).MAP(t),'k');
    axis equal;    axis([0 50 0 50]);
    line(xlim, ylim, 'color', 'k');
    
    subplot(3,3,[7 8])
    t = Aggregate(n).t_high_noMiss;
    plot(Aggregate(n).X(t),'k');
    hold on;plot(Aggregate(n).MAP(t),'r', 'linewidth',1.1);
    licks = Posterior_all(n).data.lick(t);
    X = Aggregate(n).X(t); MAP = Aggregate(n).MAP(t);
    plot(find(licks>0),MAP(licks>0), 'b+')
    subplot(3,3,9)
    plot(Aggregate(n).X(t), Aggregate(n).MAP(t),'k');
    axis equal;    axis([0 50 0 50]);
    line(xlim, ylim, 'color', 'k');
%     pause
end
%% plot posterior, actual & decoded position 

figure(1)

for n = 1:length(Posterior_all)%[1 4:8] %
    
    subplot(311)
    t = Aggregate(n).t_low_noMiss;
    hold off;
    imagesc(1:sum(t), Aggregate(n).bump_bins,Aggregate(n).Bumps(t,:)'); axis xy; RedWhiteBlue;
    hold on;
    plot(Aggregate(n).X(t),'k');
    plot(Aggregate(n).MAP(t),'r');
    set(gca,'CLim',[-5 5])
    title(num2str(n))
    
    subplot(312)
    t = Aggregate(n).t_norm_noMiss;
    hold off;
    imagesc(1:sum(t), Aggregate(n).bump_bins,Aggregate(n).Bumps(t,:)'); axis xy; RedWhiteBlue;
    hold on;
    plot(Aggregate(n).X(t),'k');
    plot(Aggregate(n).MAP(t),'r');
    set(gca,'CLim',[-5 5])
    
    subplot(313)
    t = Aggregate(n).t_high_noMiss;
    hold off;
    imagesc(1:sum(t), Aggregate(n).bump_bins, Aggregate(n).Bumps(t,:)'); axis xy; RedWhiteBlue;
    hold on;
    plot(Aggregate(n).X(t),'k');
    plot(Aggregate(n).MAP(t),'r');
    set(gca,'CLim',[-5 5])
%     pause
end

%% looking at pairwise-correlations
% obj.skippedCells = find(sum(60*responseModel')<1 | obj.skippedCells);
%             responseModel(obj.skippedCells,:) = []; 
%             
%             % Find the peak rate and peakBin of each cell
%             [obj.peakRates, obj.fieldCentres] = max(responseModel,[],2);
%             [~,obj.cellOrder] = sort(obj.fieldCentres);
figure(105);
subplot(111);
clear corr_time corr_time_*
for animalIdx = 1:length(Aggregate)%[1:5 7 8 6] %
    X_range = Aggregate(animalIdx).X>0 & ...
        Aggregate(animalIdx).X<33;
    EVthres = 0.01;
    spkTrain = Aggregate(animalIdx).es.spikeTrain;
    
    response_model = Posterior_all(animalIdx).decoder.model.meanModel;
    cellList = find(...
        nanmean(Posterior_all(animalIdx).decoder.model.EV) >EVthres...
        ...& nanmean(spkTrain,1) <0.5 ...
        );
    response_model = response_model(cellList,:);
    [~, peakPos]  = max(response_model,[],2);
    meanRate = nanmean(Aggregate(animalIdx).es.spikeTrain(:,cellList));
    [~,cellOrder] = sort(peakPos);
%     [~,cellOrder] = sort(meanRate);
    gm_rate = sqrt((meanRate'*ones(1,length(meanRate))).^2 + (meanRate'*ones(1,length(meanRate)))'.^2);
    
    t_low  = Aggregate(animalIdx).t_low_noMiss  & X_range;
    t_norm = Aggregate(animalIdx).t_norm_noMiss & X_range;
    t_high = Aggregate(animalIdx).t_high_noMiss & X_range;
    t_gray = Aggregate(animalIdx).t_gray;

    
    [starts_low , stops_low ] = getStartStop(Aggregate(animalIdx).t_low,lags);
    [starts_norm, stops_norm] = getStartStop(Aggregate(animalIdx).t_norm,lags);
    [starts_high, stops_high] = getStartStop(Aggregate(animalIdx).t_high,lags);
    [starts_gray, stops_gray] = getStartStop(Aggregate(animalIdx).t_gray,lags);

    [starts_correct, stops_correct] = getStartStop(Aggregate(animalIdx).t_correct,lags);
    [starts_wrong  , stops_wrong  ] = getStartStop(Aggregate(animalIdx).t_wrong,lags);
    
    [starts_correct_l, stops_correct_l] = getStartStop(Aggregate(animalIdx).t_correct_l,lags);
    [starts_wrong_l  , stops_wrong_l  ] = getStartStop(Aggregate(animalIdx).t_wrong_l,lags);
    [starts_correct_n, stops_correct_n] = getStartStop(Aggregate(animalIdx).t_correct_n,lags);
    [starts_wrong_n  , stops_wrong_n  ] = getStartStop(Aggregate(animalIdx).t_wrong_n,lags);
    [starts_correct_h, stops_correct_h] = getStartStop(Aggregate(animalIdx).t_correct_h,lags);
    [starts_wrong_h  , stops_wrong_h  ] = getStartStop(Aggregate(animalIdx).t_wrong_h,lags);
    
%     [starts_low_c , stops_low_c ] = getStartStop(Aggregate(animalIdx).t_low_correct,lags);
%     [starts_low_w , stops_low_w ] = getStartStop(Aggregate(animalIdx).t_low_wrong,lags);
%     [starts_norm_c , stops_norm_c ] = getStartStop(Aggregate(animalIdx).t_norm_correct,lags);
%     [starts_norm_w , stops_norm_w ] = getStartStop(Aggregate(animalIdx).t_norm_wrong,lags);
%     [starts_high_c , stops_high_c ] = getStartStop(Aggregate(animalIdx).t_high_correct,lags);
%     [starts_high_w , stops_high_w ] = getStartStop(Aggregate(animalIdx).t_high_wrong,lags);
    
    corr_low  = corr(spkTrain(t_low,cellList));
    corr_norm = corr(spkTrain(t_norm,cellList));
    corr_high = corr(spkTrain(t_high,cellList));
    corr_gray = corr(spkTrain(t_gray,cellList));
    
    ncorr_low  = corr_low(cellOrder, cellOrder)./gm_rate(cellOrder, cellOrder);
    ncorr_norm = corr_norm(cellOrder,cellOrder)./gm_rate(cellOrder, cellOrder);
    ncorr_high = corr_high(cellOrder,cellOrder)./gm_rate(cellOrder, cellOrder);
    
    numCells = length(cellList);
    normMatrix = ones(numCells);
    normMatrix(eye(numCells)==1) = nan;
    
    corr_low  = corr_low.*normMatrix;
    corr_norm = corr_norm.*normMatrix;
    corr_high = corr_high.*normMatrix;
    corr_gray = corr_gray.*normMatrix;
    
    marginal_low(animalIdx ,:) = get45Marginal(corr_low(cellOrder, cellOrder),20);
    marginal_norm(animalIdx,:) = get45Marginal(corr_norm(cellOrder,cellOrder),20);
    marginal_high(animalIdx,:) = get45Marginal(corr_high(cellOrder,cellOrder),20);
    marginal_gray(animalIdx,:) = get45Marginal(corr_gray(cellOrder,cellOrder),20);
    
    figure;
    subplot(236)
    plot(get45Marginal(corr_low(cellOrder,cellOrder),20),'b')
    hold on;
    plot(get45Marginal(corr_high(cellOrder,cellOrder),20),'r')
    plot(get45Marginal(corr_norm(cellOrder,cellOrder),20),'k')
    plot(get45Marginal(corr_gray(cellOrder,cellOrder),20),'color',[.5 .5 .5])
%     subplot(234)
%     imagesc(corr_low(cellOrder,cellOrder)./gm_rate(cellOrder, cellOrder))
%     title('low/gm(rate)')
%     c_l = get(gca,'Clim');
%     
%     subplot(235)
%     imagesc(corr_norm(cellOrder,cellOrder)./gm_rate(cellOrder, cellOrder))
%     title('norm/gm(rate)')
%     c_n = get(gca,'Clim');
%     
%     subplot(236)
%     imagesc(corr_high(cellOrder,cellOrder)./gm_rate(cellOrder, cellOrder))
%     title('high/gm(rate)')
%     c_h = get(gca,'Clim');
%     
%     c_all= [min([c_l c_n c_h]) max([c_l c_n c_h])];
%     for pidx = 1:3
%         subplot(2,3,pidx)
%         set(gca,'CLim', c_all);
%         axis xy; axis equal; axis tight
%         colorbar
%     end
%     figure
    subplot(231)
    imagesc(corr_low(cellOrder,cellOrder))
    title('low')
    
    subplot(232)
    imagesc(corr_norm(cellOrder,cellOrder))
    title('norm')
    
    subplot(233)
    imagesc(corr_high(cellOrder,cellOrder))
    title('high')
    
    subplot(234)
    imagesc(corr_norm(cellOrder,cellOrder) - corr_low(cellOrder,cellOrder))
    title(['norm-low ' num2str(nanmean(corr_norm(:)-corr_low(:)))])
    
    subplot(235)
    imagesc(corr_norm(cellOrder,cellOrder) - corr_high(cellOrder,cellOrder))
    title(['norm-high ' num2str(nanmean(corr_norm(:)-corr_high(:)))])
    
%     subplot(236)
%     imagesc(corr_high(cellOrder,cellOrder) - corr_low(cellOrder,cellOrder))
%     title(['high-low ' num2str(nanmean(corr_high(:)-corr_low(:)))])
    
    for asdasd = 1:6
        subplot(2,3,asdasd)
        RedWhiteBlue; colorbar; axis xy; axis square; axis tight;
        set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
        
        if asdasd>3
            set(gca,'CLim',[-0.15 0.15])
        else
            set(gca,'CLim',[-0.3 0.3])
        end
    end
    
    figure(105)
    subplot(221)
    errorbarxy(nanmean(corr_low(:)),nanmean(corr_high(:)),2*nansem(corr_low(:)),2*nansem(corr_high(:)),{'ko-','k','k'});
    hold on;
    subplot(224)
    errorbarxy(nanmean(abs(corr_norm(:)-corr_low(:))),nanmean(abs(corr_norm(:)-corr_high(:)))...
        ,2*nansem(abs(corr_norm(:)-corr_low(:))),2*nansem(abs(corr_norm(:)-corr_high(:))),{'ko-','k','k'});
    hold on;
    subplot(223)
    errorbarxy(nanmean(abs(ncorr_norm(:)-ncorr_low(:))),nanmean(abs(ncorr_norm(:)-ncorr_high(:)))...
        ,2*nansem(abs(ncorr_norm(:)-ncorr_low(:))),2*nansem(abs(ncorr_norm(:)-ncorr_high(:))),{'ko-','k','k'});
    hold on;
    
    diff_corr(animalIdx) = nanmean(corr_high(:) - corr_low(:));
    diff_corr_p(animalIdx) = signrank(corr_high(:) - corr_low(:));
    
    adiff_corr_low(animalIdx) = nanmean(abs(corr_norm(:) - corr_low(:)));
    adiff_corr_high(animalIdx) = nanmean(abs(corr_norm(:) - corr_high(:)));

    diff_ncorr(animalIdx) = nanmean(ncorr_high(:) - ncorr_low(:));
    diff_ncorr_p(animalIdx) = signrank(ncorr_high(:) - ncorr_low(:));
    
    adiff_ncorr_low(animalIdx) = nanmean(abs(ncorr_norm(:) - ncorr_low(:)));
    adiff_ncorr_high(animalIdx) = nanmean(abs(ncorr_norm(:) - ncorr_high(:)));

    %%
    [corr_time(animalIdx).low , meanSpkRate(animalIdx).low]  = calculateCorrTimeVal(spkTrain, starts_low, stops_low, cellList, lags);
    [corr_time(animalIdx).norm, meanSpkRate(animalIdx).norm] = calculateCorrTimeVal(spkTrain, starts_norm, stops_norm, cellList, lags);
    [corr_time(animalIdx).high, meanSpkRate(animalIdx).high] = calculateCorrTimeVal(spkTrain, starts_high, stops_high, cellList, lags);
    [corr_time(animalIdx).gray, meanSpkRate(animalIdx).gray] = calculateCorrTimeVal(spkTrain, starts_gray, stops_gray, cellList, lags);
    
    if isfield(Posterior_all(animalIdx).data, 'orig')
        % This needs to be plotted if needed
        dataHill = Posterior_all(animalIdx).data.orig.theta.B.hill;
        [thetaPow(animalIdx).low, thetaFreq(animalIdx).low] = getThetaPowFreq(dataHill, starts_low, stops_low);
        [thetaPow(animalIdx).norm, thetaFreq(animalIdx).norm] = getThetaPowFreq(dataHill, starts_norm, stops_norm);
        [thetaPow(animalIdx).high, thetaFreq(animalIdx).high] = getThetaPowFreq(dataHill, starts_high, stops_high);
        
        thetaPow_low(animalIdx) = nanmean(thetaPow(animalIdx).low);
        thetaFreq_low(animalIdx) = nanmean(thetaFreq(animalIdx).low);
        
        thetaPow_norm(animalIdx) = nanmean(thetaPow(animalIdx).norm);
        thetaFreq_norm(animalIdx) = nanmean(thetaFreq(animalIdx).norm);
        
        thetaPow_high(animalIdx) = nanmean(thetaPow(animalIdx).high);
        thetaFreq_high(animalIdx) = nanmean(thetaFreq(animalIdx).high);
        
    end
    corr_time(animalIdx).correct  = calculateCorrTimeVal(spkTrain, starts_correct, stops_correct, cellList, lags);
    corr_time(animalIdx).wrong    = calculateCorrTimeVal(spkTrain, starts_wrong,   stops_wrong  , cellList, lags);

    corr_time(animalIdx).correct_l  = calculateCorrTimeVal(spkTrain, starts_correct_l, stops_correct_l, cellList, lags);
    corr_time(animalIdx).wrong_l    = calculateCorrTimeVal(spkTrain, starts_wrong_l,   stops_wrong_l  , cellList, lags);
    corr_time(animalIdx).correct_n  = calculateCorrTimeVal(spkTrain, starts_correct_n, stops_correct_n, cellList, lags);
    corr_time(animalIdx).wrong_n    = calculateCorrTimeVal(spkTrain, starts_wrong_n,   stops_wrong_n  , cellList, lags);
    corr_time(animalIdx).correct_h  = calculateCorrTimeVal(spkTrain, starts_correct_h, stops_correct_h, cellList, lags);
    corr_time(animalIdx).wrong_h    = calculateCorrTimeVal(spkTrain, starts_wrong_h,   stops_wrong_h  , cellList, lags);

    %
    corr_time_diff_high_low(animalIdx,:) = nanmedian((corr_time(animalIdx).high - corr_time(animalIdx).low)./abs(corr_time(animalIdx).norm),1);
    
    corr_time_diff_wrong(animalIdx,:)  = nanmedian(corr_time(animalIdx).correct - corr_time(animalIdx).wrong,1);

    corr_time_diff_high(animalIdx,:) = nanmedian((corr_time(animalIdx).norm - corr_time(animalIdx).high)./abs(corr_time(animalIdx).norm),1);
    corr_time_diff_low(animalIdx,:)  = nanmedian((corr_time(animalIdx).norm - corr_time(animalIdx).low)./abs(corr_time(animalIdx).norm),1);
    corr_time_diff_gray(animalIdx,:)  = nanmedian((corr_time(animalIdx).norm - corr_time(animalIdx).gray)./abs(corr_time(animalIdx).norm),1);

    corr_time_high(animalIdx,:) = nanmean(corr_time(animalIdx).high,1);
    corr_time_norm(animalIdx,:) = nanmean(corr_time(animalIdx).norm,1);
    corr_time_low(animalIdx,:)  = nanmean(corr_time(animalIdx).low,1);
    corr_time_gray(animalIdx,:)  = nanmean(corr_time(animalIdx).gray,1);
    
    corr_norm_low(animalIdx) = corr(corr_time_low(:), corr_time_norm(:));
    corr_norm_high(animalIdx) = corr(corr_time_high(:), corr_time_norm(:));
    corr_norm_gray(animalIdx) = corr(corr_time_high(:), corr_time_gray(:));
    
end
corr_all.high = [];
corr_all.norm = [];
corr_all.low = [];

for animalIdx = 1:8
    corr_all.high = [corr_all.high corr_time(animalIdx).high(:,2:end)'];
    corr_all.low  = [corr_all.low  corr_time(animalIdx).low(:,2:end)'];
    corr_all.norm = [corr_all.norm corr_time(animalIdx).norm(:,2:end)'];
end
figure(99);
subplot(121);
errorarea_as(2:size(corr_time(animalIdx).high,2), nanmean(corr_all.low'),nansem(corr_all.low'),'b')
errorarea_as(2:size(corr_time(animalIdx).high,2), nanmean(corr_all.high'),nansem(corr_all.high'),'r')

subplot(122);
errorarea_as(2:size(corr_time(animalIdx).high,2), -nanmean(corr_all.low'-corr_all.high'),nansem(corr_all.low'-corr_all.high'),'k')

figure;
subplot(221)
axis square; axis equal;
axis([0 0.15 0 0.15])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, ylim, 'linestyle','--', 'color','k');
xlabel('Mean pairwise corr (low contrast)')
ylabel('Mean pairwise corr (high contrast)')

subplot(223)
axis square; axis equal;
axis([0 0.4 0 0.4])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, ylim, 'linestyle','--', 'color','k');
xlabel('Mean abs ncorr diff (norm - low)')
ylabel('Mean abs ncorr diff (norm - high)')

subplot(224)
axis square; axis equal;
axis([0.01 0.2 0.01 0.2])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, ylim, 'linestyle','--', 'color','k');
xlabel('Mean abs corr diff (norm - low)')
ylabel('Mean abs corr diff (norm - high)')

figure(105)
subplot(222)
plot(corr_norm_low, corr_norm_high,'ko')
axis square; axis equal;
axis([0.975 1 0.975 1]);
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, ylim, 'linestyle','--', 'color','k');
xlabel('Corr (corr_{low}, corr_{norm}')
ylabel('Corr (corr_{high}, corr_{norm}')

% subplot(223)
% plot(adiff_corr_low, adiff_corr_high,'ko')
% axis square; axis equal;
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% line(xlim, ylim, 'linestyle','--', 'color','k');
% xlabel('Mean(abs(corr_{norm}-corr_{low}))'); 
% ylabel('Mean(abs(corr_{norm}-corr_{high}))'); 


figure
subplot(224)
for animalIdx = 1:length(corr_time)
toplot(animalIdx,:) = ...
    nanmean(corr_time(animalIdx).correct_h - ...
    corr_time(animalIdx).wrong_h,1);
end 
errorarea_as((1:size(corr_time_diff_high_low,2)-2)*time_bin,...
    nanmean(toplot(:,3:end)), nansem(toplot(:,3:end)),'r');
% title('High contrast','fontsize',14)
% errorarea_as((1:size(corr_time_high,2)-2)*time_bin, nanmean(corr_time_diff_wrong(:,3:end)), nansem(corr_time_diff_wrong(:,3:end)),'m');
axis square;
axis tight
set(gca, 'XLim',[0 max(lags*time_bin)], 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, [0 0], 'linestyle','--', 'color','k');
xlabel('Delay (bins)')
ylabel('Autocorrelation (Correct - wrong) High')

subplot(221)
errorarea_as((1:size(corr_time_high,2)-2)*time_bin,...
    nanmean(corr_time_high(:,3:end)), nansem(corr_time_high(:,3:end)),'r');
hold on;
errorarea_as((1:size(corr_time_low,2))*time_bin,...
    nanmean(corr_time_low(:,3:end)), nansem(corr_time_low(:,3:end)),'b');
errorarea_as((1:size(corr_time_norm,2))*time_bin,...
    nanmean(corr_time_norm(:,3:end)), nansem(corr_time_norm(:,3:end)),'k');
% errorarea_as((1:size(corr_time_gray,2))*time_bin, nanmean(corr_time_gray(:,3:end)), nansem(corr_time_gray(:,3:end)),'c');
axis square;
axis tight
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% line(xlim, [0 0], 'linestyle','--', 'color','k');
xlabel('Delay (bins)')
ylabel('Autocorrelation')

subplot(223)
errorarea_as((1:size(corr_time_diff_high_low,2)-2)*time_bin,...
    mean(corr_time_diff_high_low(:,3:end)), sem(corr_time_diff_high_low(:,3:end)),'r');
axis square;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, [0 0], 'linestyle','--', 'color','k');
xlabel('Delay (bins)')
ylabel('(high - low) Autocorrelation')

subplot(222)
errorarea_as((1:size(corr_time_diff_high,2)-2)*time_bin, nanmean(corr_time_diff_high(:,3:end)), nansem(corr_time_diff_high(:,3:end)),'r');
hold on;
errorarea_as((1:size(corr_time_diff_low,2)-2)*time_bin, nanmean(corr_time_diff_low(:,3:end)), nansem(corr_time_diff_low(:,3:end)),'b');
% errorarea_as((1:size(corr_time_diff_low,2)-2)*time_bin, nanmean(corr_time_diff_gray(:,3:end)), nansem(corr_time_diff_gray(:,3:end)),'c');
axis square;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
line(xlim, [0 0], 'linestyle','--', 'color','k');
xlabel('Delay (bins)')
ylabel('(Norm - x) Autocorrelation')

%%
figure
subplot(231) % Low: w - c
for animalIdx = 1:length(corr_time)
toplot(animalIdx,:) = ...
    nanmean(corr_time(animalIdx).correct_l - ...
    corr_time(animalIdx).wrong_l,1);
end 
errorarea_as((1:size(corr_time_diff_high_low,2)-2)*time_bin,...
    nanmean(toplot(:,3:end)), nansem(toplot(:,3:end)),'b');
title('Low contrast','fontsize',14)
subplot(233)% High: w - c
for animalIdx = 1:length(corr_time)
toplot(animalIdx,:) = ...
    nanmean(corr_time(animalIdx).correct_h - ...
    corr_time(animalIdx).wrong_h,1);
end 
errorarea_as((1:size(corr_time_diff_high_low,2)-2)*time_bin,...
    nanmean(toplot(:,3:end)), nansem(toplot(:,3:end)),'r');
title('High contrast','fontsize',14)
subplot(232)% Normal: w - c
for animalIdx = 1:length(corr_time)
toplot(animalIdx,:) = ...
    nanmean(corr_time(animalIdx).correct_n - ...
    corr_time(animalIdx).wrong_n,1);
end 
errorarea_as((1:size(corr_time_diff_high_low,2)-2)*time_bin,...
    nanmean(toplot(:,3:end)), nansem(toplot(:,3:end)),'k');
title('Normal contrast','fontsize',14)

subplot(234)% All: w - c
for animalIdx = 1:length(corr_time)
toplot(animalIdx,:) = ...
    nanmean(corr_time(animalIdx).correct - ...
    corr_time(animalIdx).wrong,1);
end 
errorarea_as((1:size(corr_time_diff_high_low,2)-2)*time_bin,...
    nanmean(toplot(:,3:end)), nansem(toplot(:,3:end)),'m');
title('All contrasts','fontsize',14)

subplot(235)% lc - hc
for animalIdx = 1:length(corr_time)
toplot(animalIdx,:) = ...
    nanmean(corr_time(animalIdx).correct_h - ...
    corr_time(animalIdx).correct_l,1);
end 
errorarea_as((1:size(corr_time_diff_high_low,2)-2)*time_bin,...
    nanmean(toplot(:,3:end)), nansem(toplot(:,3:end)),'k');
ylabel('High correct - low correct','fontsize',14)
title('High correct - Low correct','fontsize',14)
subplot(236)% lc - hc
for animalIdx = 1:length(corr_time)
toplot(animalIdx,:) = ...
    nanmean(corr_time(animalIdx).high - ...
    corr_time(animalIdx).low,1);
end 
errorarea_as((1:size(corr_time_diff_high_low,2)-2)*time_bin,...
    nanmean(toplot(:,3:end)), nansem(toplot(:,3:end)),'k');
ylabel('Low - high','fontsize',14)
title('Low - high','fontsize',14)
for idx = 1:6
    subplot(2,3,idx)
    axis square;
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    line(xlim, [0 0], 'linestyle','--', 'color','k');
    xlabel('Delay (bins)')
    axis([0 (lags+1)*time_bin -0.015 0.025])
    if idx<5;ylabel('(Correct - Wrong) Autocorrelation'); end
end
%%
% figure
% errorarea_as(-9:10, mean(marginal_low),sem(marginal_low),  'b');
% errorarea_as(-9:10, mean(marginal_norm),sem(marginal_norm),'k');
% errorarea_as(-9:10, mean(marginal_high),sem(marginal_high),'r');

%% to plot a segment of the data
% n = 8;
% figure(505)
% evtId = 25;
% subplot(413)
% imagesc(log(Aggregate(n).nonNormPosterior(starts_norm(evtId):stops_norm(evtId),:))');
% colorbar
% subplot(414)
% imagesc((Aggregate(n).Posterior(starts_norm(evtId):stops_norm(evtId),:))');
% colorbar
% subplot(4,1,[1 2])
% imagesc(Aggregate(n).es.spikeTrain(starts_norm(evtId):stops_norm(evtId),:)');
% RedWhiteBlue;
% hold on;
% plot(nansum(Aggregate(n).es.spikeTrain(starts_norm(evtId):stops_norm(evtId),cellList)',1)/2,'color',[.5 .5 .5], 'linewidth',2); axis xy
% plot(nansum(Aggregate(n).es.spikeTrain(starts_norm(evtId):stops_norm(evtId),:)',1)/2,'k', 'linewidth',2); axis xy
% plot((Aggregate(n).error(starts_norm(evtId):stops_norm(evtId))'),'r', 'linewidth',2); axis xy
% colorbar

    function [corr_val, meanSpkRate] = calculateCorrTimeVal(spkTrain, starts, stops, cellList, lags)
        index = 1;
        if length(starts)>0
            for iCell = cellList
                for iSeg = 1:length(starts)
                    if sum(spkTrain(starts(iSeg):stops(iSeg),iCell))>2
                        corr_val(index,:,iSeg)...
                            = autocorr(spkTrain(starts(iSeg):stops(iSeg),iCell),lags);
                        meanSpkRate(index,iSeg)...
                            = nanmean(spkTrain(starts(iSeg):stops(iSeg),iCell));
                    else
                        corr_val(index,:,iSeg) = nan*ones(1,lags+1);
                        meanSpkRate(index,iSeg)...
                            = nanmean(spkTrain(starts(iSeg):stops(iSeg),iCell));
                    end
                end
                index = index + 1;
            end
            corr_val = squeeze(nanmean(corr_val,3));
            meanSpkRate = nanmean(meanSpkRate,2);
        else
            corr_val = nan*ones(length(cellList),lags+1);
            meanSpkRate = zeros(length(cellList));
        end
    end

    function [thetaPow, thetaFreq] = getThetaPowFreq(dataHill, starts, stops);
        for iSeg = 1:length(starts)
            thetaPow(iSeg)= mean(abs(dataHill(starts:stops)));
            thetaFreq(iSeg)= mean(diff(dataHill(starts:stops)));
        end
    end
end
