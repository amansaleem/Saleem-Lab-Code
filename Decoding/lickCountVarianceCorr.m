function [corr_val, rho_val, varWithLick] = lickCountVarianceCorr(Posterior_all, lickLims)
if nargin<2
    lickLims = [2 65];
end

%% defining the animal indices
idx = 0;
idx = idx + 1;
expt_list(idx).animal     = 'M130920_BALL';
expt_list(idx).iseries    = 1025;
expt_list(idx).expt_list  = 102:103;
idx = idx + 1;
expt_list(idx).animal     = 'M130918_BALL';
expt_list(idx).iseries    = 1030;
expt_list(idx).expt_list  = 103:105;
idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 530;
expt_list(idx).expt_list  = 104:106;
idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 531;
expt_list(idx).expt_list  = 103:106;
idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 601;
expt_list(idx).expt_list  = 103:106;
idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 602;
expt_list(idx).expt_list  = 102:106;
idx = idx + 1;
expt_list(idx).animal     = 'M140502_BALL';
expt_list(idx).iseries    = 603;
expt_list(idx).expt_list  = 107:110;
idx = idx + 1;
expt_list(idx).animal     = 'M140502_BALL';
expt_list(idx).iseries    = 604;
expt_list(idx).expt_list  = 107:110;

%% Loading the right files

type = 0;

for animalIdx = 1:length(Posterior_all)
    es = VRLoadMultipleExpts(expt_list(animalIdx).animal, ...
        expt_list(animalIdx).iseries, expt_list(animalIdx).expt_list);
    
    es = redefineOutcome(es);
    fullTrialList = unique(es.trialID);
    newTrialOrder = nan*ones(size(es.trialID));
    [trialIDs] = reorderTrialByVariance(Posterior_all(animalIdx));
    trialIDs.trial_list(isnan(trialIDs.variance_val)) = [];
    trialIDs.contrast_val(isnan(trialIDs.variance_val)) = [];
    trialIDs.variance_val(isnan(trialIDs.variance_val)) = [];
    for itrial = 1:length(trialIDs.trial_list)
        if nanmedian(es.outcome(es.trialID==trialIDs.trial_list(itrial)))<4
            trialIDs.numEarlyLicks(itrial) = nansum(es.lick(es.trialID==trialIDs.trial_list(itrial) ...
                & es.traj>lickLims(1) & es.traj<lickLims(2)));
        else
            trialIDs.numEarlyLicks(itrial) = nan;
        end
    end
    lickCounts = unique(trialIDs.numEarlyLicks(~isnan(trialIDs.numEarlyLicks)));
    subplot(3,3,animalIdx)
    plot(trialIDs.numEarlyLicks(~isnan(trialIDs.numEarlyLicks))',...
        trialIDs.variance_val(~isnan(trialIDs.numEarlyLicks))', 'ko');
    for iLick = 1:length(lickCounts)
        k = trialIDs.numEarlyLicks==lickCounts(iLick);
        varWithLick(animalIdx).numLicks(iLick) = lickCounts(iLick);
        varWithLick(animalIdx).variances{iLick} = trialIDs.variance_val(k);
    end
    subplot(3,3,9)
    plot(trialIDs.variance_val(~isnan(trialIDs.numEarlyLicks))', trialIDs.numEarlyLicks(~isnan(trialIDs.numEarlyLicks))', 'r.'); hold on;
    drawnow;
    [corr_val(animalIdx), rho_val(animalIdx)] = corr(trialIDs.variance_val(~isnan(trialIDs.numEarlyLicks) & (trialIDs.numEarlyLicks)>0)', trialIDs.numEarlyLicks(~isnan(trialIDs.numEarlyLicks) & (trialIDs.numEarlyLicks)>0)');
end
