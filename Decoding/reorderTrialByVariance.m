function [trialIDs] = reorderTrialByVariance(Posterior_all, lims)
% This takes as input the Posterior_all that is the output of the function
% checkDrift_slide
if nargin>1 && length(lims)==2
    minLim = lims(1);
    maxLim = lims(2);
else
    minLim = 10;
    maxLim = 25;
end

Posterior_all = checkDrift_slide(Posterior_all, 5, 1, 1, 0);

for animalIdx = 1:length(Posterior_all)
    clear trial_list_high trial_list_norm trial_list_low
    clear trial_variance_high trial_variance_norm trial_variance_low
    % High contrast trials
    trial_list = (Posterior_all(animalIdx).data.trialID(Posterior_all(animalIdx).t_high));
    trial_list_high = unique(trial_list);
    for itrial = 1:length(trial_list_high)
        trial_variance_high(itrial) = ...
            nanmean(Posterior_all(animalIdx).slideVar.high(trial_list == trial_list_high(itrial) & ...
            Posterior_all(animalIdx).X_high >= minLim & ...
            Posterior_all(animalIdx).X_high <= maxLim));
    end
    % Norm contrast trials
    trial_list = (Posterior_all(animalIdx).data.trialID(Posterior_all(animalIdx).t_norm));
    trial_list_norm = unique(trial_list);
    for itrial = 1:length(trial_list_norm)
        trial_variance_norm(itrial) = ...
            nanmean(Posterior_all(animalIdx).slideVar.norm(trial_list == trial_list_norm(itrial) & ...
            Posterior_all(animalIdx).X_norm >= minLim & ...
            Posterior_all(animalIdx).X_norm <= maxLim));
    end
    % Low contrast trials
    trial_list = (Posterior_all(animalIdx).data.trialID(Posterior_all(animalIdx).t_low));
    trial_list_low = unique(trial_list);
    for itrial = 1:length(trial_list_low)
        trial_variance_low(itrial) = ...
            nanmean(Posterior_all(animalIdx).slideVar.low(trial_list == trial_list_low(itrial) & ...
            Posterior_all(animalIdx).X_low >= minLim & ...
            Posterior_all(animalIdx).X_low <= maxLim));
    end
    trialIDs(animalIdx).trial_list = [trial_list_low' trial_list_norm' trial_list_high'];
    trialIDs(animalIdx).variance_val = [trial_variance_low trial_variance_norm trial_variance_high];
    trialIDs(animalIdx).contrast_val = [0.18*ones(1,length(trial_list_low)) 0.60*ones(1,length(trial_list_norm)) 0.72*ones(1,length(trial_list_high))];
    
end