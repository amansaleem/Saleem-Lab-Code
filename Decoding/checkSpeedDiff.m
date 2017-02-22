function [meanSpds, semSpds, sigs] =checkSpeedDiff(Posterior_all)

for animalIdx = 1:length(Posterior_all)
    t_low = Posterior_all(animalIdx).t_low;
    X = Posterior_all(animalIdx).data.traj;
    
    t_norm = Posterior_all(animalIdx).t_norm;
        
    t_high = Posterior_all(animalIdx).t_high;
    
    meanSpds(animalIdx,1) = nanmean(Posterior_all(animalIdx).data.smthBallSpd(t_low & X<32));
    meanSpds(animalIdx,2) = nanmean(Posterior_all(animalIdx).data.smthBallSpd(t_norm & X<32));
    meanSpds(animalIdx,3) = nanmean(Posterior_all(animalIdx).data.smthBallSpd(t_high & X<32));
    
    semSpds(animalIdx,1) = nanstd(Posterior_all(animalIdx).data.smthBallSpd(t_low & X<32));
    semSpds(animalIdx,2) = nanstd(Posterior_all(animalIdx).data.smthBallSpd(t_norm & X<32));
    semSpds(animalIdx,3) = nanstd(Posterior_all(animalIdx).data.smthBallSpd(t_high & X<32));
    
    [~,sigs(animalIdx)] = ttest2(Posterior_all(animalIdx).data.smthBallSpd(t_low & X<32), ...
        Posterior_all(animalIdx).data.smthBallSpd(t_high & X<32));
end