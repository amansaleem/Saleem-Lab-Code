function [meanPow, semPow, sigs, meanAng] =checkThetaPowDiff(Posterior_all)

for animalIdx = 1:length(Posterior_all)
    es = Posterior_all(animalIdx).data.orig;
    
    base = es.smthBallSpd>5 & es.traj<65 & es.gain==1 & es.roomLength==1 & es.traj>0 & es.contrast>0;
    t_low  = base & es.contrast<0.6;
    t_high = base & es.contrast>0.6;
    t_norm = base & es.contrast==0.6;
    
    thetaPow = abs(es.theta.B.hill);
    thetaAng = real(diff(es.theta.B.hill));
    
    meanPow(animalIdx,1) = nanmean(thetaPow(t_low));
    meanPow(animalIdx,2) = nanmean(thetaPow(t_norm));
    meanPow(animalIdx,3) = nanmean(thetaPow(t_high));
    
    meanAng(animalIdx,1) = nanmean(thetaAng(t_low));
    meanAng(animalIdx,2) = nanmean(thetaAng(t_norm));
    meanAng(animalIdx,3) = nanmean(thetaAng(t_high));
    
    semPow(animalIdx,1) = nansem(thetaPow(t_low));
    semPow(animalIdx,2) = nansem(thetaPow(t_norm));
    semPow(animalIdx,3) = nansem(thetaPow(t_high));
    
    [~,sigs(animalIdx)] = ttest2(thetaPow(t_low), ...
        thetaPow(t_high));
end
figure
subplot(121)
errorbarxy(meanPow(:,1), meanPow(:,3), semPow(:,1), semPow(:,3), {'ko', 'k', 'k'})
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis equal
axis([300 1100 300 1100])
line(xlim, ylim, 'linestyle','--', 'color','k');
xlabel('\theta Power, low contrast')
ylabel('\theta Power, High contrast')

subplot(122)
plot(meanAng(:,1), meanAng(:,3), 'ko')
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
axis equal
axis([-pi pi -pi pi])
line(xlim, ylim, 'linestyle','--', 'color','k');
xlabel('\theta Angle, low contrast')
ylabel('\theta Angle, High contrast')