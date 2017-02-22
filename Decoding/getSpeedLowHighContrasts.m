for animalIdx = 1:length(Posterior_all)
    P = Posterior_all(animalIdx);
    es = P.data;
    
    for iPos = 1:50
        posBasedSpeed.lc(animalIdx,iPos) = nanmean(es.smthBallSpd(es.contrast>0 & es.contrast<0.6 & round(es.traj/2)==iPos & es.outcome ==2));
        posBasedSpeed.hc(animalIdx,iPos) = nanmean(es.smthBallSpd(es.contrast>0 & es.contrast>0.6 & round(es.traj/2)==iPos & es.outcome==2));
        posBasedSpeed.le(animalIdx,iPos) = nanmean(es.smthBallSpd(es.contrast>0 & es.contrast<0.6 & round(es.traj/2)==iPos & es.outcome ~=2));
        posBasedSpeed.he(animalIdx,iPos) = nanmean(es.smthBallSpd(es.contrast>0 & es.contrast>0.6 & round(es.traj/2)==iPos & es.outcome~=2));
    end
    lowContrastSpeeds{animalIdx}  = es.smthBallSpd(P.t_low);
    highContrastSpeeds{animalIdx} = es.smthBallSpd(P.t_high);
    
    meanLowSpeed(animalIdx)  = nanmean(lowContrastSpeeds{animalIdx});
    meanHighSpeed(animalIdx) = nanmean(highContrastSpeeds{animalIdx});
    
    semLowSpeed(animalIdx)  = nansem(lowContrastSpeeds{animalIdx});
    semHighSpeed(animalIdx) = nansem(highContrastSpeeds{animalIdx});
    
    stdLowSpeed(animalIdx)  = nanstd(lowContrastSpeeds{animalIdx});
    stdHighSpeed(animalIdx) = nanstd(highContrastSpeeds{animalIdx});
end

% errorbarxy(meanLowSpeed, meanHighSpeed, ...
%     ...stdLowSpeed, stdHighSpeed, ...
%     semLowSpeed, semHighSpeed, ...
%     {'ko','k','k'})
figure(1)
plot(meanLowSpeed, meanHighSpeed, 'ko', 'MarkerSize', 10)
axis equal
axis([5 25 5 25])
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
line(xlim, ylim, 'color','k','linestyle','--')
set(gca, 'XTick', [10 15 20 25])
set(gca, 'YTick', [10 15 20 25])
xlabel('Mean speed at low contrast (cm/s)')
ylabel('Mean speed at high contrast (cm/s)')

%%
figure(2);
hold off
plot(nanmean(posBasedSpeed.lc,1),'b')
hold on;
plot(nanmean(posBasedSpeed.hc,1),'r')
plot(nanmean(posBasedSpeed.le,1),'b--')
plot(nanmean(posBasedSpeed.he,1),'r--')
pause

for ianimal = 1:size(posBasedSpeed.lc,1)
    hold off
    plot(posBasedSpeed.lc(ianimal,:),'b')
    hold on;
    plot(posBasedSpeed.hc(ianimal,:),'r')
    plot(posBasedSpeed.le(ianimal,:),'b--')
    plot(posBasedSpeed.he(ianimal,:),'r--')
    pause
end