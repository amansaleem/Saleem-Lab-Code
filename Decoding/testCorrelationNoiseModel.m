Ac = zeros(20,20,20);
Sc = Ac; Acov = Ac; Scov = Ac;

for iCN = 1:20
    for iX = 1:20
        for iY = 1:20
            [~, ~, ~, Ac(iCN, iX, iY), Sc(iCN, iX, iY),...
                Acov(iCN, iX, iY), Scov(iCN, iX, iY)] = ...
                test_syntheticNoise(iX-1, iY-1, iCN-1, 0);
            drawnow
        end
        display(['Done sigma_x = ' num2str(iX-1)])
    end
    display(['Done sigma_com = ' num2str(iCN-1)])
    figure(1)
    subplot(5,4,iCN)
    imagesc(0:19, 0:19, sq(Ac(iCN,:,:))); 
    axis xy; axis equal; axis square; axis tight; colorbar
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    title(['Corr, common noise = ' num2str(iCN-1)])
    
    figure(2)
    subplot(5,4,iCN)
    imagesc(0:19, 0:19, sq(Sc(iCN,:,:))); 
    axis xy; axis equal; axis square; axis tight; colorbar
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    title(['Corr (S), common noise = ' num2str(iCN-1)])
    
    figure(3)
    subplot(5,4,iCN)
    imagesc(0:19, 0:19, sq(Acov(iCN,:,:))); 
    axis xy; axis equal; axis square; axis tight; colorbar
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    title(['Cov, common noise = ' num2str(iCN-1)])
    
    figure(4)
    subplot(5,4,iCN)
    imagesc(0:19, 0:19, sq(Scov(iCN,:,:))); 
    axis xy; axis equal; axis square; axis tight; colorbar
    set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
    title(['Cov (S), common noise = ' num2str(iCN-1)])
    drawnow;
    
    figure(5)
end
