% figure
% self = 1;
thres = 0.01;
FRv = 8;

for animalIdx = 1:length(Posterior_all)
    P = Posterior_all(animalIdx);
    EV_n = nanmean(P.decoder.model.EV);
    EV_h = nanmean(P.decoder_high.model.EV);
    EV_l = nanmean(P.decoder_low.model.EV);
%     disp('Running it for shortened track');
%     normEV = nanmean(Posterior_all(animalIdx).decoder_norm_decPos.model.EV_orig,1);
%     lowEV = nanmean(Posterior_all(animalIdx).decoder_low_decPos.model.EV_orig,1);
%     highEV = nanmean(Posterior_all(animalIdx).decoder_high_decPos.model.EV_orig,1);    
    
    t = EV_n>thres & EV_h>thres & EV_l>thres;
    
    respModel_n = FRv*P.decoder.model.meanModel(t,:);
    respModel_h = FRv*P.decoder_high.model.meanModel(t,:);
    respModel_l = FRv*P.decoder_low.model.meanModel(t,:);
%     disp('Shortened track');
%     respModel_n = P.decoder_norm_decPos.model.meanModel_orig(t,:);
%     respModel_h = P.decoder_high_decPos.model.meanModel_orig(t,:);
%     respModel_l = P.decoder_low_decPos.model.meanModel_orig(t,:);
    
    [peakRate_n, peakPos_n] = max(respModel_n,[],2);
    
%     t = ones(1,length(EV_h));%;
%     respModel_h = FRv*P.decoder_high.model.meanModel(t,:);
%     [peakRate_h, peakPos_h] = max(respModel_h,[],2);
    t = EV_h>thres & EV_l>thres & EV_n>thres;% & peakPos_h'<33;
    [peakRate_h, peakPos_h] = max(respModel_h,[],2);
    
%     t = ones(1,length(EV_l));%;
%     respModel_l = FRv*P.decoder_low.model.meanModel(t,:);
%     [peakRate_l, peakPos_l] = max(respModel_l,[],2);
    t = EV_l>thres & EV_h>thres & EV_n>thres;% & peakPos_l'<33;%;
    [peakRate_l, peakPos_l] = max(respModel_l,[],2);
    
%     if ~self
%         peakRate_l = peakRate_n;
%         peakRate_h = peakRate_n;
%         
%         peakPos_l = peakPos_n;
%         peakPos_h = peakPos_n;
%     end
    %     respModel_n = na
    pModel_n = [nan*ones(size(respModel_n)) respModel_n nan*ones(size(respModel_n))];
    mid =  size(respModel_n,2)/2;
    for icell = 1:size(respModel_n,1)
        pModel_n(icell,:) = ...
            ...normalise1var...
            (pModel_n(icell,:))...
            ./(peakRate_n(icell))...
            ..../sqrt(peakRate_n(icell))...
            ;
        pModel_n(icell,:) = circshift(pModel_n(icell,:)', mid-peakPos_n(icell));
    end
    
    pModel_h = [nan*ones(size(respModel_h)) respModel_h nan*ones(size(respModel_h))];
    for icell = 1:size(respModel_h,1)
        pModel_h(icell,:) = ...
            ...normalise1var...
            (pModel_h(icell,:))...
            ./(peakRate_n(icell))...
            ..../sqrt(peakRate_h(icell))...
            ;
        pModel_h(icell,:) = circshift(pModel_h(icell,:)', mid-peakPos_h(icell));
    end
    
    pModel_l = [nan*ones(size(respModel_l)) respModel_l nan*ones(size(respModel_l))];
    for icell = 1:size(respModel_l,1)
        pModel_l(icell,:) = ...
            ...normalise1var...
            (pModel_l(icell,:))...
            ./(peakRate_n(icell))...
            ..../sqrt(peakRate_l(icell))...
            ;
        pModel_l(icell,:) = circshift(pModel_l(icell,:)', mid-peakPos_l(icell));
    end
    middle = size(pModel_h,2)/2;
    midRange = (middle-mid):(middle+mid);
    for n = 1:middle
        for icell = 1:size(pModel_h,1)
            pModel_h(icell,n) = nansum([pModel_h(icell,n), pModel_h(icell,(2*middle)-n+1)])/2;
            pModel_h(icell,(2*middle)-n+1) = pModel_h(icell,n);
        end
        for icell = 1:size(pModel_n,1)
            pModel_n(icell,n) = nansum([pModel_n(icell,n), pModel_n(icell,(middle*2)-n+1)])/2;
            pModel_n(icell,(middle*2)-n+1) = pModel_n(icell,n);
        end
        for icell = 1:size(pModel_l,1)
            pModel_l(icell,n) = nansum([pModel_l(icell,n), pModel_l(icell,(2*middle)-n+1)])/2;
            pModel_l(icell,(2*middle)-n+1) = pModel_l(icell,n);
        end
    end
    
%     for icell = 1:size(pModel_l,1)
%         fieldWidth_l(icell) = max(find(normalise1var(pModel_l(1,51:100)')>0.5))-25;
%     end
%     for icell = 1:size(pModel_h,1)
%         fieldWidth_h(icell) = max(find(normalise1var(pModel_h(1,51:100)')>0.5))-25;
%     end
%     for icell = 1:size(pModel_n,1)
%         fieldWidth_n(icell) = max(find(normalise1var(pModel_n(1,51:100)')>0.5))-25;
%     end
    
%     meanWidths(animalIdx,:) = [mean(fieldWidth_l) mean(fieldWidth_n) mean(fieldWidth_h)];
    
    subplot(3,3,animalIdx)
%     imagesc(pModel_l(:,51:100))
    hold off;
    plot(...
        ...normalise1var...
        (nanmean(pModel_l(:,midRange))),'b');
    hold on;
    plot(...
        ...normalise1var...
        (nanmean(pModel_h(:,midRange))),'r');
%     plot(nanmean(pModel_n(:,51:100)),'k');
    axis tight;
    
    popField_l(animalIdx,:) = nanmean(pModel_l(:,midRange));
    popField_n(animalIdx,:) = nanmean(pModel_n(:,midRange));
    popField_h(animalIdx,:) = nanmean(pModel_h(:,midRange));
    
    
%     pause;
end
%%
subplot(3,3,9);
hold off
plot(...
    ...normalise1var...
    (nanmean(popField_l)),'b');
hold on;
plot(...
    ...normalise1var...
    (nanmean(popField_h)),'r');
% plot(-49:2:50,normalise1var(nanmean(popField_n)),'k');
    