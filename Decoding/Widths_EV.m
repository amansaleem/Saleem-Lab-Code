thres = 0.01;
clear a_*
for animalIdx = 1:length(Posterior_all)
    figure(animalIdx+500)
    P = Posterior_all(animalIdx);
    EV_n = nanmean(P.decoder.model.EV,1);
    EV_nd = nanmean(P.decoder_norm_decPos.model.EV,1);
    
    EV_l = nanmean(P.decoder_low.model.EV,1);
    EV_ld = nanmean(P.decoder_low_decPos.model.EV,1);
    
    EV_h = nanmean(P.decoder_high.model.EV,1);
    EV_hd = nanmean(P.decoder_high_decPos.model.EV,1);
    
    t = (EV_l>thres | EV_h>thres) & EV_n>thres;
    t_d = (EV_ld>thres | EV_hd>thres);% & EV_nd>thres;
    t_all = t | t_d;
    
    diffEV(animalIdx) = nanmean(EV_h(t) - EV_l(t));
    diffEV_decP(animalIdx) = nanmean(EV_hd(t_d) - EV_ld(t_d));
    diffEVlow(animalIdx) = nanmean(EV_l(t_all) - EV_ld(t_all));
    diffEVhigh(animalIdx) = nanmean(EV_h(t_all) - EV_hd(t_all));
%     t_all = t & t_d;
%     t = t_all;
%     t_d = t_all;

%     a_h(animalIdx,:) = mean(plotPlaceFields(P.decoder_high, t));
%     a_l(animalIdx,:) = mean(plotPlaceFields(P.decoder_low, t));
%     a_n(animalIdx,:) = mean(plotPlaceFields(P.decoder, t));
%     
%     a_hd(animalIdx,:) = mean(plotPlaceFields(P.decoder_high_decPos, t_d));
%     a_ld(animalIdx,:) = mean(plotPlaceFields(P.decoder_low_decPos, t_d));
%     a_nd(animalIdx,:) = mean(plotPlaceFields(P.decoder_norm_decPos, t_d));
    
%     all_h{animalIdx} = (100/250)*(plotPlaceFields(P.decoder_high,...
%         EV_h>thres & EV_l>thres & EV_n>thres,0));%EV_l>thres & EV_h>thres & EV_n>thres));
%     all_l{animalIdx} = (100/250)*(plotPlaceFields(P.decoder_low,...
%         EV_h>thres & EV_l>thres & EV_n>thres,0));%EV_h>thres & EV_l>thres & EV_n>thres));
%     all_n{animalIdx} = (100/250)*(plotPlaceFields(P.decoder,...
%         EV_n>thres));
%     
%     all_hd{animalIdx} = (100/250)*(plotPlaceFields(P.decoder_high_decPos,...
%         EV_hd>thres & EV_ld>thres,0));
%     all_ld{animalIdx} = (100/250)*(plotPlaceFields(P.decoder_low_decPos,...
%         EV_hd>thres & EV_ld>thres,0));%EV_ld>thres & EV_hd>thres));
%     all_nd{animalIdx} = (100/250)*(plotPlaceFields(P.decoder_norm_decPos,...
%         EV_nd>thres));
    
    subplot(2,3,3)
    a_h(animalIdx,:) = (100/250)*mean(plotPlaceFields(P.decoder_high,...
        EV_h>thres & EV_l>thres & EV_n>thres,1));%EV_l>thres & EV_h>thres & EV_n>thres));
    title('High')
    subplot(2,3,2)
    a_l(animalIdx,:) = (100/250)*mean(plotPlaceFields(P.decoder_low,...
        EV_l>thres & EV_h>thres & EV_n>thres,1));%EV_l>thres & EV_ld>thres,1));%
    title('Low')
    subplot(2,3,1)
    a_n(animalIdx,:) = (100/250)*mean(plotPlaceFields(P.decoder,...
        EV_n>thres));
    title('Medium')
    
    subplot(2,3,6)
    a_hd(animalIdx,:) = (100/250)*mean(plotPlaceFields(P.decoder_high_decPos,...
        EV_hd>thres & EV_ld>thres,1));
    title('High Decoded')
    subplot(2,3,5)
    a_ld(animalIdx,:) = (100/250)*mean(plotPlaceFields(P.decoder_low_decPos,...
        EV_l>thres & EV_ld>thres,1));%EV_ld>thres & EV_hd>thres));
    title('Low Decoded')
    subplot(2,3,4)
    a_nd(animalIdx,:) = (100/250)*mean(plotPlaceFields(P.decoder_norm_decPos,...
        EV_nd>thres));
    title('Medium Decoded')
    numCells(animalIdx) = sum(t);
    numCells_d(animalIdx) = sum(t_d);
end

    