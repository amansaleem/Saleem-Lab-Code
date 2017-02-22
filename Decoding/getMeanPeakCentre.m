clear hPost* fracActv* fActv peakRate meanRate

if ~exist('start_lim')
    start_lim = 101;
end

if ~exist('samp_rate')
    samp_rate = 60;
end

allFields.n = [];
allFields.l = [];
allFields.h = [];

meanRate.g = [];

peakRate.n = [];
meanRate.n = [];
peakCent.n = [];
fracActv.n = [];

peakRate.h = [];
meanRate.h = [];
peakCent.h = [];
fracActv.h = [];

peakRate.l = [];
meanRate.l = [];
peakCent.l = [];
fracActv.l = [];

peakRate.nd = [];
peakRate.hd = [];
peakRate.ld = [];

peakRate.ld_low = [];
peakRate.l_low = [];

lowRate.n = [];
lowRate.h = [];
lowRate.l = [];
lowRate.nd = [];
lowRate.hd = [];
lowRate.ld = [];

thres = 0.01;

for animalIdx = 1:length(Posterior_all)
    %% Load the data
    P = Posterior_all(animalIdx);
    placeFields_n = P.decoder.model.meanModel;
    placeFields_h = P.decoder_high.model.meanModel;
    placeFields_l = P.decoder_low.model.meanModel;
    
    try
        placeFields_nd = P.decoder_norm_decPos.model.meanModel;
        placeFields_hd = P.decoder_high_decPos.model.meanModel;
        placeFields_ld = P.decoder_low_decPos.model.meanModel;
        decodedP = 1;
        EV_nd = nanmean(P.decoder_norm_decPos.model.EV,1);
        EV_hd = nanmean(P.decoder_high_decPos.model.EV,1);
        EV_ld = nanmean(P.decoder_low_decPos.model.EV,1);
    catch
        decodedP = 0;
    end
    EV_n = nanmean(P.decoder.model.EV,1);
    EV_h = nanmean(P.decoder_high.model.EV,1);
    EV_l = nanmean(P.decoder_low.model.EV,1);
    
    t_n = P.t_norm & P.data.traj<start_lim*2;
    t_h = P.t_high & P.data.traj<start_lim*2;
    t_l = P.t_low & P.data.traj<start_lim*2;
    t_g = P.data.contrast==0;
    
    spikeTrain = P.data.spikeTrain;
    
    subsetCells = EV_h>=thres & EV_l>=thres ...
        & EV_n>=thres...
        ...& EV_ld>=thres & EV_hd>=thres...
        ;
    
    if decodedP
        subsetCells_low = EV_ld>=thres & EV_l>=thres;
        subsetCells_dec = EV_ld>=thres ...
            & EV_hd>=thres;
        %         subsetCells = EV_ld>=thres & EV_l>=thres;
        %         subsetCells_dec = subsetCells;
    end
    %% get peak rates and centers
    [peakRate_n, peakPos_n] = max(placeFields_n(subsetCells,:),[],2);
    [peakRate_h, peakPos_h] = max(placeFields_h(subsetCells,:),[],2);
    [peakRate_l, peakPos_l] = max(placeFields_l(subsetCells,:),[],2);
    
    [lowRate_n] = min(placeFields_n(subsetCells,:),[],2);
    [lowRate_h] = min(placeFields_h(subsetCells,:),[],2);
    [lowRate_l] = min(placeFields_l(subsetCells,:),[],2);
    
        peakRate_l = samp_rate*peakRate_l;
        peakRate_h = samp_rate*peakRate_h;
        peakRate_n = samp_rate*peakRate_n;
    if decodedP
        [peakRate_nd, peakPos_nd] = max(placeFields_nd(subsetCells_dec,:),[],2);
        [peakRate_hd, peakPos_hd] = max(placeFields_hd(subsetCells_dec,:),[],2);
        [peakRate_ld, peakPos_ld] = max(placeFields_ld(subsetCells_dec,:),[],2);
        
        [lowRate_nd] = min(placeFields_nd(subsetCells_dec,:),[],2);
        [lowRate_hd] = min(placeFields_hd(subsetCells_dec,:),[],2);
        [lowRate_ld] = min(placeFields_ld(subsetCells_dec,:),[],2);
        
        [peakRate_l_low,  peakPos_l_low ] = max(placeFields_l(subsetCells_low,:),[],2);
        [peakRate_ld_low, peakPos_ld_low] = max(placeFields_ld(subsetCells_low,:),[],2);
        
        peakRate_ld = samp_rate*peakRate_ld;
        peakRate_hd = samp_rate*peakRate_hd;
        peakRate_nd = samp_rate*peakRate_nd;
        
        peakRate_l_low = samp_rate*peakRate_l_low;
        peakRate_ld_low = samp_rate*peakRate_ld_low;
        
    end
    
    %% get Posterior max - min
    if length(Posterior_all)==6
        hPost_n{animalIdx} = max(2.^(P.Posterior_norm(P.X_norm<start_lim,:))/50,[],2)...
            -min(2.^(P.Posterior_norm(P.X_norm<start_lim,:))/50,[],2);
        hPost_h{animalIdx} = max(2.^(P.Posterior_high(P.X_high<start_lim,:))/50,[],2)...
            -min(2.^(P.Posterior_high(P.X_high<start_lim,:))/50,[],2);
        hPost_l{animalIdx} = max(2.^(P.Posterior_low(P.X_low < start_lim,:))/50 ,[],2)...
            -min(2.^(P.Posterior_low(P.X_low < start_lim,:))/50,[],2);
    else
        hPost_n{animalIdx} = max((P.Posterior_norm(P.X_norm<start_lim,:)),[],2)...
            -min((P.Posterior_norm(P.X_norm<start_lim,:)),[],2);
        hPost_h{animalIdx} = max((P.Posterior_high(P.X_high<start_lim,:)),[],2)...
            -min((P.Posterior_high(P.X_high<start_lim,:)),[],2);
        hPost_l{animalIdx} = max((P.Posterior_low(P.X_low < start_lim,:)) ,[],2)...
            -min((P.Posterior_low(P.X_low < start_lim,:)),[],2);
    end
    %% get mean rates
    meanRate_n = samp_rate*nanmean(spikeTrain(t_n,subsetCells),1);
    meanRate_h = samp_rate*nanmean(spikeTrain(t_h,subsetCells),1);
    meanRate_l = samp_rate*nanmean(spikeTrain(t_l,subsetCells),1);
    meanRate_g = samp_rate*nanmean(spikeTrain(t_g,subsetCells),1);
    
    %% Get fraction of cells active
    fracActv_n{animalIdx} = sum(spikeTrain(t_n,subsetCells)>0,2)./sum(subsetCells);
    fracActv_h{animalIdx} = sum(spikeTrain(t_h,subsetCells)>0,2)./sum(subsetCells);
    fracActv_l{animalIdx} = sum(spikeTrain(t_l,subsetCells)>0,2)./sum(subsetCells);
    
    %% collecting them
    peakRate.n = [peakRate.n peakRate_n'];
    meanRate.n = [meanRate.n meanRate_n];
    peakCent.n = [peakCent.n peakPos_n'];
    fracActv.n = [fracActv.n fracActv_n{animalIdx}'];
    
    peakRate.h = [peakRate.h peakRate_h'];
    meanRate.h = [meanRate.h meanRate_h];
    peakCent.h = [peakCent.h peakPos_h'];
    fracActv.h = [fracActv.h fracActv_h{animalIdx}'];
    
    peakRate.l = [peakRate.l peakRate_l'];
    meanRate.l = [meanRate.l meanRate_l];
    peakCent.l = [peakCent.l peakPos_l'];
    fracActv.l = [fracActv.l fracActv_l{animalIdx}'];
    
    lowRate.n = [lowRate.n lowRate_n'];
    lowRate.h = [lowRate.h lowRate_h'];
    lowRate.l = [lowRate.l lowRate_l'];
    
    if decodedP
        peakRate.nd = [peakRate.nd peakRate_nd'];
        peakRate.hd = [peakRate.hd peakRate_hd'];
        peakRate.ld = [peakRate.ld peakRate_ld'];
        
        peakRate.l_low  = [peakRate.l_low peakRate_l_low'];
        peakRate.ld_low = [peakRate.ld_low peakRate_ld_low'];
        
        lowRate.nd = [lowRate.nd lowRate_nd'];
        lowRate.hd = [lowRate.hd lowRate_hd'];
        lowRate.ld = [lowRate.ld lowRate_ld'];
    end
    meanRate.g = [meanRate.g meanRate_g];
    
    mRate.h(animalIdx) = mean(meanRate_h./meanRate_n);
    fActv.h(animalIdx) = mean(fracActv_h{animalIdx});
    hPost.h(animalIdx) = mean(hPost_h{animalIdx});
    
    mRate.l(animalIdx) = mean(meanRate_l./meanRate_n);
    fActv.l(animalIdx) = mean(fracActv_l{animalIdx});
    hPost.l(animalIdx) = mean(hPost_l{animalIdx});
    
    mRate.n(animalIdx) = mean(meanRate_n./meanRate_n);
    fActv.n(animalIdx) = mean(fracActv_n{animalIdx});
    hPost.n(animalIdx) = mean(hPost_n{animalIdx});
    
    lRate.h(animalIdx) = mean(lowRate_h);
    lRate.l(animalIdx) = mean(lowRate_l);
    lRate.n(animalIdx) = mean(lowRate_n);
    
    modRate.h(animalIdx) = mean((peakRate_h-lowRate_h)./(lowRate_h+peakRate_h));
    modRate.l(animalIdx) = mean((peakRate_l-lowRate_l)./(lowRate_l+peakRate_l));
    modRate.n(animalIdx) = mean((peakRate_n-lowRate_n)./(lowRate_n+peakRate_n));
    
    pRate.h(animalIdx) = mean(peakRate_h);% mean(2*peakRate_h./(peakRate_h+peakRate_l));
    pRate.l(animalIdx) = mean(peakRate_l);%mean(2*peakRate_l./(peakRate_h+peakRate_l));
    pRate.n(animalIdx) = mean(peakRate_n);%mean(peakRate_n);
    
    pRate.h_s(animalIdx) = nansem(peakRate_h);%nansem(2*peakRate_h./(peakRate_h+peakRate_l));
    pRate.l_s(animalIdx) = nansem(peakRate_l);%nansem(2*peakRate_l./(peakRate_h+peakRate_l));
    pRate.n_s(animalIdx) = nansem(peakRate_n);%nansem(peakRate_n);
    
    if decodedP
        pRate.hd(animalIdx) = mean(peakRate_hd);%mean(2*peakRate_hd./(peakRate_hd+peakRate_ld));
        pRate.ld(animalIdx) = mean(peakRate_ld);%mean(2*peakRate_ld./(peakRate_hd+peakRate_ld));
        pRate.nd(animalIdx) = mean(peakRate_nd);%mean(peakRate_nd);
    
        pRate.hd_s(animalIdx) = nansem(peakRate_hd);%nansem(2*peakRate_hd./(peakRate_hd+peakRate_ld));
        pRate.ld_s(animalIdx) = nansem(peakRate_ld);%nansem(2*peakRate_ld./(peakRate_hd+peakRate_ld));
        pRate.nd_s(animalIdx) = nansem(peakRate_nd);
    
        pRate.l_low(animalIdx) = mean(peakRate_l_low);%mean(2*peakRate_l_low./(peakRate_l_low+peakRate_ld_low));
        pRate.ld_low(animalIdx) = mean(peakRate_ld_low);%mean(2*peakRate_ld_low./(peakRate_l_low+peakRate_ld_low));
        
        pRate.l_low_s(animalIdx) = nansem(peakRate_l_low);%nansem(2*peakRate_l_low./(peakRate_l_low+peakRate_ld_low));
        pRate.ld_low_s(animalIdx) = nansem(peakRate_ld_low);%nansem(2*peakRate_ld_low./(peakRate_l_low+peakRate_ld_low));
        
        lRate.hd(animalIdx) = mean(lowRate_hd);
        lRate.ld(animalIdx) = mean(lowRate_ld);
        lRate.nd(animalIdx) = mean(lowRate_nd);
        
        modRate.hd(animalIdx) = mean((peakRate_hd-lowRate_hd)./(lowRate_hd+peakRate_hd));
        modRate.ld(animalIdx) = mean((peakRate_ld-lowRate_ld)./(lowRate_ld+peakRate_ld));
        modRate.nd(animalIdx) = mean((peakRate_nd-lowRate_nd)./(lowRate_nd+peakRate_nd));
    end
    %%
    cellList = find(subsetCells);
    for iCell = cellList;
        allFields.n = [allFields.n normalise1var(placeFields_n(iCell,:))'];
        allFields.l = [allFields.l normalise1var(placeFields_l(iCell,:))'];
        allFields.h = [allFields.h normalise1var(placeFields_h(iCell,:))'];
    end
    numCells(animalIdx) = sum(subsetCells);
    numCells_decP(animalIdx) = sum(subsetCells_dec);
    numCells_low(animalIdx) = sum(subsetCells_low);
end

%%

clear p
for n = 1:length(hPost_h)
    [~,p(n)] = ttest2(hPost_h{n},hPost_l{n});
end
[hPost.h'-hPost.l' p']