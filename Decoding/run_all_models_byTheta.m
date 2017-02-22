% function run_all_decoding
clear expt_list
idx = 0;

idx = idx + 1;
expt_list(idx).animal     = 'M130920_BALL';
expt_list(idx).iseries    = 1025;
expt_list(idx).expt_list  = 102:103;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M130918_BALL';
expt_list(idx).iseries    = 1030;
expt_list(idx).expt_list  = 103:105;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 530;
expt_list(idx).expt_list  = 104:106;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 531;
expt_list(idx).expt_list  = 103:106;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 601;
expt_list(idx).expt_list  = 103:106;
expt_list(idx).shank_list  = [0:3 5:7];

idx = idx + 1;
expt_list(idx).animal     = 'M140501_BALL';
expt_list(idx).iseries    = 602;
expt_list(idx).expt_list  = 102:106;
expt_list(idx).shank_list  = 0:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140502_BALL';
expt_list(idx).iseries    = 603;
expt_list(idx).expt_list  = 107:110;
expt_list(idx).shank_list  = 2:7;

idx = idx + 1;
expt_list(idx).animal     = 'M140502_BALL';
expt_list(idx).iseries    = 604;
expt_list(idx).expt_list  = 107:110;
expt_list(idx).shank_list  = [1 3 5 6];

smthWin = 5;
box_filt = 1;
trainCorrect = 0;
quickProcess = 1;
loadA = 0;
list = 8;

for idx = 1:length(expt_list)
   [out] = getModels(expt_list(idx).animal, ...
        expt_list(idx).iseries, trainCorrect, smthWin, box_filt, expt_list(idx).shank_list, loadA);
   varName = ['Posterior_' num2str(expt_list(idx).iseries) '_smooth_' num2str(smthWin)];
   eval([varName ' = out;']);
   out.meta = expt_list(idx);
   models(idx) = eval(varName);
   %    close all;
   drawnow
   clear out
end
% save Posterior_all_150706_trainedCorrect_thetaBins_FineBin1000_FixedSmth20 Posterior_all -v7.3
% save Posterior_all_160526_trainedCV_thetaBins_FineBin_EVover50_quick1 -v7.3
%%
% models = models(list);
thres = 0.001;
X = [0:0.001:0.04];
Y = [-0.02:0.001:0.02];
for n = 1:length(models)
    t = nanmean(models(n).norm.EV)>thres ...
        & nanmean(models(n).low.EV)>thres ...;
        & nanmean(models(n).high.EV)>thres...
        ;
    subplot(324)
    [sort_order] = plotPlaceFields_model(models(n).norm, t);
    subplot(323)
    hist(nanmean(models(n).norm.EV),X)
    sort_order = [];
    subplot(322)
    plotPlaceFields_model(models(n).low, t, sort_order);
    subplot(321)
    hist(nanmean(models(n).low.EV),X)
    
    subplot(326)
    plotPlaceFields_model(models(n).high, t, sort_order);
    subplot(325)
    hist(nanmean(models(n).high.EV),X)
    
    %     subplot(222)
    %     hist(nanmean(models(n).norm.EV)-nanmean(models(n).low.EV),Y)
    %     hold on;
    %     y = nanmean(models(n).norm.EV(:) - nanmean(models(n).low.EV(:)));
    %     line([y y], ylim);
    %     hold off;
    %     subplot(224)
    %     hist(nanmean(models(n).high.EV)-nanmean(models(n).low.EV),Y)
    %     hold on;
    %     y = nanmean(models(n).high.EV(:) - nanmean(models(n).low.EV(:)));
    %     line([y y], ylim);
    h = nanmean(nanmean(models(n).high.EV));
    l = nanmean(nanmean(models(n).low.EV));
    m = nanmean(nanmean(models(n).norm.EV));
    display([num2str(n) '   ' num2str(loadA)])
    [l m h]
    pause
    %     hold off
end