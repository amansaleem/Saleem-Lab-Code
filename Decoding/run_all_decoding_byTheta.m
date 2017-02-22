% function run_all_decoding
idx = 0;
smthWin = 0;
box_filt = 1;
trainCorrect = 0;
quickProcess = 0;


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

for idx = [6 4 5 7 8 1:3] %1:length(expt_list)
   [out] = plot_some_decoding_byTheta(expt_list(idx).animal, ...
        expt_list(idx).iseries, expt_list(idx).expt_list, 'contrast', trainCorrect, smthWin, box_filt, quickProcess, expt_list(idx).shank_list);
   varName = ['Posterior_' num2str(expt_list(idx).iseries) '_smooth_' num2str(smthWin)];
   eval([varName ' = out;']);
   Posterior_all(idx) = eval(varName);
   %    close all;
   drawnow
   clear out
end
for idx = 1:length(expt_list)
    Posterior_all(idx).meta = expt_list(idx);
end

save Posterior_all_150623_trainedCV_thetaBins_noInter_100bestDec Posterior_all -v7.3
