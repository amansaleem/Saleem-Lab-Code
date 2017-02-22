% function run_all_decoding
idx = 0;
smthWin = 100;
box_filt = 1;
trainCorrect = 0;
numTheta = 8;
trainPhases = 1:8;

% [Posterior_all] = popThetaPrecisionplot_some_decoding_new(animal, iseries, exp_list, type, only_correct, smthWin, box_filt, shank_list, numTheta)

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
expt_list(idx).numTheta    = numTheta;

for idx = 1:length(expt_list)
   [out] = popThetaPrecision(expt_list(idx).animal, ...
        expt_list(idx).iseries, expt_list(idx).expt_list, 'contrast', trainCorrect, smthWin, box_filt, expt_list(idx).shank_list, numTheta, trainPhases);
%     [Posterior_all] = popThetaPrecisionplot_some_decoding_new(animal, iseries, exp_list, type, only_correct, smthWin, box_filt, shank_list, numTheta)
   varName = ['Posterior_' num2str(expt_list(idx).iseries) '_smooth_' num2str(smthWin)];
   eval([varName ' = out;']);
   Posterior_all(idx) = eval(varName);
   drawnow; pause(5)
   clear out
end
for idx = 1:length(expt_list)
    Posterior_all(idx).meta = expt_list(idx);
end
save ./Data/Posterior_popTheta_150106_trainedCV_100 -v7.3
