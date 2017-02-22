% function run_all_decoding
idx = 0;
smthWin = 50;
box_filt = 1;
trainCorrect = 0;
quickProcess = 0;
% type = 'common'; 
type = 'indep';

for idx = 1:8
   es = [];
   load(['./Data/es_model_' type '_exp_' num2str(idx)]);
   [out] = plot_some_decoding_newer_sim(...
       es.animal, es.iseries, es.expt_list, es, ...
        'contrast', trainCorrect, smthWin, box_filt, quickProcess);
   varName = ['Posterior_' num2str(es.iseries) '_smooth_' num2str(smthWin)];
   eval([varName ' = out;']);
   Posterior_all(idx) = eval(varName);
%    close all; 
   
   
%    figure;
%    subplot(122)
%    plot(cumsum(hist(Posterior_all(idx).confidence.low,50))./sum(~isnan(Posterior_all(idx).confidence.low)))
%    hold on;
%    plot(cumsum(hist(Posterior_all(idx).confidence.high,50))./sum(~isnan(Posterior_all(idx).confidence.high)),'r')
%    title('Confidence')
%    subplot(121)
%    plot(cumsum(hist(Posterior_all(idx).error.high,50))./sum(~isnan(Posterior_all(idx).error.high)),'r')
%    hold on;
%    plot(cumsum(hist(Posterior_all(idx).error.low,50))./sum(~isnan(Posterior_all(idx).error.low)))
%    title('Error')
   pause(2)
   clear out
end
% for idx = 1:length(expt_list)
%     Posterior_all(idx).meta = es.expt_list(idx);
% end

save Posterior_all_150120_sim_indep_50 -v7.3
