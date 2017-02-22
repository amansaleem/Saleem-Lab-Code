%%
high_distr = (max(Posterior_all.Posterior_high')-min(Posterior_all.Posterior_high'));
low_distr = (max(Posterior_all.Posterior_low')-min(Posterior_all.Posterior_low'));
mean(high_distr) - mean(low_distr)
[h,p] = kstest2(high_distr, low_distr)
%%
% 
% maxInTime_n = max(P_918.Posterior_norm');
% maxInTime_l = max(P_918.Posterior_low');
% maxInTime_h = max(P_918.Posterior_high');
% minInTime_n = min(P_918.Posterior_norm');
% minInTime_l = min(P_918.Posterior_low');
% minInTime_h = min(P_918.Posterior_high');
% 
% confidence_l = maxInTime_l - minInTime_l;
% confidence_n = maxInTime_n - minInTime_n;
% confidence_h = maxInTime_h - minInTime_h;
% 
% for t = 1:size(P_918.Posterior_norm,1)
%     accuracy_n(t) = P_918.Posterior_norm(t,P_918.X_norm(t));
% end
% for t = 1:size(P_918.Posterior_high,1)
%     accuracy_h(t) = P_918.Posterior_high(t,P_918.X_high(t));
% end
% for t = 1:size(P_918.Posterior_low,1)
%     accuracy_l(t) = P_918.Posterior_low(t,P_918.X_low(t));
% end
%    
% precision_n = std(P_918.Posterior_norm');
% precision_l = std(P_918.Posterior_low');
% precision_h = std(P_918.Posterior_high');
% 
% figure
% %%
% [conf_distr,bins]    = hist(confidence_n,20);
% [conf_distr_h]       = hist(confidence_h,bins);
% [conf_distr_l]       = hist(confidence_l,bins);
% 
% conf_distr = conf_distr./sum(conf_distr);
% conf_distr_l = conf_distr_l./sum(conf_distr_l);
% conf_distr_h = conf_distr_h./sum(conf_distr_h);
% subplot(311)
% plot(bins,(conf_distr),'k','linewidth',1.5)
% hold on;
% plot(bins,(conf_distr_l),'b','linewidth',1.5)
% plot(bins,(conf_distr_h),'r','linewidth',1.5)
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight
% axis square
% xlabel('Confidence')
% %%
% [accu_distr,bins]    = hist(accuracy_n,20);
% [accu_distr_h]       = hist(accuracy_h,bins);
% [accu_distr_l]       = hist(accuracy_l,bins);
% 
% accu_distr = accu_distr./sum(accu_distr);
% accu_distr_l = accu_distr_l./sum(accu_distr_l);
% accu_distr_h = accu_distr_h./sum(accu_distr_h);
% subplot(312)
% plot(bins,(accu_distr),'k','linewidth',1.5)
% hold on;
% plot(bins,(accu_distr_l),'b','linewidth',1.5)
% plot(bins,(accu_distr_h),'r','linewidth',1.5)
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis tight;
% axis square
% xlabel('Accuracy')
% line([0 0],ylim)
% %%
% [prec_distr,bins]    = hist(precision_n,20);
% [prec_distr_h]       = hist(precision_h,bins);
% [prec_distr_l]       = hist(precision_l,bins);
% 
% prec_distr = prec_distr./sum(prec_distr);
% prec_distr_l = prec_distr_l./sum(prec_distr_l);
% prec_distr_h = prec_distr_h./sum(prec_distr_h);
% subplot(313)
% plot(bins,cumsum(accu_distr),'k','linewidth',1.5)
% hold on;
% plot(bins,cumsum(accu_distr_l),'b','linewidth',1.5)
% plot(bins,cumsum(accu_distr_h),'r','linewidth',1.5)
% set(gca, 'box','off','TickDir','out','fontsize',14,'color','none')
% axis square
% xlabel('Precision');
