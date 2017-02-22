run_trials = find(balldata.ismoving);
stop_trials = find(~balldata.ismoving);

WL = 30000;
nO = round(WL*0.8);
% [pow_run,f_run] = pwelch(run,WL,nO,[],30000);
% [pow_stop,f_stop] = pwelch(stop,WL,nO,[],30000);
% f_stop
% range = find(abs(f-120) == min(abs(f-120)));

for n = 1:length(run_trials)
 run = expt.data{2}{run_trials(n)};
 run = run - mean(run);
 [pow_run(n,:),f] = pwelch(run,WL,nO,[],WL);
end

for n = 1:length(stop_trials)
 stop = expt.data{2}{stop_trials(n)};
 stop = stop - mean(stop);
 [pow_stop(n,:),f] = pwelch(stop,WL,nO,[],WL);
end

range = find(abs(f-100) == min(abs(f-100)));

plot(f(1:range), pow_stop(1:range), 'k', f(1:range), pow_run(1:range), 'r')