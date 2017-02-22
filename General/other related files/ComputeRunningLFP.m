function ComputeRunningLFP(animal,iseries,iexp,ichan)

SetDefaultDirs

% load data
p = ProtocolLoad(animal,iseries,iexp);
expt = ExptLoadCerebusTraces(animal,iseries,iexp);
balldata = getBallTuning(animal,iseries,iexp);


% determine LFP power spectrum of each stimulus
for istim = 1:p.nstim
    for irep = 1:p.nrepeats
        data = expt.data{ichan}{istim,irep} - mean(expt.data{ichan}{istim,irep});
        WindowLength = expt.samplerate/2; %0.5 stimulus-length
        nOverlap = round(WindowLength*0.8);
        [pow,f] = pwelch(data,WindowLength,nOverlap,[],expt.samplerate);
        power(istim,irep,:) = pow;
    end
end


% divide into running and no-running trials
run = power(repmat(balldata.ismoving,[1 1 size(power,3)])==1);
run = reshape(run,[sum(sum(balldata.ismoving)) size(power,3)]);

sit = power(repmat(balldata.ismoving,[1 1 size(power,3)])==0);
sit = reshape(sit,[p.nstim*p.nrepeats-sum(sum(balldata.ismoving)) size(power,3)]);

nrun = size(run,1);
nsit = size(sit,1);

meanrun = nanmean(run,1);
meansit = nanmean(sit,1);


% plot power spectra
figure 
semilogy(f',meanrun,'r','LineWidth',2), hold on
semilogy(f',meansit,'b','LineWidth',2)
legend('power when running','power when sitting')
xlabel('frequency (Hz)'), ylabel('power')