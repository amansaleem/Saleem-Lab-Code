%... datafile001.ns5 = zirkus monitor
%... datafile002.ns5 = zupervision monitor
%... datafile003.ns5 = matrox monitor (blue)
%... datafile004.ns5 = grating in vs, tf=4 Hz
%... datafile005.ns5 = virtual reality, flickering photodiode square
%... datafile006.ns5 = virtual reality, white wall
%... datafile007.ns5 = table

global pepNEV

[~,SamplingRateInKHZ,nChans] = nsopen('\\ZSERVER\Data\Cerebus\MARIEKE_TEST\6\M110125_BALL_spont-refScrew.ns5');

for chan=1:16
data = single(pepNEV.ns.Data.data(chan,:));

data = data - mean(data);
WL = SamplingRateInKHZ*1000;
nO = round(WL*0.8);
[pow(chan,:),f] = pwelch(data,WL,nO,[],SamplingRateInKHZ*1000);

range = find(f<100);

end

figure;
for chan=1:16
semilogy(f(range),pow(chan,range),'Color',rand(1,3)); hold on
pause
end



