function [expt2 VRdata chans] = getVRspikeTimes(animal,series,iexp)

%% initialization
global DIRS
global REPLAY
% global pepNEV

SetDefaultDirs

if ~REPLAY
    fname = [animal '_' num2str(series) '_' num2str(iexp)];
    dDIRname = [DIRS.data filesep,animal filesep num2str(series) filesep num2str(iexp)];
else
    fname = [animal '_' num2str(series) '_R' num2str(iexp)];
    dDIRname = [DIRS.data filesep,animal filesep num2str(series) filesep 'R' num2str(iexp)];    
end

load([dDIRname filesep fname '_contspikes']);

%% get all spikeTimes
for ichan = 1:length(chans)
    for iprot = 1:chans(ichan).nprots
        n = 1;
        for ispike = 1:length(chans(ichan).spikes)
            if chans(ichan).cc(ispike) == iprot
                expt2.spikeTimes{ichan,iprot}(n) = chans(ichan).spikes(ispike).t;
                n = n + 1;
            end
        end
        if n == 1
            expt2.spikeTimes{ichan,iprot} = [];
            expt2.spikeRate{ichan,iprot} = [];
        else
            for iscreen = 2:length(expt2.screenTimes)
                nspike = sum((expt2.spikeTimes{ichan,iprot} > expt2.screenTimes(iscreen-1)) & (expt2.spikeTimes{ichan,iprot} <= expt2.screenTimes(iscreen)));
                expt2.spikeRate{ichan,iprot}(iscreen-1) = nspike/(expt2.screenTimes(iscreen)-expt2.screenTimes(iscreen-1));
            end
        end
    end
end

spikeRate = expt2.spikeRate;
% save them to expt
save([dDIRname filesep fname '_contspikes'],'expt2','VRdata','chans');
% a = 1
%% get spikerate

