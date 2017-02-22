function [expt,VRdata,chans] = LoadContinuousVRData(animal,series,iexp,replay_in,flag_a)


%% initialization
global DIRS
global pepNEV
global REPLAY
% profile on
SetDefaultDirs


if isempty(REPLAY)
    if nargin<4
        REPLAY = 0;
    else
        REPLAY = replay_in;
    end
end

if nargin<5
    flag_a = 0;
end

if ~REPLAY
    DIRname = [DIRS.Cerebus filesep,animal filesep num2str(series) filesep num2str(iexp)];
    fname = [animal '_' num2str(series) '_' num2str(iexp)];
    dname = [DIRname filesep fname '.ns5'];
    
    dDIRname = [DIRS.data filesep,animal filesep num2str(series) filesep num2str(iexp)];
    sDIRname = [DIRS.spikes filesep,animal filesep num2str(series)];
else
    DIRname = [DIRS.Cerebus filesep,animal filesep num2str(series) filesep 'R' num2str(iexp)];
    fname = [animal '_' num2str(series) '_R' num2str(iexp)];
    dname = [DIRname filesep fname '.ns5'];
    
    dDIRname = [DIRS.data filesep,animal filesep num2str(series) filesep 'R' num2str(iexp)];
    sDIRname = [DIRS.spikes filesep,animal filesep num2str(series)];
end

if exist([dDIRname filesep fname '_contspikes.mat'])
    if ~flag_a
        load([dDIRname filesep fname '_continuous'])
    end
    load([dDIRname filesep fname '_contspikes'])
    try
        load([sDIRname filesep 'DiscrimParsCerebusTraces_' num2str(iexp)])
        if ~flag_a
            save([dDIRname filesep fname '_continuous'],'expt','VRdata','chans');
        end
        save([dDIRname filesep fname '_contspikes'],'expt2','VRdata','chans');
    catch
        display('Cannot find discrimination parameters');
    end
else
    if ~REPLAY
        VRDIRname = [DIRS.ball filesep,animal filesep num2str(series)];
        VRfname = [animal '_' num2str(series) '_session_' num2str(iexp) '_trial001'];
    else
        VRDIRname = [DIRS.ball filesep,animal filesep num2str(series) filesep '_R'];
        VRfname = [animal '_' num2str(series) '_session_' num2str(iexp) 'R_trial001'];
    end
    
    VRdname = [VRDIRname filesep VRfname '.mat'];
    
    VRdata = load(VRdname);
    
    %% open protocol
    % p = ProtocolLoad(animal,series,iexp);
    
    
    %% open nev file and load data
    [nevopen_outcome,~,~] = nevopen([dname(1:end-2) 'ev']);
    if nevopen_outcome < 0
        error('ExptLoadBlackRock:NevOpenError','Could not open %s.',dname);
    end
    
    % photodiode times in units of clock ticks
    screenTimes  = pepNEV.sync.timestamps;
    screenTimes(find(diff(screenTimes)<100) + 1) = [];
    screenTimes(1) = [];
    
    % stimOffTimes = pepNEV.sync.timestamps(2:2:end);
    
    % error checking
    % if diff([length(stimOnTimes) length(stimOffTimes) max(p.seqnums(:))]) ~= 0
    %     error('umatched number of start/stop sync times or number of trials for %s', dname);
    % end
    
    
    %% open .ns* file and load data
    [~,SamplingRateInKHZ,nchan] = nsopen(dname);
    for ichan = 1:nchan
        data{ichan} = single(pepNEV.ns.Data.data(ichan,:));
    end
    
    
    %% get Michigan Layout and rearrange
    [arrayLayout] = MichiganGetLayout(animal,series,nchan);
    data = data(arrayLayout);
    nevclose;
    
    
    %% get speed of ball
    load([DIRname filesep fname '_ballCoords']);
    
    % write ball coordinates in different format
    ballcoords = zeros(5,length(XY));
    for tpt = 1:length(XY)
        if ~isempty(XY{tpt})
            space = find(isspace(XY{tpt})==1);
            ballcoords(1,tpt) = str2num(XY{tpt}((space(1)+1):space(2)));
            ballcoords(2,tpt) = str2num(XY{tpt}((space(2)+1):space(3)));
            ballcoords(3,tpt) = str2num(XY{tpt}((space(3)+1):space(4)));
            ballcoords(4,tpt) = str2num(XY{tpt}((space(4)+1):space(5)));
            ballcoords(5,tpt) = str2num(XY{tpt}((space(5)):end));
        else
            break
        end
    end
    ballcoords = ballcoords(:,1:(tpt-1));
    
    % determine ball speed
    movement = sqrt(mean(ballcoords(1:4,:),1).^2);
    
    % write 'movement' in sample rate of neural data
    ballFs = median(diff(ballcoords(5,:)));             % sample rate of ball
    freq = round(ballFs * SamplingRateInKHZ * 1000);    % how often the ball is sampled compared to neural data
    datapoints = [freq:freq:freq*size(ballcoords,2)];   % write vector of datapoints in neural data were ball data were sampled
    ballspeed = zeros(1,length(data{1}));
    ballspeed(datapoints) = movement;
    
    % Write the times of the screen refresh
    screenRefresh = zeros(1,length(data{1}));
    screenRefresh(screenTimes) = 1;
    
    %% ns5 file, originally sampled at 30000 s/sec
    filter_bounds = [500 7000];
    [b,a] = butter(3,2*filter_bounds/SamplingRateInKHZ/1000);
    
    %%
    for n = 1:length(data)
        chanNames(n) = (500+n);
        threshvalues(n) = NaN;
        threshsigns(n) = -1;
        expt.data{n,1}{1} = data{1,n};
    end
    
    %% write expt struct for output
    expt.animal           = animal;
    expt.iseries          = series;
    expt.iexp             = iexp;
    expt.nrepeats         = 1;
    expt.nchans           = n;
    expt.Source           = 'CerebusTraces';
    expt.channames        = chanNames;
    expt.nstim             = 1;
    expt.threshsigns      = threshsigns;
    expt.threshvalues     = threshvalues;
    expt.filtercoeffs     = {b,a};
    expt.screenRefresh    = screenRefresh;
    expt.ballspeed        = ballspeed;
    expt.samplerate       = SamplingRateInKHZ * 1000;
    expt.screenTimes      = double(screenTimes)/expt.samplerate;
    
    expt.stimdurs         = NaN; % p.pfiledurs;
    expt.prestimdur       = 0;
    expt.poststimdur      = 0;
    expt.timestamp        = 'xxx';
    expt.unitspervolt     = NaN;
    expt.ChanDescriptions = {};
    expt.timestamps       = {};
    expt.ecg              = [];
    expt.eeg              = [];
    
    expt2.animal           = animal;
    expt2.iseries          = series;
    expt2.iexp             = iexp;
    expt2.nrepeats         = 1;
    expt2.nchans           = n;
    expt2.Source           = 'CerebusTraces';
    expt2.channames        = chanNames;
    expt2.nstim             = 1;
    expt2.threshsigns      = threshsigns;
    expt2.threshvalues     = threshvalues;
    expt2.filtercoeffs     = {b,a};
    expt2.screenRefresh    = screenRefresh;
    expt2.ballspeed        = ballspeed;
    expt2.samplerate       = SamplingRateInKHZ * 1000;
    expt2.screenTimes      = double(screenTimes)/expt.samplerate;
    
    expt2.stimdurs         = NaN; % p.pfiledurs;
    expt2.prestimdur       = 0;
    expt2.poststimdur      = 0;
    expt2.timestamp        = 'xxx';
    expt2.unitspervolt     = NaN;
    expt2.ChanDescriptions = {};
    expt2.timestamps       = {};
    expt2.ecg              = [];
    expt2.eeg              = [];
    
    clear global pepNEV
    try
        if ~REPLAY
            load([sDIRname filesep 'DiscrimParsCerebusTraces_' num2str(iexp)]);
        else
            load([sDIRname filesep 'DiscrimParsCerebusTraces_R' num2str(iexp)]);
        end
        
        if ~exist([DIRS.data filesep animal]); mkdir([DIRS.data filesep animal]); end;
        if ~exist([DIRS.data filesep,animal filesep num2str(series)]); mkdir([DIRS.data filesep,animal filesep num2str(series)]); end;
        
        if~REPLAY
                if ~exist([DIRS.data filesep,animal filesep num2str(series) filesep num2str(iexp)]); mkdir([DIRS.data filesep,animal filesep num2str(series) filesep num2str(iexp)]); end
        else
                if ~exist([DIRS.data filesep,animal filesep num2str(series) filesep num2str(iexp)]); mkdir([DIRS.data filesep,animal filesep num2str(series) filesep 'R' num2str(iexp)]); end
        end
        
        save([dDIRname filesep fname '_continuous'],'expt','VRdata','chans');
        save([dDIRname filesep fname '_contspikes'],'expt2','VRdata','chans');
    catch   
        %% save experiment file
        % cd('C:\Marieke\analysis mouse\CerebusData')
        if ~exist([DIRS.data filesep animal]); mkdir([DIRS.data filesep animal]); end;
        if ~exist([DIRS.data filesep,animal filesep num2str(series)]); mkdir([DIRS.data filesep,animal filesep num2str(series)]); end;
        
        if ~REPLAY
             if ~exist([DIRS.data filesep,animal filesep num2str(series) filesep num2str(iexp)]); mkdir([DIRS.data filesep,animal filesep num2str(series) filesep num2str(iexp)]); end
        else
            if ~exist([DIRS.data filesep,animal filesep num2str(series) filesep 'R' num2str(iexp)]); mkdir([DIRS.data filesep,animal filesep num2str(series) filesep 'R' num2str(iexp)]); end
        end
        
        save([dDIRname filesep fname '_continuous'],'expt','VRdata');
        save([dDIRname filesep fname '_contspikes'],'expt2','VRdata');
    end
    
end
% profile viewer
% p = profile('info');
% profsave(p, 'profile_results');
% 
% profile off
 
clear global REPLAY
