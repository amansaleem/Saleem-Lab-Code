function expt = LoadContinuousData(animal,iseries,iexp)


%% initialization
global DIRS
global pepNEV

SetDefaultDirs

DIRname = [DIRS.Cerebus filesep,animal filesep num2str(iseries) filesep num2str(iexp)];
fname = [animal '_' num2str(iseries) '_' num2str(iexp)];
dname = [DIRname filesep fname '.ns5'];


%% open protocol
p = ProtocolLoad(animal,iseries,iexp);


%% open nev file and load data
[nevopen_outcome,~,~] = nevopen([dname(1:end-2) 'ev']);
if nevopen_outcome < 0
    error('ExptLoadBlackRock:NevOpenError','Could not open %s.',dname);
end

% photodiode times in units of clock ticks
stimOnTimes  = pepNEV.sync.timestamps(1:2:end);
stimOffTimes = pepNEV.sync.timestamps(2:2:end);

% error checking
if diff([length(stimOnTimes) length(stimOffTimes) max(p.seqnums(:))]) ~= 0
    error('umatched number of start/stop sync times or number of trials for %s', dname);
end


%% open .ns* file and load data
[~,SamplingRateInKHZ,nchan] = nsopen(dname);
for ichan = 1:nchan
    expt.data{ichan} = single(pepNEV.ns.Data.data(ichan,:));
end


%% get Michigan Layout and rearrange
[arrayLayout] = MichiganGetLayout(animal,iseries,nchan);
expt.data = expt.data(arrayLayout);
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
ballspeed = zeros(1,length(expt.data{1}));
ballspeed(datapoints) = movement;


%% write expt struct for output
expt.animal           = animal;
expt.iseries          = iseries;
expt.iexp             = iexp;
expt.data             = expt.data';
expt.stimOnTimes      = stimOnTimes;
expt.stimOffTimes     = stimOffTimes;
expt.ballspeed        = ballspeed;
expt.samplerate       = SamplingRateInKHZ * 1000;


%% save experiment file
cd('\\Zserver\Lab\Tmp\marieke\continuous data')
save([animal '_' int2str(iseries) '_' int2str(iexp) '_continuous'],'expt');
