function [balldata] = getNewBallData(animal, iseries, iexp)
% MS 2010-07-21 created
% AS 2010-12-10 two major changes made. File loading and OnTimes and
%               OffTimes detection amended to be slightly more general and not affected
%               by photodiode level
% AA 2010-12-14
% AS 2010-12-18 Added event detection and loading of existing ball data files, rather than processing
% EVENT DETECTION CAN BECOME MAIN

SetDefaultDirs

AnimalDir = fullfile(DIRS.ball,'Recording',animal);
if ~isdir(AnimalDir), mkdir(AnimalDir); end

SeriesDir = fullfile(AnimalDir,num2str(iseries));
if ~isdir(SeriesDir), mkdir(SeriesDir); end

ExpDir    = fullfile(SeriesDir,num2str(iexp));
if ~isdir(ExpDir), mkdir(ExpDir); end

FileName  = fullfile(ExpDir, [num2str(animal) '_' num2str(iseries) '_' num2str(iexp) '_BallData.mat']);
% AS: to avoid re-processing each time... load existing file
if exist(FileName)==2
    fprintf('Loading existing Ball data file ...')
    load(FileName);
    fprintf('done\n');
    return
end

% load data
DIRname = [DIRS.Cerebus filesep,animal filesep num2str(iseries) filesep num2str(iexp) filesep];
fname = [animal '_' num2str(iseries) '_' num2str(iexp)];
load([DIRname filesep fname '_photodiode']);
load([DIRname filesep fname '_ballCoords']);
fprintf('Processing ball data...');
prot = ProtocolLoad(animal,iseries,iexp);

% find start times of stimuli (negative deflection in photodiode signal)
OnThres = min(photodiode.values)/5; % AS: added this to avoid problems when changing the photodiode level
OnTimes = find(([0 diff(photodiode.values)']')<=OnThres & (photodiode.values)<=OnThres);
idx = find(diff(OnTimes)<1000)+1; % AS: added this to avoid multiple detecting within short periods
OnTimes(idx)=[];

% find stop times of stimuli (positive deflection in photodiode signal)
OffThres = max(photodiode.values)/2;
OffTimes = find(([0 diff(photodiode.values)']')>=OffThres & (photodiode.values)>=OffThres);
idx = find(diff(OffTimes)<1000)+1;
OffTimes(idx)=[];

if diff([length(OnTimes) length(OffTimes)]) ~= 0
    clear OnTimes OffTimes
    % AS: To detect events, rather than On and Off times, similar to a digital event detection
    eventThres = max(abs(photodiode.values))/4.5;
    events = find(([0 abs(diff((photodiode.values)))']')>eventThres & abs(photodiode.values)>eventThres);
    idx = find(diff(events)<1000)+1;
    events(idx)=[];
    OnTimes = events(1:2:end);
    OffTimes = events(2:2:end);
    if diff([length(OnTimes) length(OffTimes)]) ~= 0
        error('unmatched number of stimulus on- and off-times for photodiode signal');
    end
end

% align photodiode and ball traces
startball = str2num(XY{1}((end-8):end));
timediff = startball - photodiode.starttime(2); % in sec
photoSampleRatio = (photodiode.stoptime(2)-photodiode.starttime(2))./length(photodiode.values); %(prot.pars(1)/10)./mean(OffTimes-OnTimes);
% write ball coordinates in different format
ballcoords = zeros(6,length(XY));
for tpt = 1:length(XY)
    if ~isempty(XY{tpt})
        space = find(isspace(XY{tpt})==1);
        ballcoords(1,tpt) = str2num(XY{tpt}((space(1)+1):space(2)));
        ballcoords(2,tpt) = str2num(XY{tpt}((space(2)+1):space(3)));
        ballcoords(3,tpt) = str2num(XY{tpt}((space(3)+1):space(4)));
        ballcoords(4,tpt) = str2num(XY{tpt}((space(4)+1):space(5)));
        ballcoords(5,tpt) = str2num(XY{tpt}((space(5)):end));
        ballcoords(6,tpt) = str2num(XY{tpt}(1:(space(1))));
    else
        break
    end
end
ballcoords = ballcoords(:,1:(tpt-1));

% get snippets of ball data
for stim = 1:length(OnTimes)
    stimOnBall = timediff+(OnTimes(stim).*photoSampleRatio);
    stimOffBall = timediff+(OffTimes(stim).*photoSampleRatio);
    [temp istimOn] = min(abs(ballcoords(5,:)-stimOnBall));
    [temp istimOff] = min(abs(ballcoords(5,:)-stimOffBall));
    balldata.data{stim} = ballcoords(:,istimOn:istimOff);
end

% save ball data in same directory as 'expt' structure
%save(FileName, 'balldata');
%fprintf(' done\n');