function getBallData_try(animal, iseries, iexp)
% MS 2010-07-21 created
% AS 2010-12-10 two major changes made. File loading and OnTimes and
%               OffTimes detection amended to be slightly more general and not affected
%               by photodiode level

SetDefaultDirs

% load data
% AS: commented this out and added the lower two lines because of problems loading the Mouse files
% fname = sprintf('%s%su%03i_%03i',[DIRS.Cerebus filesep],[animal filesep num2str(iseries) filesep num2str(iexp) filesep],iseries,iexp);
DIRname = [DIRS.Cerebus filesep,animal filesep num2str(iseries) filesep num2str(iexp) filesep];
fname = [animal '_' num2str(iseries) '_' num2str(iexp)];
load([DIRname filesep fname '_photodiode']);
load([DIRname filesep fname '_ballCoords']);

% find start times of stimuli (negative deflection in photodiode signal)
OnThres = min(photodiode.values)/5; % AS: added this to avoid problems when changing the photodiode level
OnTimes = find(([0 diff(photodiode.values)']')<=OnThres & (photodiode.values)<=OnThres);
idx = find(diff(OnTimes)<100)+1; % AS: added this to avoid multiple detecting within short periods
OnTimes(idx)=[];

% find stop times of stimuli (positive deflection in photodiode signal)
OffThres = max(photodiode.values)/5;
OffTimes = find(([0 diff(photodiode.values)']')>=OffThres & (photodiode.values)>=OffThres);
idx = find(diff(OffTimes)<100)+1;
OffTimes(idx)=[];

if diff([length(OnTimes) length(OffTimes)]) ~= 0
    error('unmatched number of stimulus on- and off-times for photodiode signal');
end

% align photodiode and ball traces
startball = str2num(XY{1}((end-8):end));
timediff = (startball - photodiode.starttime)*1000; % in ms
sampleRatio = length(XY) / length(photodiode.values);

% write ball coordinates in different format
for tpt = 1:length(XY)
    space = find(isspace(XY{tpt})==1);
    ballcoords(1,tpt) = str2num(XY{tpt}((space(1)+1):space(2)));
    ballcoords(2,tpt) = str2num(XY{tpt}((space(2)+1):space(3)));
    ballcoords(3,tpt) = str2num(XY{tpt}((space(3)+1):space(4)));
    ballcoords(4,tpt) = str2num(XY{tpt}((space(4)+1):space(5)));
end

% get snippets of ball data
for stim = 1:length(OnTimes)
    stimonball = round((OnTimes(stim)+timediff)*sampleRatio);
    stimoffball = round((OffTimes(stim)+timediff)*sampleRatio);
    balldata{stim} = ballcoords(:,stimonball:stimoffball);
end

% save ball data in same directory as 'expt' structure
AnimalDir = fullfile(DIRS.data,animal);
SeriesDir = fullfile(AnimalDir,num2str(iseries));
ExpDir    = fullfile(SeriesDir,num2str(iexp));
FileName  = fullfile(ExpDir, [num2str(animal) '_' num2str(iseries) '_' num2str(iexp) '_BallData.mat']);
save(FileName, 'balldata');
