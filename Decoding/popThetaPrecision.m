function [Posterior_all] = popThetaPrecision(...
    animal, iseries, exp_list, type, only_correct, smthWin, box_filt, shank_list, numTheta, train_bins, area)

%% Parameters that can be changed
test_normal = 1; % Learning/training on normal condition
plot_all = 1;

display(['Smoothing window = ' num2str(smthWin)]);
if nargin<4
    type = 'contrast';
end
if nargin<5
    only_correct = false;
end
if nargin<6
    smthWin = 250;
end
if nargin<7
    box_filt = 0;
end
if nargin<8
    shank_list = false;
end
if nargin<9
    numTheta = 8;
end
if nargin<10
    trainBins = 1:train_bins;
end
if nargin<11
    area = 'CA1';
end

switch animal
    case 'M130918_BALL'
        iseries = 1030;
        switch area
            case 'V1'
                cell_list = 190:321;
            case 'CA1'
                cell_list = 1:189;
        end
    case 'M130920_BALL'
        iseries = 1025;
        switch area
            case 'V1'
                cell_list = 153:249;
            case 'CA1'
                cell_list = 1:152;
        end
    otherwise
        cell_list = [];
end

%% Load the spiking and theta data
% SetDirs;
es = [];
% VRLoadMultipleExpts(animal, iseries, exp_list,'SPIKES_THETA',[18 22],shank_list);
% if ~isempty(cell_list)
%     es.spikeTrain = es.spikeTrain(:,cell_list);
% end
load(['Data' filesep 'es_' animal '_' num2str(iseries)]);
display(['Processing: ' animal '_' num2str(iseries)]);

%% Discretising the theta phase
[thetaPhase, thetaBins] = normalise1var(angle(es.theta.B.hill),numTheta);

%% Getting conditions

Posterior_all.animal = animal;
Posterior_all.iseries = iseries;
Posterior_all.exp_list = exp_list;
Posterior_all.smthWin = smthWin;
Posterior_all.only_correct = only_correct;
Posterior_all.box_filt = box_filt;
Posterior_all.numTheta = numTheta;
Posterior_all.train_bins = train_bins;

Posterior_all.thetaPhase = thetaPhase;
Posterior_all.thetaBins = thetaBins;
Posterior_all.thetaPow   = abs(es.theta.B.hill);

trainPhases = zeros(size(es.traj));
for iTrainBins = train_bins
    trainPhases(Posterior_all.thetaPhase==iTrainBins) = 1;
end

switch type
    case 'contrast'
        base = es.traj~=0 & es.contrast~=0 & ~isnan(es.traj) & es.gain==1 & es.roomLength==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.contrast==0.6 & base;
        t_low  = es.contrast==0.18 & base;
        t_high = es.contrast==0.72 & base;
        t_gray = es.traj~=0 & es.contrast==0 & ~isnan(es.traj) & es.smthBallSpd>5;
    case 'gain'
        base = es.traj~=0 & ~isnan(es.traj) & es.contrast==0.6 & es.roomLength==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.gain==1 & base;
        t_low  = es.gain<1 & base;
        t_high = es.gain>1 & base;
    case 'roomlength'
        base = es.traj~=0 & ~isnan(es.traj) & es.contrast==0.6 & es.gain==1 & es.smthBallSpd>5 & es.trajspeed>=0;
        t      = es.roomLength==1 & base;
        t_low  = es.roomLength<1 & base;
        t_high = es.roomLength>1 & base;
end

spkRate = zeros(size(es.spikeTrain));
es.spikeTrain(es.traj==0,:) = nan; % to see if there is an effect of the gray screen.

for icell = 1:size(es.spikeTrain,2);
    if box_filt
        spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 60, smthWin,'box');
    else
        spkRate(:,icell) = smthInTime(es.spikeTrain(:,icell), 60, smthWin); %*(smthWin*60./1000);
    end
end

%% Decoding all the conditions
dec = bayesDecoder;
dec.numBins = 50;

dec_l = bayesDecoder;%dec;
dec_l.numBins = 50;
dec_h = bayesDecoder;%dec;
dec_h.numBins = 50;

if only_correct
    dec2 = bayesDecoder;%dec;
    dec2.numBins = 50;
    [dec, ~, ~, ~] = ...
        dec.trainDecoder( es.traj(t & trainPhases & es.outcome==2), spkRate(t & es.outcome==2 & trainPhases,:), 0);
    [~,  Posterior_norm] = dec.predictBayesDecoder(spkRate(t,:), 0,'best');
    [dec2, ~, X_norm, ~] = ...
        dec2.trainDecoder( es.traj(t), spkRate(t,:), 0);
else
    dec2 = bayesDecoder;%dec;
    dec2.numBins = 50;
    [dec, ~, ~, ~] = ...
        dec.trainDecoder( es.traj(t & trainPhases),...
        spkRate(t & trainPhases,:), 0);
    [~,  Posterior_norm] = dec.predictBayesDecoder(spkRate(t,:), 0,'best');
    [dec2, ~, X_norm, ~] = ...
        dec2.trainDecoder( es.traj(t), spkRate(t,:), 0);
end

% [dec_l, ~, X_low_orig, Posterior_low_orig] = ...
%     dec_l.trainDecoder( es.traj(t_low), spkRate(t_low,:), 0);
% [dec_h, ~, X_high_orig, Posterior_high_orig] = ...
%     dec_h.trainDecoder( es.traj(t_high), spkRate(t_high,:), 0);
[ML_low,  Posterior_low] = dec.predictBayesDecoder(spkRate(t_low,:), 0,'best');
[~,  Posterior_gray] = dec.predictBayesDecoder(spkRate(t_gray,:), 0,'best');
[ML_high,  Posterior_high] = dec.predictBayesDecoder(spkRate(t_high,:), 0,'best');

%% Cleaning up the sizes (to account for cross-validation division rounding)
[X_low, bins_low] = normalise1var(es.traj(t_low), dec.numBins);
[X_high,bins_high] = normalise1var(es.traj(t_high), dec.numBins);

% X_high_orig = X_high(1:size(Posterior_high_orig,1));
% X_low_orig  = X_low(1:size(Posterior_low_orig,1));
X_norm = X_norm(1:size(Posterior_norm,1));

timePs = find(t);
timePs = timePs(1:size(Posterior_norm,1));
t = false(size(t));
t(timePs) = true;

timePs = find(t_low);
timePs = timePs(1:size(Posterior_low,1));
t_low = false(size(t_low));
t_low(timePs) = true;

timePs = find(t_high);
timePs = timePs(1:size(Posterior_high,1));
t_high = false(size(t_high));
t_high(timePs) = true;

%% Defining the outputs
Posterior_all.Posterior_norm = Posterior_norm;
Posterior_all.Posterior_high = Posterior_high;
Posterior_all.Posterior_low  = Posterior_low ;

Posterior_all.Posterior_gray  = Posterior_gray;
Posterior_all.only_correct = only_correct;

Posterior_all.X_norm = X_norm;
Posterior_all.X_high = X_high;
Posterior_all.X_low  = X_low;

Posterior_all.data      = es;
Posterior_all.decoder   = dec;
% Posterior_all.decoder_low   = dec_l;
% Posterior_all.decoder_high   = dec_h;

Posterior_all.t_norm    = t;
Posterior_all.t_low     = t_low;
Posterior_all.t_high    = t_high;

%% Getting the adjusted confusion matrixes
Posterior_diff_norm = nan*ones(size(Posterior_norm,1),99);
Posterior_diff_low  = nan*ones(size(Posterior_low,1),99);
Posterior_diff_high = nan*ones(size(Posterior_high,1),99);

for idx = 1:length(X_norm)
    Posterior_diff_norm(idx,(51-X_norm(idx)):(100-X_norm(idx))) = Posterior_norm(idx,:);
end
for idx = 1:length(X_low)
    Posterior_diff_low(idx,(51-X_low(idx)):(100-X_low(idx))) = Posterior_low(idx,:);
end
for idx = 1:length(X_high)
    Posterior_diff_high(idx,(51-X_high(idx)):(100-X_high(idx))) = Posterior_high(idx,:);
end

Posterior_all.Posterior_diff_norm = Posterior_diff_norm;
Posterior_all.Posterior_diff_low  = Posterior_diff_low;
Posterior_all.Posterior_diff_high = Posterior_diff_high;

%% Calculating the error, confidence and Accuracy
t_n = zeros(size(Posterior_all.Posterior_norm));
t_l = zeros(size(Posterior_all.Posterior_low));
t_h = zeros(size(Posterior_all.Posterior_high));
for n = 1:size(t_l,1)
    t_l(n,Posterior_all.X_low(n)) = 1;
end
for n = 1:size(t_h,1)
    t_h(n,Posterior_all.X_high(n)) = 1;
end
for n = 1:size(t_n,1)
    t_n(n,Posterior_all.X_norm(n)) = 1;
end
Posterior_all.confidence.low = (max(Posterior_all.Posterior_low') - min(Posterior_all.Posterior_low'));
Posterior_all.confidence.norm = (max(Posterior_all.Posterior_norm')- min(Posterior_all.Posterior_norm'));
Posterior_all.confidence.high = (max(Posterior_all.Posterior_high') - min(Posterior_all.Posterior_high'));
Posterior_all.confidence.gray = (max(Posterior_all.Posterior_gray') - min(Posterior_all.Posterior_gray'));

[~,X_ML_low] = max(Posterior_all.Posterior_low');
[~,X_ML_norm] = max(Posterior_all.Posterior_norm');
[~,X_ML_high] = max(Posterior_all.Posterior_high');
[~,X_ML_gray] = max(Posterior_all.Posterior_gray');

Posterior_all.MAP.low = X_ML_low;
Posterior_all.MAP.norm = X_ML_norm;
Posterior_all.MAP.high = X_ML_high;
Posterior_all.MAP.gray = X_ML_gray;

Posterior_all.error.low  = ((X_ML_low - X_low'));
Posterior_all.error.norm = ((X_ML_norm - X_norm'));
Posterior_all.error.high = ((X_ML_high - X_high'));

%% Discretising the theta phase
[thetaPhaseB, thetaBinsB] = normalise1var(angle(Posterior_all.data.theta.B.hill),numTheta);
[thetaPhase, thetaBins] = normalise1var(angle(Posterior_all.data.theta.A.hill),numTheta);

Posterior_all.thetaPhase = thetaPhase;
Posterior_all.thetaPhaseB = thetaPhaseB;

Posterior_all.thetaBins = thetaBins;
Posterior_all.thetaPow   = abs(Posterior_all.data.theta.A.hill);
Posterior_all.thetaPowB   = abs(Posterior_all.data.theta.B.hill);

%% Getting the confusion matrixes
for n = 1:dec.numBins
    %
    timePts = (X_norm == n);
    meanPost_norm(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_high == n);
    meanPost_high(:,n) = nanmean(Posterior_high(timePts, :),1);
    %
    timePts = (X_low == n);
    meanPost_low(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
%     timePts = (X_high_orig == n);
%     meanPost_high_orig(:,n) = nanmean(Posterior_high_orig(timePts, :),1);
%     %
%     timePts = (X_low_orig == n);
%     meanPost_low_orig(:,n) = nanmean(Posterior_low_orig(timePts, :),1);
    %
    timePts = (X_norm == n & (es.outcome(t) ~= 2));
    meanPost_norm_pass(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_norm == n & (es.outcome(t) == 2));
    meanPost_norm_actv(:,n) = nanmean(Posterior_norm(timePts, :),1);
    %
    timePts = (X_low == n & (es.outcome(t_low) == 2));
    meanPost_low_actv(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
    timePts = (X_low == n & (es.outcome(t_low) ~= 2));
    meanPost_low_pass(:,n) = nanmean(Posterior_low(timePts, :),1);
    %
    timePts = (X_high == n & (es.outcome(t_high) ~= 2));
    meanPost_high_pass(:,n) = nanmean(Posterior_high(timePts, :),1);
    %
    timePts = (X_high == n & (es.outcome(t_high) == 2));
    meanPost_high_actv(:,n) = nanmean(Posterior_high(timePts, :),1);
end

Posterior_all.meanPost.norm = meanPost_norm;
Posterior_all.meanPost.high = meanPost_high;
Posterior_all.meanPost.low =  meanPost_low;

% Posterior_all.meanPost.high_orig = meanPost_high_orig;
% Posterior_all.meanPost.low_orig =  meanPost_low_orig;

%% Getting phase vs. stuff
Posterior_all = phaseVsStuff(Posterior_all, 65, 1, 25);
end