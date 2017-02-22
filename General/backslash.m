function backslash(ianimal, iseries)
% close all;
% clear;
% clc;

addpath('\\ZSERVER\Lab\Share\Saleem\V1 gamma\Data');

load([num2str(ianimal) '_' num2str(iseries) '_' 'freq_V1_CA1']);
es = eval(['es_' num2str(ianimal) '_' num2str(iseries)]);
eval(['clear es_' num2str(ianimal) '_' num2str(iseries)]);

% load 501_531_freq_V1_CA1
% es = es_501_531;
% clear es_501_531;

% load 502_603_freq_V1_CA1
% es = es_502_603;
% clear es_502_603;

inputs2use = [1:6];
% {1}= lickMatrix;   {2}= rewardMatrix;   {3}= vrStartMatrix;
% {4}= vrEndMatrix;  {5}= runStartMatrix; {6}= runEndMatrix;
addBiasTerm = 1;
freqMax = 75;
freqMin = 55;
nLags2omit = 2;
maxlag = 10;
order = 2*(maxlag+nLags2omit); % order of the kernels length = order+1
tol = 0;

t = es.contrast~=0;

grays = es.traj;
grays(es.traj==0) = 0;
grays(es.traj~=0) = 1;
diffGrays = ([0 diff(grays)']');

running = zeros(size(es.traj));
running(es.smthBallSpd>2) = 1;
diffRun = ([0 diff(running)']');

VR_starts   = find(diffGrays==1 & t);
VR_ends     = find(diffGrays==-1 & circshift(t,2));
lickTimes   = find(es.lick & t);
rewTimes    = find(~isnan(es.reward) & t);
run_starts  = find(diffRun==1 & (es.lick==0) & isnan(es.reward));
run_ends    = find(diffRun==-1 & (es.lick==0) & isnan(es.reward));

% building the matrices for the backslash operation later


lickMatrix = toeplitz(single(es.lick & t), zeros(1, order+1));
% lickMatrix = toeplitz(single(es.lick), zeros(1, order+1));
rewardMatrix = toeplitz(single(~isnan(es.reward)), zeros(1, order+1));
vrStartMatrix = toeplitz(single(diffGrays==1 & t), zeros(1, order+1));
vrEndMatrix = toeplitz(single(diffGrays==-1 & circshift(t,2)), zeros(1, order+1));
runStartMatrix = toeplitz(single(diffRun==1), zeros(1, order+1));
runEndMatrix = toeplitz(single(diffRun==-1), zeros(1, order+1));

inputs{1}= lickMatrix;
inputs{2}= rewardMatrix;
inputs{3}= vrStartMatrix;
inputs{4}= vrEndMatrix;
inputs{5}= runStartMatrix;
inputs{6}= runEndMatrix;

% inputs{1}= fliplr(lickMatrix);
% inputs{2}= fliplr(rewardMatrix);
% inputs{3}= fliplr(vrStartMatrix);
% inputs{4}= fliplr(vrEndMatrix);
% inputs{5}= fliplr(runStartMatrix);
% inputs{6}= fliplr(runEndMatrix);

labels{1} = 'lick';
labels{2} = 'reward';
labels{3} = 'vrStart';
labels{4} = 'vrEnd';
labels{5} = 'runStart';
labels{6} = 'runEnd';

% eventMatrix = [lickMatrix, rewardMatrix, vrStartMatrix, vrEndMatrix, runStartMatrix, runEndMatrix];
eventMatrix = cell2mat(inputs(inputs2use));
% eventMatrix = fliplr(eventMatrix);
% adding a column of ones to account for dc bias
if addBiasTerm
    eventMatrix = [eventMatrix, ones(size(eventMatrix, 1), 1)];
end

% imagesc(eventMatrix);

freqIdx = es.freq>freqMin & es.freq<freqMax;
freqAxis = es.freq(freqIdx);
spMatrix = es.powB(:, freqIdx);%.*repmat(freqAxis, size(es.powA, 1), 1);
% imagesc(1:size(spMatrix, 1), freqAxis, log(spMatrix)');
% set(gca, 'YDir', 'normal')

% hGauss = fspecial('gauss', [1 11], 3);
% spMatrix = imfilter(spMatrix, hGauss);
% hGauss = fspecial('gauss', 3);
% spMatrix = imfilter(spMatrix, hGauss);

%% backslash

% kernels = pinv(eventMatrix)*log(spMatrix);
if tol>0
    kernels = pinv(eventMatrix,tol)*(spMatrix);
else
    kernels = pinv(eventMatrix)*(spMatrix);
end
%
%% plotting
figure

climVals = [min(kernels(:)), max(kernels(:))];
if addBiasTerm
    outputs = mat2cell(kernels, [ones(length(inputs2use), 1)*(order+1); 1], size(kernels, 2));
else
    outputs = mat2cell(kernels, [ones(length(inputs2use), 1)*(order+1)], size(kernels, 2));
end

tauAxis = -order/2:order/2;

nOutputs = length(outputs);
nRows = floor(sqrt(nOutputs));
nColumns = ceil(nOutputs/nRows);
for iOutput = 1:nOutputs
    subplot(nRows, nColumns, iOutput)
    if size(outputs{iOutput}, 1)==1
        imagesc(0, freqAxis, outputs{iOutput}');
    else
        imagesc(tauAxis(nLags2omit+1:end-nLags2omit), freqAxis, outputs{iOutput}(nLags2omit+1:end-nLags2omit, :)');
    end
        set(gca, 'Clim', climVals);
        colorbar;
    set(gca, 'YDir', 'normal')
    try
        title(labels{inputs2use(iOutput)});
    end
end

%% forward model check

% we might consider using kernels without the omitted lags, but then we
% need to regenerate eventatrix and crop the kernels matrix smartly
spEst = eventMatrix*kernels;
climVals = [min(spMatrix(:)), max(spMatrix(:))];
figure
subplot(2, 1, 1)
imagesc(size(spEst,1), freqAxis,spMatrix');
axis xy
set(gca, 'Clim', climVals);
title('The original spectrogram');
subplot(2, 1, 2)
imagesc(size(spEst,1), freqAxis,spEst')
axis xy
set(gca, 'Clim', climVals);
title('The estimated spectrogram');