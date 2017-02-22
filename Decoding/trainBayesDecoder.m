function [obj, prediction, V, fPosterior, fnPosterior, fnonNormPosterior] = trainBayesDecoder(obj, V, Y, smth_win, subset)
% program to implement a 1D linear decoder (using ridge regression) given training, cv and test
% sets , delay
% Usage: [Performance, Prediction, Model] = linearDecoder(Y, V, CVO, train_mean)
% Inputs: Y - firing rate
%         V - Variable being coded/decoded
%         CVO(optional) - crossvalidation object
%         train_mean(optional) - mean over the training set, special cases only
%
% Outputs:Performance - mean fraction of explained variance predicted by linear ridge decoder
%         Prediction  - prediction of the variable values by the decoder
%         Model       - The underlying model parameters (Useful for decoding novel stimuli)
%             WV - weight array (actual model)
%             Performance - same as performance above, but for individual iterations
%             train_mean  - needed for decoding (especially to calculate EV of decoded quantity)
%
% Aman Saleem
% Oct 2013

if isempty(obj.kfold)
    obj.kfold = 5;
end

if sum(isnan(V))>0 | sum(isnan(Y(:)))>0
    display('WARNING!!! Nans in the data: making a temp fix');
    t = ones(size(V));
    t(isnan(V)) = 0;
    t(isnan(sum(Y,2))) = 0;
    V = V(t>0);
    Y = Y(t>0,:);
end

if isempty(obj.CVO)
    obj.CVO = crossValPartition(ones(1,length(V)),obj.kfold);
end

if obj.kfold == 1
    obj.CVO.kfold = 1;
    obj.CVO.train{1}= ones(1,length(V));
    obj.CVO.cv{1}   = obj.CVO.train{1};
    obj.CVO.test{1} = obj.CVO.train{1};
end

if isempty(obj.train_mean)
    calTrainMean = 1;
else
    calTrainMean = 0;
end

for icell = 1:size(Y, 2);
    Y(:,icell) = smthInTime(Y(:,icell), obj.sampleRate, smth_win, 'same',subset);
end

obj.performance = zeros(1,obj.CVO.kfold);
Prediction = [];
% Bayes decoder needs discretization
[V obj.bins] = normalise1var(V, obj.numBins);
fPosterior = [];
fnPosterior = [];
fnonNormPosterior = [];

% Get the 1D maps
allMap = oneDimMap;
allMap.CVO = obj.CVO;
allMap.kfold = obj.kfold;
allMap.bins  = obj.bins; 
allMap.numBins = obj.numBins;
allMap = allMap.trainSpikeMap(V, Y, smth_win, obj.fixedSmth);

for iter = 1:obj.CVO.kfold
    
    Ytrain  = Y(obj.CVO.train{iter},:);
    Ycv     = Y(obj.CVO.cv{iter},:);
    Ytest   = Y(obj.CVO.test{iter},:);
    
    Vtrain  = V(obj.CVO.train{iter},:);
    Vcv     = V(obj.CVO.cv{iter},:);
    Vtest   = V(obj.CVO.test{iter},:);
    
    if calTrainMean
        train_mean = mean(Vtrain,1);
    else
        train_mean = obj.train_mean;
    end
    
    %% the main section of bayes decoder
    % Getting the 1D map for each neuron (calculating the place fields)
    % ...the main training component
    respModel = zeros(size(Y,2),obj.numBins);
    for icell = 1:size(Y,2)
%         [model(icell)] = get1Dmap(Y(:,icell), V', obj.numBins, obj.bins, obj.CVO, iter, obj.sampleRate, obj.smth_win);
        obj.model.trained(iter).respModel_orig(icell,:) = allMap.model.tuning(icell).respModel(iter,:);
        %         obj.model.trained(iter).respModel(icell,:) = model(icell).tuning./sum(model(icell).tuning);
        obj.model.trained(iter).respModel(icell,:) = allMap.model.tuning(icell).respModel(iter,:);
        obj.model.EV(iter,icell) = allMap.model.EV(iter,icell);
%         obj.model.L(iter,icell) = allMap.model.L(iter,icell);
%         obj.model.Q(iter,icell) = allMap.model.Q(iter,icell);
        tempModel(iter,icell,:) = obj.model.trained(iter).respModel(icell,:);
        minRate(iter,icell) = min(obj.model.trained(iter).respModel_orig(icell,:));
        maxRate(iter,icell) = max(obj.model.trained(iter).respModel_orig(icell,:));
    end
    obj.model.EV(iter,obj.model.EV(iter,:)<0) = 0;
%     goodCells = find(obj.model.EV(iter,:)>=0);
%     goodCells = find(obj.model.EV(iter,:)>=0 & (maxRate(iter,:)>2*minRate(iter,:)));
    goodCells = (max(obj.model.trained(iter).respModel'*obj.sampleRate)>=0);
    % Calculate the performance of the response model
    [Posterior, nonNormPosterior] = calcPosterior(obj.model.trained(iter).respModel(goodCells,:), Ytest(:,goodCells)); % Get the posterior estimates
    nPosterior = Posterior./repmat(median(Posterior),size(Posterior,1),1);
    nPosterior = nPosterior./(sum(nPosterior,2)*ones(1,size(nPosterior,2)));
    
    % Converting to log likelihood, which is independent of number of bins
    baseline    = 1./obj.numBins;
    Posterior   = log2(Posterior/baseline);
    nPosterior   = log2(nPosterior/baseline);
    
    [maxInTime] = max(Posterior,[],2);
    [minInTime] = min(Posterior,[],2);
    
    % Relative performance
    for t = 1:length(Vtest)
        probPeak(t,1) = Posterior(t,Vtest(t));
    end
%     obj.relPerformance(iter) = nanmean((maxInTime-probPeak)./(maxInTime-minInTime));
    obj.relPerformance(iter) = nanmean(probPeak);
    obj.confidence(iter)     = nanmean(maxInTime);
        
    [~, VFit] = max(Posterior,[],2);
    VFit = VFit';
    
    % Normal performance
    Perf(iter) = 1 - ((nansum((VFit' - Vtest).^2)))./((nansum((Vtest - train_mean).^2)));
    
    %
    obj.model.train_mean(iter) = train_mean;
    %%
    fPosterior = [fPosterior' Posterior']';
    fnPosterior = [fnPosterior' nPosterior']';
    fnonNormPosterior = [fnonNormPosterior' nonNormPosterior']';
    
    Prediction = [Prediction VFit];
end
Perf(Perf<0) = 0;
obj.performance = Perf;
obj.meanPerformance = nanmean(Perf);

% obj.model.meanModel = squeeze(nanmedian(tempModel,1)); % problem if only
% one cell;
obj.model.meanModel = (nanmean(tempModel,1));
obj.model.meanModel = reshape(obj.model.meanModel,size(obj.model.meanModel,2),size(obj.model.meanModel,3));

obj.model.bestModel = zeros(size(obj.model.meanModel));
[~,ibestPerf] = max(obj.relPerformance);
% [~,ibestPerf] = max(Perf);

for icell = 1:size(Y,2)
        obj.model.bestModel(icell,:) = obj.model.trained(ibestPerf).respModel(icell,:);
end

prediction  = Prediction';