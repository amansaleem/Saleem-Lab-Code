function [obj, prediction, V] = trainLinearDecoder(obj, V, Y, smth_win)
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
% Based on the methods used in Saleem AB et al, 2013
% Oct 2013

if isempty(obj.kfold)
    obj.kfold = 5;
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

if sum(isnan(V))>0 | sum(isnan(Y(:)))>0
    display('WARNING!!! Nans in the data: making a temp fix');
    t = ones(size(V));
    t(isnan(V)) = 0;
    t(isnan(sum(Y,2))) = 0;
    V = V(t>0);
    Y = Y(t>0,:);
end
    
for icell = 1:size(Y, 2);
        Y(:,icell) = smthInTime(Y(:,icell), obj.sampleRate, smth_win);
end

obj.performance = zeros(1,obj.CVO.kfold);
Prediction = [];

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
    
    %% the main section of ridge regression
    lambda = ([-1000 -500 -200 -70 -50 -10 -5 -1 0 1 5 10 50 70 200 500 1000]);
    EV = zeros(1,length(lambda));
    for i = 1:length(lambda)
        WV = ((Vtrain)'*Ytrain*pinv(Ytrain'*Ytrain - lambda(i)*eye(size(Ytrain,2))));
        VFit = WV*Ycv';
        EV(i) = 1 - ((sum((VFit' - Vcv).^2)))./((sum((Vcv - mean(Vtrain)).^2)));
    end
    [~, i] = max(EV);
    
    WV = ((Vtrain)'*Ytrain*pinv(Ytrain'*Ytrain - lambda(i)*eye(size(Ytrain,2))));
    VFit = WV*Ytest';
    
    Perf(iter) = 1 - ((sum((VFit' - Vtest).^2)))./((sum((Vtest - train_mean).^2)));
    %%
    Prediction = [Prediction VFit];
    
    obj.model.WV{iter} = WV;
    overall(iter,:) = WV;
    obj.model.train_mean(iter) = train_mean;
end
Perf(Perf<0) = 0;
obj.performance = Perf;
obj.meanPerformance = nanmean(Perf);
obj.model.meanModel = nanmean(overall,1);

prediction  = Prediction';