function [obj, prediction, V] = trainTreeDecoder(obj, V, Y, smth_win)
% program to implement a tree bagger regression based decoder given
% training, cv (used to find optimal tree size) and test sets
% Usage: [Performance, Prediction, Model] = treeDecoder(Y, V, CVO, train_mean)
% Inputs: Y - firing rate
%         V - Variable being coded/decoded
%         CVO(optional) - crossvalidation object
%         maxNumTrees(optional) - maximum no.of trees to use in the
%         treebagger prog, less trees for faster implementation, more for
%         performance. Overfitting not much an issue because find the best
%         no.of trees through cross-validation.
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

if isempty(obj.maxNumTrees)
    obj.maxNumTrees = 4;
end

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
    
    %% the main section of tree bagger decoding
    b = TreeBagger(obj.maxNumTrees,[Ytrain],Vtrain,'oobpred','on','OOBVarImp', 'on', 'Method','regression');
    display(num2str(iter))
    parfor numTrees = 1:obj.maxNumTrees
        VFit      = predict(b,Ycv,'trees',[1:numTrees]);
        EV(numTrees) = 1 - ((sum((VFit - Vcv).^2)))./((sum((Vcv - mean(Vtrain)).^2)));
    end
    [~, numTrees] = max(EV);
    VFit      = predict(b,Ytest,'trees',[1:numTrees]);
    
    Perf(iter) = 1 - ((sum((VFit - Vtest).^2)))./((sum((Vtest - mean(Vtrain)).^2)));
    
    %%
    Prediction = [Prediction' VFit']';
    
    obj.model.tree{iter}.model = b;
    obj.model.tree{iter}.numTrees = numTrees;
    obj.model.train_mean(iter) = train_mean;
    nt(iter) = numTrees;
end
Perf(Perf<0) = 0;
obj.meanPerformance = nanmean(Perf);
obj.performance = Perf;
obj.model.numTrees = nanmean(nt);
obj.model.meanModel = TreeBagger(obj.maxNumTrees,[Y],V,'oobpred','on','OOBVarImp', 'on', 'Method','regression');

prediction  = Prediction';