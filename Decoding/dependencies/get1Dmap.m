function [model, spred] = get1Dmap(zua, variable, nGrid, bins, CVO, iter, sampFreq, win)
% function to fit and test a 1D model, find the optimal smoothing
% function too

X_train = variable(CVO.train{iter}>0);
X_test  = variable(CVO.test{iter}>0);
X_cv    = variable(CVO.cv{iter}>0);

% for icell = 1:numCells
if nargin<7
    sampFreq = 200;
    win = 10;
    smth = 0;
else 
    smth = 1;
end

if smth
    zua_train = smthInTime(zua(CVO.train{iter}>0), sampFreq, win);
    zua_test  = smthInTime(zua(CVO.test{iter}>0), sampFreq, win);
    zua_cv    = smthInTime(zua(CVO.cv{iter}>0), sampFreq, win);
else
    zua_train = zua(CVO.train{iter}>0);
    zua_test  = zua(CVO.test{iter}>0);
    zua_cv    = zua(CVO.cv{iter}>0);
end
test = zua(CVO.test{iter});

% get spike count map
scMap = full(sparse(X_train, 1, (zua_train), nGrid, 1));
% get occupancy map
occMap = full(sparse(X_train, 1, 1, nGrid, 1));

n1 = 10;
grids = 1:n1;

EV = zeros(1,length(grids));
% special_smooth_1d(input, win, bins, nGrid)
%get the FR map by smoothing
for win=1:length(grids)
    FRMap = special_smooth_1d(scMap, 1./grids(win), bins, n1)...
        ./special_smooth_1d(occMap, 1./grids(win), bins, n1);
    pred  = (FRMap(X_cv));
    if smth
        spred    = smthInTime(pred, sampFreq, win);
    else
        spred = pred;
    end
    EV(win) = calCrossValExpVar(zua_train, zua_cv, spred);
end
EV = round(100*EV);
[~, idx] = max(EV);
model.swin = grids(idx);
model.bins = bins;
model.tuning = special_smooth_1d(scMap, 1./model.swin, bins, n1)...
    ./special_smooth_1d(occMap, 1./model.swin, bins, n1);
pred  = (model.tuning(X_test));

if smth
    spred    = smthInTime(pred, sampFreq, win);
else
    spred = pred;
end

[model.EV model.corr model.L model.Q model.train_mean] = calCrossValExpVar(zua_train, zua_test, spred, test, pred);

end