function [Y, AUC,  meanCumSum, topCumSum, botCumSum] = booterror(distr, iter)

if nargin<2
    iter = 10000;
end

% sampleSize
ss = length(distr);

resample_idx = ceil(rand(iter,ss)*ss);

resample_distr = distr(resample_idx); %generates a matrix with c=no shuffles, r=no. iter

X = min([min(distr) 0]):.1:1;
% X=0:0.1:1;
allHist = hist(resample_distr', X);


Y = cumsum(allHist);
AUC = sum(Y)/length(distr);
% yTop = cumsum(allHist(10:11,:));
if nargout>2
    meanCumSum  =  mean(Y');
    topCumSum   =  meanCumSum + std(Y');
    botCumSum   =  meanCumSum - std(Y');
    meanCumSum  =  meanCumSum/length(distr);
    topCumSum   =  topCumSum/length(distr);
    botCumSum   =  botCumSum/length(distr);
end