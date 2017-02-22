function [Posterior, nonNormPosterior] = calcPosterior(respModel,Y)
% Function to calculate the posterior probability of a stimulus set (Posterior), given the
% population firing rate (Y) and response model (respModel)
% Adapted from code from Daniel Bendor (Bendor & Wilson, Nat Neuro, 2013)
% 
% Usage: [Posterior] = calcPosterior(respModel,Y)
% 
% Aman Saleem
% October 2013
T       = size(Y,1);        %no.of time bins
numCells= size(Y,2);        %no.of Cells
M       = size(respModel,2);%no.of response model bins
if min(respModel(:)) < (max(respModel(:))/(100*M))
    respModel = respModel + (max(respModel(:))/(100*M));
end
% respModel(sum(isnan(respModel'))>1,:)=0;
ProdFields = log(ones(T,M));
sumFields = zeros(1,M);
for icell = 1:numCells
    if sum(isnan(respModel(icell,:)))>1 % | max(respModel(icell,:)*60)<2
        continue
    end
    spkCount = Y(:,icell);
    ProdFields = ProdFields + log((ones(T,1)*respModel(icell,:)).^(spkCount*ones(1,M)));
    sumFields = sumFields + respModel(icell,:);
end
% sumFields  = nansum(respModel,1);

nonNormPosterior   = ProdFields + log((ones(T,1)*exp(-sumFields)));
nonNormPosterior   = exp(nonNormPosterior);
Posterior  = nonNormPosterior./(sum(nonNormPosterior,2)*ones(1,M));
end