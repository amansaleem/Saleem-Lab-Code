function [autocorr_val, shiftedMeanModel, meanModel] = plotPlaceFields_new(model, t, sortOrNot)

if nargin<2
    t = true(size(nanmean(model.EV,1)));
end

if nargin<3
    sortOrNot = 1;
end

EV = nanmean(model.EV(:,t),1);
% EV = EV(t);
cellList = t;


cellIDs = 1:length(t);
cellIDs = cellIDs(t);

meanModel = model.meanModel;
meanModel = meanModel(cellList,:);

[~,maxPos] = max(meanModel');
if sortOrNot
    if sortOrNot==2
        [~,sort_order] = sort(EV);
    else
        [~,sort_order] = sort(maxPos);
    end
else
    sort_order = 1:length(maxPos);
end

    autocorr_val = [];
for icell = 1:sum(cellList)
%     [autocorr_v(icell,:)] = autocorr(meanModel(icell,:),size(meanModel,2)-1);
    meanModel(icell,:) = normalise1var(meanModel(icell,:));

%     meanModel(icell,:) = meanModel(icell,:)  ./  max(meanModel(icell,:));

%     meanModel(icell,:) = meanModel(icell,:)  -  mean(meanModel(icell,:));
%     meanModel(icell,:) = meanModel(icell,:) ./ std(meanModel(icell,:));
%     meanModel(icell,:) = meanModel(icell,:)  ./  (mean(meanModel(icell,:)).^2);
%     autocorr_val(icell,:) = find(autocorr_v(icell,:)<0,1,'first');
end

% autocorr_val = autocorr_val(sort_order);
% b = 0;
% b(1) = 0;
% autocorr_val = sum(autocorr_v(sort_order,:)>b(1),2);
% cross50 = sum(meanModel(sort_order,:)>0.5,2);
[~,idx] = max(meanModel');
shiftedMeanModel = nan*ones(size(meanModel,1),100);
shiftedMeanModel(:,51:100) = meanModel;
% [~,sort_idx] = sort(idx);
% meanModel = meanModel(sort_idx,:);
numCells = 0;
for n = 1:size(meanModel,1)
    if idx(n)<40
        shiftedMeanModel(n,:) = circshift(shiftedMeanModel(n,:),[0 -idx(n)]);
        numCells = numCells +1;
    end
end
% imagesc(shiftedMeanModel,[-3 3])
errorarea_as(-97:2:102, nanmean(shiftedMeanModel,1), nansem(shiftedMeanModel))
% imagesc(meanModel(sort_order,:),[-3 3]);
% imagesc(meanModel(sort_order,:),[0 1]);
axis xy;
% RedWhiteBlue;
set(gca, 'box','off','TickDir','out','fontsize',14,'color','none');
% set(gca, 'YTick',1:length(sort_order), 'YTickLabel', cellIDs(sort_order));
% set(gca, 'YTick',1:length(sort_order), 'YTickLabel', (100/250)*autocorr_val);
set(gca, 'YTick',1:length(sort_order), 'YTickLabel', round(1000*EV(sort_order))/1000);
title(num2str(numCells));