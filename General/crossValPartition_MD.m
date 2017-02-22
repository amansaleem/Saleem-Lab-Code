function [CVO] = crossValPartition_MD(CVO, input, trialInd, kfold, sepCV )

%%modified by Mika to account for shuffling of trials

% Usage: [CVO] = crossValPartition(input, kfold, sepCV);
% 
% This is a function that creates a kfold partition of the array 'input'
% input - binary array. Only the 1s are used to create the crossvalidation set
% kfold - no.of divisions
% sepCV - create a separate cross-validation set, independent of test set
% 
% CVO - output stucture
% CVO.input = input;
% CVO.kfold = kfold;
% CVO.separateCV = true/false;
% CVO.train{n};
% CVO.test{n};
% CVO.cv{n};

in = find(input);
if nargin<4
    kfold = 5;
end
if nargin<5
    sepCV = 1;
end

CVO.input = input;
CVO.kfold = kfold;
CVO.separateCV = sepCV;

%do this to make sure that both 2d and 1d model
%use the same crossval samples
if ~isfield(CVO, 'trialPerms')  
%shuffle trial indices
perms = randperm(trialInd(end));
else
perms = CVO.trialPerms;    
end
% starting the iterations of cross-validation
kidx = floor(length(perms)./kfold).*([0:kfold]);
kidx(end) = length(perms);



for iter = 1:kfold
    % setting up the training, cross-validation and test sets
    train = false(size(in));
    cv = train;
    test = cv;
    
    iperm = perms(kidx(iter)+1:kidx(iter+1));
    for ip = iperm
    test(trialInd == ip) = true;
    end
    %test(perms(kidx(iter)+1:kidx(iter+1)))         = true;
    
    if sepCV
        if iter == kfold
            iperm = perms(kidx(1)+1:kidx(2));
            for ip = iperm
            cv(trialInd == ip) = true;
            end
            %cv(kidx(1)+1:kidx(2))           = true;
        else
            iperm = perms(kidx(iter+1)+1:kidx(iter+2));
            for ip = iperm
            cv(trialInd == ip) = true;
            end
            %cv(kidx(iter+1)+1:kidx(iter+2)) = true;
        end
    else
        cv = test;
    end
    train(~cv & ~test)   = true;
    
    
     CVO.train{iter}     = train;
    CVO.test{iter}      = test;
    CVO.cv{iter}        = cv;
    
%     CVO.train{iter}     = false(size(input));
%     CVO.test{iter}      = false(size(input));
%     CVO.cv{iter}        = false(size(input));
% 
%     CVO.train{iter}(in(train)) = true;
%     CVO.test{iter}(in(test)) = true;
%     CVO.cv{iter}(in(cv)) = true;
end
CVO.trialPerms = perms;