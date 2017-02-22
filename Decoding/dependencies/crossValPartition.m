function [CVO] = crossValPartition(input, kfold, sepCV)

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
if nargin<2
    kfold = 5;
end
if nargin<3
    sepCV = 1;
end

CVO.input = input;
CVO.kfold = kfold;
CVO.separateCV = sepCV;

% starting the iterations of cross-validation
kidx = floor(length(in)./kfold).*([0:kfold]);

% CVO = cvpartition(ones(size(normR)),'k',kfold);
for iter = 1:kfold
    % setting up the training, cross-validation and test sets
    train = false(size(in));
    cv = train;
    test = cv;
        
    test(kidx(iter)+1:kidx(iter+1))         = true;
    if sepCV
        if iter == kfold
            cv(kidx(1)+1:kidx(2))           = true;
        else
            cv(kidx(iter+1)+1:kidx(iter+2)) = true;
        end
    else
        cv = test;
    end
    train(~cv & ~test)   = true;
    
    CVO.train{iter}     = false(size(input));
    CVO.test{iter}      = false(size(input));
    CVO.cv{iter}        = false(size(input));

    CVO.train{iter}(in(train)) = true;
    CVO.test{iter}(in(test)) = true;
    CVO.cv{iter}(in(cv)) = true;
end