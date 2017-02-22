classdef treeDecoder < Decoder
    %     Create a Decoder class object of type 'tree'. The decoder can be trained using the
    %     function 'trainDecoder' and used to decode with the function
    %     'predictVariable'
    %     Usage:
    % obj = treeDecoder;
    % optional arguments
    %     'variable': variable being coded/decoded
    %     'variable_range': min and max values of variable during training
    %     'smth_win': smth_win used for training
    %     'kfold': number of cross-validation windows, 1 fors no cross-validation
    %     'CVO': cross-validation structure
    %     'sampleRate': the sampling rate of the data.
    %
    % [obj, prediction, X] = trainDecoder(obj, X, Y, smth_win);
    % inputs: obj: object
    %         X: variable encoded/decoded
    %         Y: firing rate / signal
    %         smth_win(optional): the smoothing window to use on the signal (dfl: 150ms)
    % Ouputs: obj
    %         cross-validated prediction
    %         X (if there are NaN's in the data, these time points are taken out
    % [X] = predictVariable(obj, Y, smth_win);
    % inputs: obj: object
    %         Y: firing rate / signal
    %         smth_win(optional): the smoothing window to use on the signal (dfl: 150ms)
    % output: X: variable encoded/decoded
    %
    % Aman Saleem
    % October 2013
    
    properties
        kfold;       % number of partitions made to train the decoder
        train_mean;  % mean response used for training
        CVO;    % cross-validation structure
        maxNumTrees; % max number of tree for decoding, (dfl: 5)
    end
    
    methods
        function obj = treeDecoder(varargin)
            
            pnames = {'type' 'variable' 'variable_range'...
                'smth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'maxNumTrees'};
            dflts  = {'tree' 'P' []...
                150 []...
                5 []...
                [] [] 60 5};
            [obj.type, obj.variable, obj.variable_range, ...
                obj.smth_win, obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate obj.maxNumTrees] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
        function [obj, prediction, X] = trainDecoder(obj, X, Y, smth_win);
            % [obj, prediction, X] = trainDecoder(obj, X, Y, smth_win);
            % inputs: obj: object
            %         X: variable encoded/decoded
            %         Y: firing rate / signal
            %         smth_win(optional): the smoothing window to use on the signal (dfl: 150ms)
            % Ouputs: obj
            %         cross-validated prediction
            %         X (if there are NaN's in the data, these time points are taken out
            if nargin<4
                smth_win = obj.smth_win;
            end
            [obj, prediction, X] = trainTreeDecoder(obj, X, Y, smth_win);
        end
        
        function [X] = predictVariable(obj, Y, smth_win);
            % [X] = predictVariable(obj, Y, smth_win);
            % inputs: obj: object
            %         Y: firing rate / signal
            %         smth_win(optional): the smoothing window to use on the signal (dfl: 150ms)
            % output: X: variable encoded/decoded
            if nargin<3
                smth_win = obj.smth_win;
            end
            
            t = ones(size(Y,1),1);
            t(isnan(sum(Y,2))) = 0;
            t = t>0;
            X = NaN*ones(size(Y,1),1);
            
            [X(t)] = predictTreeDecoder(obj, Y(t,:), smth_win);
        end
        
        function [prediction] = predictLinearDecoder(obj, Y, smth_win);
            for icell = 1:size(Y, 2);
                Y(:,icell) = smthInTime(Y(:,icell), obj.sampleRate, smth_win);
            end
            prediction = predict(obj.model.meanModel, Y, 'trees', [1:obj.model.numTrees]);
        end
    end
end

