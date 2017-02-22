classdef oneDimMap < spikeMap
    %     Create a spikeMap class object of type '1D'. The decoder can be trained using the
    %     function 'trainDecoder' and used to decode with the function
    %     'predictVariable'
    %     
    % Aman Saleem
    % Jan 2014
    properties
        kfold;       % number of partitions made to train the decoder
        train_mean;  % mean response used for training
        CVO;         % cross-validation structure
        bins;        % the values of the bin limits
        numBins;     % number of bins by which the variable is discretized
    end
    
    methods
        function obj = oneDimMap(varargin)
            
            pnames = {'dimensionality' 'variable' 'variable_range'...
                'smth_win' 'model'...
                'kfold' 'performance'...
                'train_mean' 'CVO' 'sampleRate' 'numBins'};
            dflts  = {'1D' 'P' []...
                150 []...
                5 []...
                [] [] 60 50};
            [obj.dimensionality, obj.variable, obj.variable_range, ...
                obj.smth_win, obj.model,...
                obj.kfold, obj.performance,...
                obj.train_mean, obj.CVO, obj.sampleRate obj.numBins] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
        
        function [obj, prediction, X] = trainSpikeMap(obj, X, Y, smth_win);
            if nargin<4
                smth_win = obj.smth_win;
            else
                obj.smth_win = smth_win;
            end
            [obj, prediction, X] = trainOneDimMap(obj, X, Y, smth_win);
        end
        
        function plotOnTraj(obj, X, Y);
            
        end
        
        function plotMap(obj, Y, smth_win);
            
        end
        
        function [X] = testMap(obj, Y);
        
        end
    end
end

