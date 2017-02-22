classdef Decoder
    %     And abstract class that can decode neural data
    %       Methods include: train, predict (decode), ...
    %       Current decoders include: linear, treebagger and bayes
    % Aman Saleem
    % Oct 2013
    
    properties
        type;       % linear / treebagger / bayes
        variable;   % the description of the variable decoder
        variable_range; % min and max values of the variable used while training 
        smth_win = 150;
        model;
        performance; % cross-validated (if kfold>1) performance at training
        meanPerformance; % mean performance
        sampleRate = 60;
    end
    
    methods (Abstract)
        [obj, prediction] = trainDecoder(obj, X, Y, smth_win);
        [X] = predictVariable(obj, Y, smth_win);
    end
    methods
        function obj = Decoder(varargin)
            
            pnames = {'type' 'variable' 'variable_range'...
                'smth_win' 'model' 'sampleRate'};
            dflts  = {'linear' 'P' []...
                150 [] 60};
            [obj.type, obj.variable, obj.variable_range, ...
                obj.smth_win, obj.model obj.sampleRate] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
    end
end