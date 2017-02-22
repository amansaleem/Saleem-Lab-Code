classdef spikeMap < handle
    %     And abstract class that can get maps from neural data
    %       Methods include: train, plotonTraj, plotMap, testNew, ...
    %       Current version: 1D and 2D ready
    % Aman Saleem
    % Jan 2014
    
    properties
        dimensionality;       % 1D / 2D / 3D
        variable;   % the description of the variable decoder
        variable_range; % min and max values of the variable used while training 
        smth_win;
        model;
        performance; % cross-validated (if kfold>1) performance at training
        meanPerformance; % mean performance
        sampleRate;
        fullModel;
    end
    
    methods (Abstract)
        [obj]       = trainSpikeMap(obj, X, Y);
        [X]         = testMap(obj, Y);
        plotOnTraj(obj);
        plotMap(obj);
    end
    
    methods
        function obj = spikeMap(varargin)
            
            pnames = {'dimensionality' 'variable' 'variable_range'...
                'smth_win' 'model' 'sampleRate'};
            dflts  = {'1' 'P' []...
                150 [] 60};
            [obj.dimensionality, obj.variable, obj.variable_range, ...
                obj.smth_win, obj.model obj.sampleRate] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
        end
    end
end