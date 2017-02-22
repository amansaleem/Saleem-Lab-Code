classdef kernel < handle
    
    %   Object based code to calculate Kernels
    %___________________________________________
    %   Usage:
    %   obj = kernel;
    %    optional arguments
    %     'smth_win': smoothing to be applied to the stimulus and response
    %     (def: 0 or no smoothing)
    %     'downSampleRate': the rate at which the kernel is based on (dfl:
    %     5 Hz)
    %     'sampleRate': original sample rate of data (dfl: 60 Hz)
    %     'tSamples': time of samples in time (seconds)
    %     'tFutureSamples': section of tSamples that is devoted to the future
    %     'CVO': cross-validation structure
    %     'zscore': whether or not to zscore the stimulus
    %     'ridgeValue': the value of lambda for ridge regression
    % %%%%%%%%
    % The functions to call:
    %___________________________________________________________
    % [obj, kernels, kStim, kResp] = obj.buildKernel(stim, resp);
    % inputs: stim: stimulus (nT x nStim)
    %         resp: firing rate / signal (nT x nCells)
    %
    % Ouputs:
    %         kernels: cross-validated prediction
    %         kStim: the processed (downsampled) stimulus
    %         kResp: the processed (downsampled) resp
    %________________________________________________
    % [kStim, kResp] = obj.processInputs(stim, resp);
    % inputs: same at previous, This is so that it can be run on a new
    % section / sample of the data, to implement cross-validation
    %
    %______________________________________
    % [kResp] = obj.predictResponses(kStim);
    % This is to predict the responses given the stimulus in the right
    % format. Again seperated so it can be used for cross-validation
    %%%%%%%%%%%%%
    % Example usage:
    %     stim = rand(50000,1);
    %     resp(:,1) = conv(stim, sin(0:(pi/50):50),'same');
    %     resp(:,2) = conv(stim, sin(0:(pi/100):100),'same');
    %     resp(:,3) = conv(stim, sin(0:(pi/200):200),'same');
    %
    %     k = kernel('sampleRate', 100);
    %     k.tSamples = 6;
    %     k.tFutureSamples = 3;
    %     [k,kernels, kstim, kresp] = k.buildKernel(stim, resp);
    %
    %     tresp = k.predictResponses(kstim);
    %
    %     figure;
    %     for n = 1:3
    %         subplot(3,4,(n-1)*4+1)
    %         plot(k.shiftValues, k.kernels(:,n))
    %
    %         subplot(3,4,(n-1)*4+2:n*4);
    %         plot(kresp(:,n),'k'); hold on; plot(tresp(:,n),'r'); hold off;
    %     end
    %%%%%%%%%%%%%%%%%%%
    % Aman Saleem
    % March 2014
    
    % Edits: 03/14: AS added rounding the rates, so resample doesn't crash.
    
    properties
        sampleRate;       % number of partitions made to train the decoder
        downSampleRate;  % mean response used for training
        CVO;    % cross-validation structure (unused for now_
        tSamples; % (in sec)
        tFutureSamples; % (in sec)
        tHistorySamples; % (in sec)
        smth_win; % smoothing window of input and responses (in ms)
        zscore; % whether or not to zscore the stim
        ridgeValue;
        kernels;
        nStim = 1;
        numShifts;
        numCells;
        shiftValues = [];
    end
    
    methods
        function obj = kernel(varargin)
            
            pnames = {'smth_win'...
                'sampleRate' 'downSampleRate'...
                'tSamples' 'tFutureSamples' 'CVO' ...
                'zscore' 'ridgeValue'};
            dflts  = {0 ...
                60 5 ...
                5 1 [] ...
                false 0.001 ...
                };
            [obj.smth_win, obj.sampleRate, obj.downSampleRate, ...
                obj.tSamples, obj.tFutureSamples,...
                obj.CVO, obj.zscore obj.ridgeValue...
                ] = ...
                internal.stats.parseArgs(pnames,dflts,varargin{:});
            obj.tHistorySamples = obj.tSamples - obj.tFutureSamples;
        end
        
        function [obj, kernels, kStim, kResp] = buildKernel(obj, stim, resp)
            obj.tHistorySamples = obj.tSamples - obj.tFutureSamples;
            if ~isinteger(obj.downSampleRate) | ~isinteger(obj.downSampleRate)
                display('WARNING! Sampling rate not integer, rounding it. Check the rates in the output object.')
            end
            obj.downSampleRate = round(10*obj.downSampleRate)/10;
            obj.sampleRate = round(10*obj.sampleRate)/10;
            [kStim, kResp] = obj.processInputs(stim, resp);
            [obj, kernels] = obj.trainKernel(kStim, kResp);
        end
        
        function [kStim, kResp] = processInputs(obj, stim, resp)
            %% Handle the responses
            smth_resp = zeros(size(resp));
            
            for icell = 1:size(resp,2)
                smth_resp(:,icell) = obj.smthInTime_local(resp(:,icell), obj.sampleRate, obj.smth_win);
                % 'dss' stands for down-sampled smooth
                dss_resp(:,icell) = resample(smth_resp(:,icell), round(10*obj.downSampleRate), round(10*obj.sampleRate));
            end
            
            kResp = dss_resp;
            
            %% Handle the stimulus
            smth_stim = zeros(size(stim));
            nStim = size(stim,2);
            obj.nStim = nStim;
            
            for istim = 1:nStim
                smth_stim(:,istim) = obj.smthInTime_local(stim(:,istim), obj.sampleRate, obj.smth_win);
                dss_stim(:,istim) = resample(smth_stim(:,istim), round(10*obj.downSampleRate), round(10*obj.sampleRate));
                if obj.zscore
                    dss_stim(:,istim) = (dss_stim(:,istim) - nanmean(dss_stim(:,istim))) ./ std(dss_stim(:,istim));
                end
            end
            
            ShiftDs = -floor(obj.tHistorySamples*obj.downSampleRate) : ceil(obj.tFutureSamples*obj.downSampleRate);
            numShifts = length(ShiftDs);
            obj.numShifts = numShifts;
            obj.shiftValues = ShiftDs*(1./obj.downSampleRate);
            
            kStim = zeros(size(dss_stim,1), numShifts*nStim);
            for iShift = 1:numShifts
                kStim(:, (iShift-1)*nStim+(1:nStim)) = circshift(dss_stim, ShiftDs(iShift));
            end
        end
        
        function [obj, kernels] = trainKernel(obj, X, Y)
            %            % Solving for Y = xA
            %            % Inversion is A = pinv(X'X+lambda) X'Y
            lambda = obj.ridgeValue;
            obj.kernels = pinv(X'*X - lambda*eye(size(X,2)))*(X'*Y);
            %             obj.kernels = (X + normrnd(0, lambda, size(X))) \ Y;
            
            ShiftDs = -floor(obj.tHistorySamples*obj.downSampleRate) : ceil(obj.tFutureSamples*obj.downSampleRate);
            numShifts = length(ShiftDs);
            numCells = size(Y,2);
            obj.numCells = numCells;
            
            kernels = zeros(obj.nStim, numShifts, numCells);
            for icell = 1:numCells
                for iShift = 1:numShifts
                    kernels(:,iShift,icell) = obj.kernels((iShift-1)*obj.nStim+(1:obj.nStim),icell);
                end
            end
        end
        
        function [resp] = predictResponses(obj, stim)
            
            resp = reshape((stim * obj.kernels),[],obj.numCells);
            
        end
        
        function out = smthInTime_local(obj, in, sampFreq, smthWin)
            %
            % Usage: out = smthInTime(in, sampFreq, smthWin, type subset)
            % to smooth in time the signal 'in'
            %
            % sampFreq = the sampling freq
            % smthWin  = the half width of the smoothing window (in ms), 0 gives no
            % smoothing
            % type     = Shape argument of conv (optional), 'same' is default
            % subset   = logical of region to be considered same size as input, everything else comes out
            % as NaNs
            %
            % Aman Saleem
            % November 2013
            
            
            subset = true(size(in));
            
            type = 'same';
            
            
            if smthWin==0
                out = in;
            else
                nanentries = isnan(in);
                subset = ~nanentries & subset;
                
                in(nanentries) = 0;
                samp_int = 1000/sampFreq;
                
                win = round(smthWin/samp_int);
                
                sGrid = max([100 win*5]);
                s = (-sGrid:sGrid);
                
                sfilt = (1./(win*(sqrt(2*pi)))).*exp(-(s.^2)./(2*(win^2)));
                
                % padding the inputs
                pad = zeros(size(sfilt));
                pad_length = length(pad);
                if size(in,1)~=1
                    flip_again = true;
                    in = in';
                    subset = subset';
                else
                    flip_again = false;
                end
                pad_in = [pad in pad];
                pad_subset = [pad double(subset) pad];
                
                pad_out         = conv(pad_in, sfilt, type);
                norm_subset     = conv(pad_subset, sfilt, type);
                
                pad_out = pad_out./norm_subset;
                out = pad_out((pad_length+1) : (end-pad_length));
                out(~subset) = NaN;
                
                if flip_again
                    out = out';
                end
            end
            
            
        end
    end
end
    
