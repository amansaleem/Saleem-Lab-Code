function [Pred] = runModelsDecPred(...
    area, istheta, spike_smth_win, ...
    box, numSampleBins, smth_win, trainCorrect, trainAllCorrect)
display(['Running area: ' area ' Theta: ' num2str(istheta) ...
    ' Spike smoothing: ' num2str(spike_smth_win) ...
    ' smoothing: ' num2str(smth_win) ' Train correct: ' num2str(trainCorrect)]) 
idx = 1;
for iseries = [1025, 1030, 530, 531, 601, 602, 603, 604]
    es = loadCA1V1data(iseries, istheta, area, spike_smth_win, box);
    %% Remove interneurons from the CA1 recordings
    if strcmp(area,'CA1')
        switch iseries
            case 1025
                % display('Warning!!! Removing interneurons on 1025');
                interneurons = zeros(1,size(es.spikeTrain,2));
                interneurons([17, 18,  9,11,12, 44, 45, 41, 78, 102, 104, 129, 113, 140, 145]) = 1;
                es.spikeIDs = es.spikeIDs(~interneurons);
                es.spikeTrain = es.spikeTrain(:,~interneurons);
            case 1030
                % display('Warning!!! Removing interneurons on 1030');
                interneurons = zeros(1,size(es.spikeTrain,2));
                interneurons([15, 20, 30, 39, 71, 132]) = 1;
                es.spikeIDs = es.spikeIDs(~interneurons);
                es.spikeTrain = es.spikeTrain(:,~interneurons);
            case 602
                % display('Warning!!! Removing interneurons on 602');
                interneurons = zeros(1,length(es.spikeIDs));
                interneurons([5 8 10 18 24 29 35 38 40 52 54 59 60 61 77 90 92 93 107 109 127]) = 1;
                es.spikeIDs = es.spikeIDs(~interneurons);
                es.spikeTrain = es.spikeTrain(:,~interneurons);
            case 604
                % display('Warning!!! Removing interneurons on 604');
                interneurons = zeros(1,length(es.spikeIDs));
                interneurons([8 22 25 26 19 31 40 52 53 66 91 97 109]) = 1;
                es.spikeIDs = es.spikeIDs(~interneurons);
                es.spikeTrain = es.spikeTrain(:,~interneurons);
            case 531
                % display('Warning!!! Removing interneurons on 531');
                interneurons = zeros(1,length(es.spikeIDs));
                interneurons([11 15 16 44 53 86 94 98]) = 1;
                es.spikeIDs = es.spikeIDs(~interneurons);
                es.spikeTrain = es.spikeTrain(:,~interneurons);
            case 603
                % display('Warning!!! Removing interneurons on 603');
                interneurons = zeros(1,length(es.spikeIDs));
                interneurons([3 8 9 15 57 63 64 66 74 116 131 138 171 177 182 185 188]) = 1;
                es.spikeIDs = es.spikeIDs(~interneurons);
                es.spikeTrain = es.spikeTrain(:,~interneurons);
            case 601
                % display('Warning!!! Removing interneurons on 601');
                interneurons = zeros(1,length(es.spikeIDs));
                interneurons([4 9 10 11 12 20 24 29 47 50 52 60 64 71 84]) = 1;
                es.spikeIDs = es.spikeIDs(~interneurons);
                es.spikeTrain = es.spikeTrain(:,~interneurons);
            case 530
                % display('Warning!!! Removing interneurons on 530');
                interneurons = zeros(1,length(es.spikeIDs));
                interneurons([4 8 36 39 56 62 72 78 84 95 98 101 99 102 103 105 131 155 163 167 207 217 222 229 231 266 270]) = 1;
                es.spikeIDs = es.spikeIDs(~interneurons);
                es.spikeTrain = es.spikeTrain(:,~interneurons);
        end
    end
    %% Get the conditions
    base = es.traj>1 & es.contrast~=0 ...
        & ~isnan(es.traj) & round(100*es.gain)/100==1 ...
        & round(100*es.roomLength)/100==1 & es.smthBallSpd>5 ...
        & es.trajspeed>=0;
    t_med  = round(es.contrast*100)/100==0.6 & base;
    t_low  = round(es.contrast*100)/100==0.18 & base;
    t_high = round(es.contrast*100)/100==0.72 & base;
    
    % Setting the firing rate in the inter trial interval to nan
    es.spikeTrain(es.traj==0,:) = nan;
    %% Get outcomes
    [es] = redefineOutcome(es);
    %% Initialise decoders (low, med, high)
    dec_m = bayesDecoder;
    dec_m.numBins = numSampleBins;
    dec_m.fixedSmth = smth_win;
    
    dec_l = bayesDecoder;
    dec_l.numBins = numSampleBins;
    dec_l.fixedSmth = smth_win;
    
    dec_h = bayesDecoder;
    dec_h.numBins = numSampleBins;
    dec_h.fixedSmth = smth_win;
    %% Train and run the decoders within condition
    spkRate = es.spikeTrain;
    if trainCorrect
        % temp to get the X
        dec = bayesDecoder;
        dec.numBins = numSampleBins; dec.fixedSmth = smth_win;
        [~, ~, X_med, ~] = dec.trainDecoder(es.traj(t_med), spkRate(t_med,:), 0);
        clear dec
        [dec_m, ~, ~, ~] = dec_m.trainDecoder( ...
            es.traj(t_med & es.outcome==2), spkRate(t_med & es.outcome==2,:), 0);
        [~,  Posterior_med] = dec_m.predictBayesDecoder(spkRate(t_med,:), 0,'mean');
        if trainAllCorrect
            % Low C
            dec = bayesDecoder;
            dec.numBins = numSampleBins; dec.fixedSmth = smth_win;
            [~, ~, X_low, ~] = dec.trainDecoder(es.traj(t_low), spkRate(t_low,:), 0);
            clear dec
            [dec_l, ~, ~, ~] = dec_l.trainDecoder( ...
                es.traj(t_low & es.outcome==2), spkRate(t_low & es.outcome==2,:), 0);
            [~,  Posterior_low] = dec_l.predictBayesDecoder(spkRate(t_low,:), 0,'mean');
            % High C
            dec = bayesDecoder;
            dec.numBins = numSampleBins; dec.fixedSmth = smth_win;
            [~, ~, X_high, ~] = dec.trainDecoder(es.traj(t_high), spkRate(t_high,:), 0);
            clear dec
            [dec_h, ~, ~, ~] = dec_h.trainDecoder( ...
                es.traj(t_high & es.outcome==2), spkRate(t_high & es.outcome==2,:), 0);
            [~,  Posterior_high] = dec_h.predictBayesDecoder(spkRate(t_high,:), 0,'mean');
        else
            [dec_l, ~, X_low, Posterior_low] = ...
                dec_l.trainDecoder( es.traj(t_low), spkRate(t_low,:), 0);
            [dec_h, ~, X_high, Posterior_high] = ...
                dec_h.trainDecoder( es.traj(t_high), spkRate(t_high,:), 0);
        end
    else
        [dec_m, ~, X_med, Posterior_med] = ...
            dec_m.trainDecoder( es.traj(t_med), spkRate(t_med,:), 0);
        [dec_l, ~, X_low, Posterior_low] = ...
            dec_l.trainDecoder( es.traj(t_low), spkRate(t_low,:), 0);
        [dec_h, ~, X_high, Posterior_high] = ...
            dec_h.trainDecoder( es.traj(t_high), spkRate(t_high,:), 0);
    end
    %% Test the decoders across condition
    % w.r.t. med
    [~,  Posterior_high_tmed] = dec_m.predictBayesDecoder(spkRate(t_high,:), 0,'mean');
    [~,  Posterior_low_tmed]  = dec_m.predictBayesDecoder(spkRate(t_low,:), 0,'mean');
    % w.r.t. high
    [~,  Posterior_low_thigh]  = dec_h.predictBayesDecoder(spkRate(t_low,:), 0,'mean');
    % w.r.t. low
    [~,  Posterior_high_tlow]  = dec_l.predictBayesDecoder(spkRate(t_high,:), 0,'mean');
    %% Getting the ML estimates
    % Within
    [~,X_ML_low]  = max(Posterior_low');
    [~,X_ML_med ] = max(Posterior_med');
    [~,X_ML_high] = max(Posterior_high');
    % Across
    [~,X_ML_high_m] = max(Posterior_high_tmed');
    [~,X_ML_low_m ] = max(Posterior_low_tmed');
    [~,X_ML_low_h ] = max(Posterior_low_thigh');
    [~,X_ML_high_l] = max(Posterior_high_tlow');
    %% Format output
    % Inputs
    Pred(idx).Inputs.area = area;
    Pred(idx).Inputs.numSampleBins = numSampleBins;
    Pred(idx).Inputs.smth_win = smth_win;
    Pred(idx).Inputs.trainAllCorrect = trainAllCorrect;
    % Original data
    Pred(idx).es = es;
    Pred(idx).t_med = t_med;
    Pred(idx).t_low = t_low;
    Pred(idx).t_high = t_high;
    % Decoders
    Pred(idx).decoders.med = dec_m;
    Pred(idx).decoders.low = dec_l;
    Pred(idx).decoders.high = dec_h;
    % Posteriors
    Pred(idx).posteriors.med = Posterior_med;
    Pred(idx).posteriors.low = Posterior_low;
    Pred(idx).posteriors.high = Posterior_high;
    Pred(idx).posteriors.lowTmed = Posterior_low_tmed;
    Pred(idx).posteriors.highTmed = Posterior_high_tmed;
    Pred(idx).posteriors.lowThigh = Posterior_low_thigh;
    Pred(idx).posteriors.highTlow = Posterior_high_tlow;
    % X
    Pred(idx).X.med = X_med;
    Pred(idx).X.low = X_low;
    Pred(idx).X.high= X_high;
    % ML
    Pred(idx).ML.med = X_ML_med;
    Pred(idx).ML.low = X_ML_low;
    Pred(idx).ML.high = X_ML_high;
    Pred(idx).ML.lowTmed = X_ML_low_m;
    Pred(idx).ML.highTmed = X_ML_high_m;
    Pred(idx).ML.lowThigh = X_ML_low_h;
    Pred(idx).ML.highTlow = X_ML_high_l;
    idx = idx + 1;
end
%% Save the data
[~,k] = system('hostname');
if strcmp('Aman-PC',k(1:end-1))
    DataDir = 'C:\Users\Aman\Dropbox\Work\Data\CA1V1\modelsAndDecoding';
else
    DataDir = 'E:\Dropbox\Work\Data\CA1V1\modelsAndDecoding';
end
% DataDir = 'C:\Users\aman\Dropbox\Work\Data\CA1V1\modelsAndDecoding';
if smth_win==5
    DataDir = [DataDir filesep 'Fixed5'];
elseif smth_win==10
    DataDir = [DataDir filesep 'Fixed10'];
elseif smth_win==100
    DataDir = [DataDir filesep 'Fixed100'];
elseif smth_win==8
    DataDir = [DataDir filesep 'Fixed8'];
elseif isempty(smth_win)
    DataDir = [DataDir filesep 'Optimal'];
end
if istheta
    DataDir = [DataDir filesep 'Theta'];
elseif spike_smth_win==250 & box==1
    DataDir = [DataDir filesep '250'];
elseif spike_smth_win==250 & box==0
    DataDir = [DataDir filesep '250gauss'];
end
if trainCorrect
    DataDir = [DataDir filesep 'Correct'];
else
    DataDir = [DataDir filesep 'CV'];
end
fileName = ['modelsAndPreds_' area '_Theta_' num2str(istheta) ...
    '_Samples_' num2str(numSampleBins) '_trainAllCorrect_' num2str(trainAllCorrect)];
save([DataDir filesep fileName],'Pred','-v7.3');    