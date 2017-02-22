% [Pred] = runModelsDecPred_decodedPos(...
%     area, istheta, spike_smth_win, ...
%     box, numSampleBins, smth_win, ...
%     trainCorrect, trainAllCorrect)
for numBins = [200]
    %% Theta,
    for istheta = 1:2
        spike_smth_win = 0;
        box = 1;
        for smth_win = [5 10]
            for trainCorrect = 0:1
                runModelsDecPred_decodedPos('CA1', istheta, spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 0);
                runModelsDecPred_decodedPos('V1', istheta,  spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 0);
                
                runModelsDecPred_decodedPos('CA1', istheta, spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 1);
                runModelsDecPred_decodedPos('V1', istheta,  spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 1);
            end
        end
    end    
    %% 250ms windows
    istheta = 0;
    spike_smth_win = 250;
    for box = 0:1
        for smth_win = [5 10]
            for trainCorrect = 0:1
                runModelsDecPred_decodedPos('CA1', istheta, spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 0);
                runModelsDecPred_decodedPos('V1', istheta,  spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 0);
                
                runModelsDecPred_decodedPos('CA1', istheta, spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 1);
                runModelsDecPred_decodedPos('V1', istheta,  spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 1);
            end
        end
    end
    
%     %% Optimal smoothing
%     for istheta = 1:2
%         spike_smth_win = 0;
%         box = 1;
%         smth_win = [];
%         for trainCorrect = 0:1
%             runModelsDecPred_decodedPos('CA1', istheta, spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 0);
%             runModelsDecPred_decodedPos('V1', istheta,  spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 0);
%             
%             runModelsDecPred_decodedPos('CA1', istheta, spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 1);
%             runModelsDecPred_decodedPos('V1', istheta,  spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 1);
%         end
%     end
%     istheta = 0;
%     spike_smth_win = 250;
%     for box = 0:1
%         smth_win = [];
%         for trainCorrect = 0:1
%             runModelsDecPred_decodedPos('CA1', istheta, spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 0);
%             runModelsDecPred_decodedPos('V1', istheta,  spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 0);
%             
%             runModelsDecPred_decodedPos('CA1', istheta, spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 1);
%             runModelsDecPred_decodedPos('V1', istheta,  spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 1);
%         end
%     end
end