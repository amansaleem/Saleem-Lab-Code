% [Pred] = runModelsDecPred(...
%     area, istheta, spike_smth_win, ...
%     box, numSampleBins, smth_win, ...
%     trainCorrect, trainAllCorrect)
for numBins = 250%[200 50]
    %% Theta,
    for istheta = [2 1]%:2
        spike_smth_win = 0;
        box = 1;
        for smth_win = [8 10]%[5 10]
            for trainCorrect = 0%:1
                runModelsDecPred('CA1', istheta, spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 0);
                runModelsDecPred('V1', istheta,  spike_smth_win,...
                    box, numBins, smth_win, trainCorrect, 0);
            end
        end
    end
end
% for numBins = [50 200]
%     %% Theta,
%     for istheta = 1:2
%         spike_smth_win = 0;
%         box = 1;
%         for smth_win = [5 10]
%             for trainCorrect = 0:1
%                 runModelsDecPred('CA1', istheta, spike_smth_win,...
%                     box, numBins, smth_win, trainCorrect, 0);
%                 runModelsDecPred('V1', istheta,  spike_smth_win,...
%                     box, numBins, smth_win, trainCorrect, 0);
%                 
%                 runModelsDecPred('CA1', istheta, spike_smth_win,...
%                     box, numBins, smth_win, trainCorrect, 1);
%                 runModelsDecPred('V1', istheta,  spike_smth_win,...
%                     box, numBins, smth_win, trainCorrect, 1);
%             end
%         end
%     end
%     
%     %% 250ms windows
%     istheta = 0;
%     spike_smth_win = 250;
%     for box = 0:1
%         for smth_win = [5 10]
%             for trainCorrect = 0:1
%                 runModelsDecPred('CA1', istheta, spike_smth_win,...
%                     box, numBins, smth_win, trainCorrect, 0);
%                             runModelsDecPred('V1', istheta,  spike_smth_win,...
%                                 box, numBins, smth_win, trainCorrect, 0);
%                 
%                 runModelsDecPred('CA1', istheta, spike_smth_win,...
%                     box, numBins, smth_win, trainCorrect, 1);
%                             runModelsDecPred('V1', istheta,  spike_smth_win,...
%                                 box, numBins, smth_win, trainCorrect, 1);
%             end
%         end
%     end
%     
%     %% Optimal smoothing
%     for istheta = 1:2
%         spike_smth_win = 0;
%         box = 1;
%         smth_win = [];
%         for trainCorrect = 0:1
%             runModelsDecPred('CA1', istheta, spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 0);
%                     runModelsDecPred('V1', istheta,  spike_smth_win,...
%                         box, numBins, smth_win, trainCorrect, 0);
%             
%             runModelsDecPred('CA1', istheta, spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 1);
%                     runModelsDecPred('V1', istheta,  spike_smth_win,...
%                         box, numBins, smth_win, trainCorrect, 1);
%         end
%     end
%     istheta = 0;
%     spike_smth_win = 250;
%     for box = 0:1
%         smth_win = [];
%         for trainCorrect = 0:1
%             runModelsDecPred('CA1', istheta, spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 0);
%                     runModelsDecPred('V1', istheta,  spike_smth_win,...
%                         box, numBins, smth_win, trainCorrect, 0);
%             
%             runModelsDecPred('CA1', istheta, spike_smth_win,...
%                 box, numBins, smth_win, trainCorrect, 1);
%                     runModelsDecPred('V1', istheta,  spike_smth_win,...
%                         box, numBins, smth_win, trainCorrect, 1);
%         end
%     end
% end