function [rateMap I_skaggs R traj spikeRate normTraj posXY] = getVRPlaceMaps_fromRates(rateMap, discreet_steps, replay_in)

% global REPLAY
% global DIRS
% 
% SetDefaultDirs
% 
if nargin<3
    REPLAY = 0;
else
    REPLAY = replay_in;
end

if nargin<2
    discreet_steps = 39;
else
    discreet_steps = discreet_steps - 1;
end
% 
% [expt VRdata chans] = getVRspikeTimes(animal,series,iexp);
% 
% if REPLAY
%     dDIRname = [DIRS.ball filesep animal filesep num2str(series)];
%     [VRdata_o] = load([dDIRname filesep VRdata.EXP.replayed_filename]);
% end
% 
% spikeRate   = expt.spikeRate;
% traj        = VRdata.TRIAL.traj(1:max(find(VRdata.TRIAL.traj))-1);
% 
% maxTime = min(length(traj),length(spikeRate{1,1}));
% 
% posX        = VRdata.TRIAL.posdata(1:maxTime,1);
% posY        = VRdata.TRIAL.posdata(1:maxTime,3);
% 
% if REPLAY
%     trialIdx    = VRdata_o.TRIAL.trialIdx(1:maxTime);
% else
%     trialIdx    = VRdata.TRIAL.trialIdx(1:maxTime);
% end
% 
% if length(traj) ~= length(spikeRate{1,1})
%     disp(['!!!!!WARNING!!!!!!: trajectory has ' num2str(length(traj) - length(spikeRate{1,1})) ' more frames than screen refreshes'])
% end
% 
% maxTime = min(length(traj),length(spikeRate{1,1}));
% 
% normTraj = double((traj - min(traj))./max(traj));
% normTraj = round(normTraj*discreet_steps)./discreet_steps;
% 
% normTraj = normTraj(1:maxTime);
% places = unique(normTraj);
% trials = unique(trialIdx);
% 
for icorner = 1:4
    corIdx(icorner) = ceil(icorner*(discreet_steps+1)/4);
end

pXY = zeros(discreet_steps+1,2);

pXY(1:corIdx(1),1) = 1;
pXY(1:corIdx(1),2) = 1:corIdx(1);

pXY(corIdx(1)+1:corIdx(2),2) = corIdx(1);
pXY(corIdx(1)+1:corIdx(2),1) = 1:corIdx(1);

pXY(corIdx(2)+1:corIdx(3),1) = corIdx(1);
pXY(corIdx(2)+1:corIdx(3),2) = corIdx(1):-1:1;

pXY(corIdx(3)+1:corIdx(4),2) = 1;
pXY(corIdx(3)+1:corIdx(4),1) = corIdx(1):-1:1;

map = ones(corIdx(1), corIdx(1),3);
map_exist = zeros(corIdx(1), corIdx(1));

ifigure = 1;

for ichan = 1:size(rateMap,1)
    for iprot = 1:size(rateMap,2)
        if ~isempty(rateMap(ichan,iprot))
            lambda_x = rateMap(ichan,iprot).meanRate;
            lambda = rateMap(ichan,iprot).Lambda;
            log_lambda = log2(lambda_x./lambda);
            
            %             lambda
            I_skaggs(ifigure) = sum( lambda_x(~isinf(log_lambda)) .* log2(lambda_x(~isinf(log_lambda))./lambda) ) / lambda;
            
            rateMap(ichan,iprot).I_skaggs(ifigure) = I_skaggs(ifigure) ;
            R(ifigure) =  nanmean(nansem(rateMap(ichan,iprot).rate'));
            
            rateMap(ichan,iprot).R = R(ifigure);
            
%             if isfinite(real(chans(ichan).cellids(iprot)))
                figure(ifigure);
                subplot(2,1,1);
                ifigure = ifigure + 1;
                % % %
                hold on;
                %                 for n = 1:4
                %                     fill([round(n*discreet_steps/4)-1 round(n*discreet_steps/4)+1 round(n*discreet_steps/4)+1 round(n*discreet_steps/4)-1],...
                %                         [0 0 max(mean(rateMap(ichan,iprot).rate))*1.2 max(mean(rateMap(ichan,iprot).rate))*1.2],[0.7 0.7 0.7]);
                %                 end
                if REPLAY
                    errorbar(nanmean(rateMap(ichan,iprot).rate,2),nansem(rateMap(ichan,iprot).rate'),'r');
                else
                    errorbar(nanmean(rateMap(ichan,iprot).rate,2),nansem(rateMap(ichan,iprot).rate'),'b');
                end
                hold on;
                
                axis tight
                ylims = get(gca, 'YLim');
                for icorner = 1:4
                    cline =  line([ceil(icorner*(discreet_steps+1)/4) ceil(icorner*(discreet_steps+1)/4)], ylims);
                    set(cline, 'Color', 'b', 'linestyle', '--');
                end
                
%                 title(['Animal: ' num2str(animal) ' Date: ' num2str(series) ' Session: ' num2str(iexp) ' Channel: ' num2str(chans(ichan).channame) ' Prot: ' num2str(chans(ichan).cellids(iprot))], 'fontsize', 14)
                xlabel('Position on track' , 'fontsize', 14);
                ylabel('Firing rate (spks/s)', 'fontsize', 14);
                
                %                 figure(ifigure);
%                 subplot(1,2,2);
%                 
%                 % % %
%                 hold on;
%                 %                 for n = 1:4
%                 %                     fill([round(n*discreet_steps/4)-1 round(n*discreet_steps/4)+1 round(n*discreet_steps/4)+1 round(n*discreet_steps/4)-1],...
%                 %                         [0 0 max(mean(rateMap(ichan,iprot).rate))*1.2 max(mean(rateMap(ichan,iprot).rate))*1.2],[0.7 0.7 0.7]);
%                 %                 end
%                 for itrial = 1:length(trials)
%                     if REPLAY
%                         plot(rateMap(ichan,iprot).rate(:,itrial),'color',[1 0.7 0.7]);
%                     else
%                         plot(rateMap(ichan,iprot).rate(:,itrial),'color',[0.7 0.7 1]);
%                     end
%                 end
%                 
%                 if REPLAY
%                     plot(nanmean(rateMap(ichan,iprot).rate,2),'r', 'linewidth', 2);
%                 else
%                     plot(nanmean(rateMap(ichan,iprot).rate,2),'b', 'linewidth', 2);
%                 end
%                 
%                 hold on;
%                 
%                 axis tight
%                 ylims = get(gca, 'YLim');
%                 for icorner = 1:4
%                     cline =  line([ceil(icorner*(discreet_steps+1)/4) ceil(icorner*(discreet_steps+1)/4)], ylims);
%                     set(cline, 'Color', 'b', 'linestyle', '--');
%                 end
%                 
%                 title(['Animal: ' num2str(animal) ' Date: ' num2str(series) ' Session: ' num2str(iexp) ' Channel: ' num2str(chans(ichan).channame) ' Prot: ' num2str(chans(ichan).cellids(iprot))], 'fontsize', 14)
%                 xlabel('Position on track' , 'fontsize', 14);
%                 ylabel('Firing rate (spks/s)', 'fontsize', 14);
                                subplot(2,1,2);
                                mRate = nanmean(rateMap(ichan,iprot).rate,2);
                                maxRate = max(mRate);
                
                                mRate = mRate./maxRate;
                                mRate = round(mRate*63)+1;
                
                                k = jet;
                                m3Rate = k(mRate,:);
                
                                for iplace = 1:length(mRate)
                                    map(pXY(iplace,1), pXY(iplace,2),:) =  m3Rate(iplace,:);
                                end
                
                                imagesc(map); colormap(jet);
                                axis equal
                                axis off
                % %                 hold on;
                
%             end
        end
        
    end
end
% figure;
% % % %
% hold on;
% %                 for n = 1:4
% %                     fill([round(n*discreet_steps/4)-1 round(n*discreet_steps/4)+1 round(n*discreet_steps/4)+1 round(n*discreet_steps/4)-1],...
% %                         [0 0 max(mean(rateMap(ichan,iprot).rate))*1.2 max(mean(rateMap(ichan,iprot).rate))*1.2],[0.7 0.7 0.7]);
% %                 end
% errorbar(nanmean(rateMap(ichan,iprot).P_place,2),nansem(rateMap(ichan,iprot).P_place'));
% hold on;

axis tight
ylims = get(gca, 'YLim');
for icorner = 1:4
    cline =  line([ceil(icorner*(discreet_steps+1)/4) ceil(icorner*(discreet_steps+1)/4)], ylims);
    set(cline, 'Color', 'b', 'linestyle', '--');
end

title(['Animal: ' num2str(animal) ' Date: ' num2str(series) ' Session: ' num2str(iexp) ' Prob of place'])
a = 1;


clear global REPLAY