function balldata = getBallTuning(animal, iseries, iexp)
% program to read ball data and the data from units and start to compare
% different recordings

setDefaultDirs;
global DIRS


prot = ProtocolLoad(animal,iseries,iexp);
expt = ExptLoadCerebusTraces(animal,iseries, iexp);
units = UnitLoad(DIRS.spikes, animal,iseries, iexp);

balldata = getNewBallData(animal,iseries,iexp);

balldata.movement = cell(size(prot.seqnums));
balldata.meanMovement = zeros(size(prot.seqnums));
idxSize = size(prot.seqnums,1)*size(prot.seqnums,2);
mm = zeros(1,idxSize);

for idx = 1:idxSize
    balldata.movement{idx} = balldata.data{prot.seqnums(idx)}(1:4,:);
    balldata.movement{idx}(5,:) = round(sqrt(mean(balldata.movement{idx}.^2)));
    balldata.meanMovement(idx) = round(mean(balldata.movement{idx}(5,:)));
    balldata.mm(idx) = round(mean(balldata.movement{idx}(5,:)));
end

figure;
hist(balldata.mm, 15); 
xlabel('Movement','fontsize',12);
ylabel('#','fontsize',12);
title('Histogram: to pick a threshold','fontsize',14);

balldata.thres = input('Enter the threshold for movement:  ');
balldata.ismoving = (balldata.meanMovement > balldata.thres);
close

units_run = units;
units_stop = units;

for unitIdx = 1:length(units)
    if (units(unitIdx).ichan < 10000)
        units_run(unitIdx).ichan = str2num([num2str(units(unitIdx).ichan) '11']);
        units_stop(unitIdx).ichan = str2num([num2str(units(unitIdx).ichan) '00']);
        
        for n = 1:size(units(unitIdx).stimdurs,1)
            for m = 1:size(units(unitIdx).stimdurs,2)
                if balldata.ismoving(n,m) == 1
                    units_run(unitIdx).stimdurs(n,m) = units(unitIdx).stimdurs(n,m);
                    units_run(unitIdx).spiketimes{n,m} = units(unitIdx).spiketimes{n,m};
                    units_stop(unitIdx).spiketimes{n,m} = NaN;
                    units_stop(unitIdx).stimdurs(n,m) = NaN;
                else
                    units_stop(unitIdx).stimdurs(n,m) = units(unitIdx).stimdurs(n,m);
                    units_stop(unitIdx).spiketimes{n,m} = units(unitIdx).spiketimes{n,m};
                    units_run(unitIdx).spiketimes{n,m} = NaN;
                    units_run(unitIdx).stimdurs(n,m) = NaN;
                end
            end
        end
        unitSave(units_run(unitIdx), DIRS.spikes)
        unitSave(units_stop(unitIdx), DIRS.spikes)
    end
end



save([DIRS.ball filesep 'Recording' filesep animal filesep num2str(iseries) filesep num2str(iexp) filesep ...
    animal '_' num2str(iseries) '_' num2str(iexp) '_BallData.mat'], 'balldata');

