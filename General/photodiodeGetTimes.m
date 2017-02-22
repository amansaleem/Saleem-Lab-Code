function [dig_pho,n_pho]=photodiodeGetTimes(animal,iseries,iexp)

% animal='M120307_BALL';
% iseries=322;
% iexp=106;

addpath \\zserver\Code\Spikes
global pepNEV;
global DIRS; 

SetDefaultDirs;


rsamp=1;

CHANNELS_ORDER = MichiganGetLayout(animal,iseries,33);
chan=find(CHANNELS_ORDER==33);


expName1 = [animal '_' num2str(iseries) '_' num2str(iexp) '.ns5'];
ns5file1 = ['\\ZSERVER\Data\Cerebus\' animal filesep num2str(iseries) filesep num2str(iexp) filesep expName1];

[~,SamplingRateInKHZ,nchan] = nsopen2(ns5file1);



pho = resample(double(pepNEV.ns.Data.data(chan,:)),1,rsamp);

%digitalise photodiode signal
n_pho = zeros(size(pho));
n_pho(pho>mean(pho)) = 1;
dig_pho = find(abs(diff(n_pho))>0.5)+1;

%phoTimes = dig_pho'./30000;

NameTimes = [animal '_' num2str(iseries) '_' num2str(iexp) '_screenTimes.mat'];
filenameScreentime=['\\ZSERVER\Data\multichanspikes\' animal filesep num2str(iseries) filesep NameTimes];
save(filenameScreentime, 'dig_pho')

end