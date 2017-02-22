function inspectOutput(expt,fname,chan,thresh)
% inspectOutput is used to inspect output data from cerebus (nev files, ns* files, and photodiode output)
% inspectOutput(expt,fname,chan,type,thresh)
% expt = e.g. 'M100311_0'
% fname = e.g. 'u001_001.nev'
% chan = e.g. 10
% thresh = (optional) threshold to rethreshold nev data to; 32 is a good default
% 2010-03 MS

global pepNEV

addpath('\\Zserver\Code\CerebusTools')
fileName = ['C:\cerebus_data\' expt '\' fname];
%fileName = ['\\zserver\Code\MouseRoom\NeuroBall\' expt '\' fname];
%fileName = ['\\zserver\Data\Cerebus\' expt '\' fname];

if nargin< 4, thresh = []; end

if ~isempty(strfind(fname,'nev'))
    
    nevopen_outcome = nevopen(fileName); %% open nev file
    PD_times  = pepNEV.sync.timestamps;  %% photodiode timestamps
    PD_values = pepNEV.sync.digital;     %% photodiode values
    [Timestamps, waveforms] = nevwaves(chan); %% read nev file
    nevclose;
    
    if thresh %% rethreshold and realign 
    Xraw = waveforms;
    [Xrealigned,iXraw] = changeThreshold(thresh,Xraw,[],nevopen_outcome);
     % Xrealigned is Xraw(:,~iXraw) combined with realigned Xraw(:,iXraw)
    % Xrealigned(:,iXraw) are all the traces (realigned) above threshold
    % Xraw(:,iXraw) are all the raw traces (not realigned) above threshold
    end
    
    % plot spike times
    figure(1); 
    if thresh, plot(Xrealigned(chan,:));
    else stem(Timestamps,ones(1,length(Timestamps))), ylim([-1 2]);
    end
    title('nev data');
    
    % plot photodiode signal
    figure(2);
    stem(PD_times,PD_values,'.');
    title('photodiode');
    
elseif ~isempty(strfind(fname,'ns'))
    
    nsopen_outcome = nsopen(fileName); %% open ns* file
    waveform = pepNEV.ns.Data.data(:,chan);
    nevclose;
    
    % plot continuous data
    figure; plot(waveform);
    title('ns data');
       
end





