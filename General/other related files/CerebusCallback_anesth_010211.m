function CerebusCallback_anesth(varargin)

SetDefaultDirs

global XY running counter Time DIRS 
persistent dirroot zdirroot nevroot nevfile resp pd starttime stimcounter stimTimes
local.Cerebus = 'C:\cerebus_data';

WinsockReceiver = varargin{1};

[foo,ReceivedString] = WinsockReceiver.GetData('',[],[]);

fprintf(1,'%s\n', ReceivedString);

[header,rest] = strtok(ReceivedString,' ');
[animal,rest] = strtok(rest);
[series,rest] = strtok(rest);
[expt,rest] = strtok(rest);

% create directories
dirroot  = fullfile(local.Cerebus,animal,series,expt);
zdirroot = fullfile(DIRS.Cerebus,animal,series,expt); 
if ~isdir(dirroot), mkdir(dirroot); end
if ~isdir(zdirroot), mkdir(zdirroot); end
nevroot = BuildFileName(dirroot,animal,series,expt);
nevfile = [nevroot '.nev'];

switch header
    case 'ExpStart' 
        % set up stuff
        resp = [];
        pd = audiorecorder(8000,8,1); % for photodiode
        Time = tic;
        stimcounter = 0;
        XY = cell(200000,1); % empty cell array for mouse counter
        %XY = zeros(200000,5); % empty array for mouse counter
        % start sampling
        cbmex('fileconfig',nevfile,'',1); % blackrock traces 
        running = 1; % ball coordinates
        starttime(1) = toc(Time);
        record(pd);  % photodiode
        starttime(2) = toc(Time);
    case 'BlockStart'
    case 'StimStart'
        stimcounter = stimcounter + 1;
        stimTimes(1,stimcounter) = toc(Time);
    case 'StimEnd'
        if(running)
            [rpt,rest] = strtok(rest);
            [stim,rest] = strtok(rest);
            stim = str2num(stim);
            rpt  = str2num(rpt);
            resp{stim,rpt} = cbmex('trialdata');
            stimTimes(2,stimcounter) = toc(Time);
        end
    case 'BlockEnd'
    case {'ExpEnd','ExpInterrupt'}
        if(~isempty(resp))
            % stop sampling
            cbmex('fileconfig',nevfile,'',0); % blackrock traces
            running = 0; counter = []; % ball coordinates
            stoptime(1) = toc(Time);
            stop(pd); % photodiode
            stoptime(2) = toc(Time);
            pdData = getaudiodata(pd,'int16');
            photodiode.starttime = starttime;
            photodiode.stoptime = stoptime;
            photodiode.stimTimes = stimTimes;
            photodiode.values = pdData;
            % remove unnecessary entries in XY and clear it
            for counter = length(XY):-1:1
               if isempty(XY{counter}), XY{counter} = []; end
               %if XY(counter,5)==0, XY(counter,:) = []; end % in case XY was not a cell array
            end
            % save data
            cd(dirroot)
            save([nevroot '.mat'],'resp');
            save([nevroot '_photodiode.mat'],'photodiode');
            save([nevroot '_ballCoords.mat'], 'XY');
            % copy all files onto zserver
            sourcefiles = [nevroot '*'];
            copyfile(sourcefiles,zdirroot,'f'); 
        end
%         cbmex('close')
end






    


