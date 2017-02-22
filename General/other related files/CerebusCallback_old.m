function CerebusCallback(varargin)

global nevroot nevfile resp running XY

WinsockReceiver = varargin{1};

% BytesTotal = varargin{3};
% TheEvent = varargin{end-1};
% EventName = varargin{end};

[foo,ReceivedString] = WinsockReceiver.GetData('',[],[]);

fprintf(1,'%s\n', ReceivedString);

[header,rest] = strtok(ReceivedString,' '); % find token in string
[date_animal,rest] = strtok(rest);
[iexp,rest] = strtok(rest);
[iblock,rest] = strtok(rest);

switch header
    case 'ExpStart' 
        nevroot = ['C:\cerebus_data\' date_animal '\u' sprintf('%03d',str2num(iexp))];
        mkdir('C:\cerebus_data\',date_animal)
    case 'BlockStart'
        running = true;
        resp = [];
        nevname = [nevroot '_'  sprintf('%03d',str2num(iblock))];
        nevfile = [nevname '.nev'];
        cbmex('fileconfig',nevfile,'',1);   %% start recording
        %running = cbmex('trialconfig',1); %% empty buffer and start recording
    case 'TrialStart'
    case 'TrialEnd'
        if(running)
            [istim,rest] = strtok(rest);
            istim = str2num(istim);
            resp{istim} = cbmex('trialdata',1); %% empty buffer and return recorded data
        end
    case 'BlockEnd'
        running = false;
        nevname = [nevroot '_'  sprintf('%03d',str2num(iblock))];
        if(~isempty(resp))
            cbmex('fileconfig',nevfile,'',0); %% stop sampling
            save([nevname '.mat'],'resp');
        end
        % save ball data
        eval(['save ' nevname '_ballCoords XY']);
    case {'ExpEnd','ExpInterrupt'}
        cbmex('close')
end




    


