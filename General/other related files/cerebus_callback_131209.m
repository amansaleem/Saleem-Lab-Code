function cerebus_callback(varargin)

global nevroot nevfile resp running;

WinsockReceiver = varargin{1};

BytesTotal = varargin{3};
TheEvent = varargin{end-1};
EventName = varargin{end};

[foo, ReceivedString] = WinsockReceiver.GetData('',[],[]);

fprintf(1,'%s\n', ReceivedString);

[header,rest] = strtok(ReceivedString,' '); % find token in string
[animal,rest] = strtok(rest);
[series,rest] = strtok(rest);
[expt,rest] = strtok(rest);

switch header
    case 'ExpStart'
        running = 1;
        resp = [];
        nevroot = ['c:\cerebus_data\' animal '\u' sprintf('%03d',str2num(series)) '_'  sprintf('%03d',str2num(expt))];
        nevfile = [nevroot '.nev'];
        cbmex('fileconfig',nevfile,'',1);   %% start sampling
    case 'BlockStart'
    case 'StimStart'
    case 'StimEnd'
        if(running)
            [rpt,rest] = strtok(rest);
            [stim,rest] = strtok(rest);
            stim = str2num(stim);
            rpt  = str2num(rpt);
            resp{stim,rpt} = cbmex('trialdata');
            plot_raster(resp(stim,:));
        end
    case 'BlockEnd'
    case {'ExpEnd', 'ExpInterrupt'}
        running = 0;
        if(~isempty(resp))
            cbmex('fileconfig',nevfile,'',0);
            save([nevroot '.mat'],'resp');
        end
end


    

