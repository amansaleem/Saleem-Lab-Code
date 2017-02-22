
%% Preplay bayes decoding analysis
%requires rates matrix - with no. rows equal to no. of ms in sleep session,
%no. columns equal to no of cells recorded on track. These matrices are
%obtained from Sleep_rates folder in the Bayesing decoding folder

%requires objects for each template - e.g. UCA, DCA, UNA, DNA, obtained
%from Object subfolder in Bayesing Decoding folder

%Require subject structures, e.g. R504, R505, R584, R1838, these contain
%times of spiking events (in ms), as well as filters for which events
%belong to which templates, moreover these structure contain filters for
%which cells out of all cells recorded on track are in each template, these
%structures can be obtained from Data subfolder of the Bayes decoding
%folder.
subjects = {'R505'};
sleepSess = {'s2','s1'};

for r = 1:length(subjects);
    animal = subjects{r};
    for s = 1:length(sleepSess);
        sleep = sleepSess{s};

    
%         animal = 'R504';
%         sleep = 's1';
        disp(['Rat: ',animal,' sleep session: ',sleep]);
        
        %load in sleep rates
        cd 'C:\Users\pals\Dropbox\Decoding data\sleep_rates'
        load([animal,'_',sleep]);
        
        %load in template objects (i.e. trained decoder)
        cd(['C:\Users\pals\Dropbox\Decoding data\Object\',animal]);
        load([animal,'UCA']);
        load([animal,'DCA']);
        load([animal,'UNA']);
        load([animal,'DNA']);
        
        %load in time indices of events
        cd 'C:\Users\pals\Dropbox\Decoding data\Data'
        load(animal,[sleep,'_events']);
        
        %% define key variables
        
        silence = 50;%length of pre and post event silence
        event_length = 300; %max event length
        smoothing_window = 100;%corresponds to 5ms
        numBins = 100; %number of spatial bins for decoder
        tBinSz = 0.001; %temporal bin size during sleep in seconds
        sBinSz = 1; %spatial bin size in cm
        if strmatch('R1838',animal)%min cells/event different for R1838
            min_cells = 4;
        else
            min_cells = 5;%min. number of cells in an event
        end
        
        if strmatch(sleep,'s1')
            events = s1_events;
        else
            events = s2_events;
        end
        
        %% Main bit
        %UCA
        UCApred = zeros(event_length,length(events));
        UCApost = zeros(event_length,numBins,length(events));
        UCAnPost = zeros(event_length,numBins,length(events));
        UCAbestRes = zeros(size(events));
        UCAbestSpd = zeros(size(events));
        UCAbestY = zeros(size(events));
        UCAevent = zeros(event_length,length(events));
        h=waitbar(0, 'Doing offline decoding for UCA ...');
        for eIdx = 1:length(events);
            event_start = events(eIdx) + (silence - 1);
            event_end = event_start + (event_length -1);
            UCAevent(:,eIdx) = event_start: event_end;
            [r c] = find(rates(UCAevent(:,eIdx),:));%find matrix indices of spikes in event
            if length(unique(c)) >= min_cells;%exclude events with less than 5 cells in it
                %do decoding
                UCAcells(:,eIdx) = length(unique(c));% #cells in an event
                [UCApred(:,eIdx),UCApost(:,:,eIdx),UCAnPost(:,:,eIdx)] = UCA.predictBayesDecoder(rates(UCAevent(:,eIdx),:),smoothing_window,'best');
                %do line-fitting
                [ UCAbestRes(eIdx), UCAbestSpd(eIdx), UCAbestY(eIdx) ] = fast_best_line3( UCAnPost(:,:,eIdx)' ); %faster version of line-fitting code
            else
                UCAevent(:,eIdx) = NaN;
                UCAcells(:,eIdx) = NaN;
                UCApred(:,eIdx) = NaN;
                UCApost(:,:,eIdx) = NaN;
                UCAnPost(:,:,eIdx) = NaN;
                UCAbestRes(eIdx) = NaN;
                UCAbestSpd(eIdx) = NaN;
                UCAbestY(eIdx) = NaN;
                
            end
            waitbar(eIdx/length(events),h);
        end
        close(h);
        
        results.UCA.event = UCAevent;
        results.UCA.cells = UCAcells;
        results.UCA.pred = UCApred;
        results.UCA.post = UCApost;
        results.UCA.nPost = UCAnPost;
        results.UCA.bestRes = UCAbestRes;
        results.UCA.bestGrd = UCAbestSpd;
        results.UCA.bestY = UCAbestY;
        clear UCA*;
        pack;
        
        %DCA
        DCApred = zeros(event_length,length(events));
        DCApost = zeros(event_length,numBins,length(events));
        DCAnPost = zeros(event_length,numBins,length(events));
        DCAbestRes = zeros(size(events));
        DCAbestSpd = zeros(size(events));
        DCAbestY = zeros(size(events));
        DCAevent = zeros(event_length,length(events));
        h=waitbar(0, 'Doing offline decoding for DCA ...');
        for eIdx = 1:length(events);
            event_start = events(eIdx) + (silence - 1);
            event_end = event_start + (event_length -1);
            DCAevent(:,eIdx) = event_start: event_end;
            [r c] = find(rates(DCAevent(:,eIdx),:));%find matrix indices of spikes in event
            if length(unique(c)) >= min_cells;%exclude events with less than 5 cells in it
                %do decoding
                DCAcells(:,eIdx) = length(unique(c));
                [DCApred(:,eIdx),DCApost(:,:,eIdx),DCAnPost(:,:,eIdx)] = DCA.predictBayesDecoder(rates(DCAevent(:,eIdx),:),smoothing_window,'best');
                %do line-fitting
                [ DCAbestRes(eIdx), DCAbestSpd(eIdx), DCAbestY(eIdx) ] = fast_best_line3( DCAnPost(:,:,eIdx)' ); %faster version of line-fitting code
            else
                DCAevent(:,eIdx) = NaN;
                DCAcells(:,eIdx) = NaN;
                DCApred(:,eIdx) = NaN;
                DCApost(:,:,eIdx) = NaN;
                DCAnPost(:,:,eIdx) = NaN;
                DCAbestRes(eIdx) = NaN;
                DCAbestSpd(eIdx) = NaN;
                DCAbestY(eIdx) = NaN;
                
            end
            waitbar(eIdx/length(events),h);
        end
        close(h);
        
        results.DCA.event = DCAevent;
        results.DCA.cells = DCAcells;
        results.DCA.pred = DCApred;
        results.DCA.post = DCApost;
        results.DCA.nPost = DCAnPost;
        results.DCA.bestRes = DCAbestRes;
        results.DCA.bestGrd = DCAbestSpd;
        results.DCA.bestY = DCAbestY;
        clear DCA*;
        pack;
        
        %UNA
        UNApred = zeros(event_length,length(events));
        UNApost = zeros(event_length,numBins,length(events));
        UNAnPost = zeros(event_length,numBins,length(events));
        UNAbestRes = zeros(size(events));
        UNAbestSpd = zeros(size(events));
        UNAbestY = zeros(size(events));
        UNAevent = zeros(event_length,length(events));
        h=waitbar(0, 'Doing offline decoding for UNA ...');
        for eIdx = 1:length(events);
            event_start = events(eIdx) + (silence - 1);
            event_end = event_start + (event_length -1);
            UNAevent(:,eIdx) = event_start: event_end;
            [r c] = find(rates(UNAevent(:,eIdx),:));%find matrix indices of spikes in event
            if length(unique(c)) >= min_cells;%exclude events with less than 5 cells in it
                UNAcells(:,eIdx) = length(unique(c));
                %do decoding
                [UNApred(:,eIdx),UNApost(:,:,eIdx),UNAnPost(:,:,eIdx)] = UNA.predictBayesDecoder(rates(UNAevent(:,eIdx),:),smoothing_window,'best');
                %do line-fitting
                [ UNAbestRes(eIdx), UNAbestSpd(eIdx), UNAbestY(eIdx) ] = fast_best_line3( UNAnPost(:,:,eIdx)' ); %faster version of line-fitting code
            else
                UNAevent(:,eIdx) = NaN;
                UNAcells(:,eIdx) = NaN;
                UNApred(:,eIdx) = NaN;
                UNApost(:,:,eIdx) = NaN;
                UNAnPost(:,:,eIdx) = NaN;
                UNAbestRes(eIdx) = NaN;
                UNAbestSpd(eIdx) = NaN;
                UNAbestY(eIdx) = NaN;
                
            end
            waitbar(eIdx/length(events),h);
        end
        close(h);
        
        results.UNA.event = UNAevent;
        results.UNA.cells = UNAcells;
        results.UNA.pred = UNApred;
        results.UNA.post = UNApost;
        results.UNA.nPost = UNAnPost;
        results.UNA.bestRes = UNAbestRes;
        results.UNA.bestGrd = UNAbestSpd;
        results.UNA.bestY = UNAbestY;
        clear UNA*;
        pack;
        
        %DNA
        DNApred = zeros(event_length,length(events));
        DNApost = zeros(event_length,numBins,length(events));
        DNAnPost = zeros(event_length,numBins,length(events));
        DNAbestRes = zeros(size(events));
        DNAbestSpd = zeros(size(events));
        DNAbestY = zeros(size(events));
        DNAevent = zeros(event_length,length(events));
        h=waitbar(0, 'Doing offline decoding for DNA ...');
        for eIdx = 1:length(events);
            event_start = events(eIdx) + (silence - 1);
            event_end = event_start + (event_length -1);
            DNAevent(:,eIdx) = event_start: event_end;
            [r c] = find(rates(DNAevent(:,eIdx),:));%find matrix indices of spikes in event
            if length(unique(c)) >= min_cells;%exclude events with less than 5 cells in it
                DNAcells(:,eIdx) = length(unique(c));
                %do decoding
                [DNApred(:,eIdx),DNApost(:,:,eIdx),DNAnPost(:,:,eIdx)] = DNA.predictBayesDecoder(rates(DNAevent(:,eIdx),:),smoothing_window,'best');
                %do line-fitting
                [ DNAbestRes(eIdx), DNAbestSpd(eIdx), DNAbestY(eIdx) ] = fast_best_line3( DNAnPost(:,:,eIdx)' ); %faster version of line-fitting code
            else
                DNAevent(:,eIdx) = NaN;
                DNAcells(:,eIdx) = NaN;
                DNApred(:,eIdx) = NaN;
                DNApost(:,:,eIdx) = NaN;
                DNAnPost(:,:,eIdx) = NaN;
                DNAbestRes(eIdx) = NaN;
                DNAbestSpd(eIdx) = NaN;
                DNAbestY(eIdx) = NaN;
                
            end
            waitbar(eIdx/length(events),h);
        end
        close(h);
        
        results.DNA.event = DNAevent;
        results.DNA.cells = DNAcells;
        results.DNA.pred = DNApred;
        results.DNA.post = DNApost;
        results.DNA.nPost = DNAnPost;
        results.DNA.bestRes = DNAbestRes;
        results.DNA.bestGrd = DNAbestSpd;
        results.DNA.bestY = DNAbestY;
        
        cd 'C:\Users\pals\Dropbox\Decoding data\Results'
        save([animal,'_',sleep],'results','-v7.3');
    end
end
