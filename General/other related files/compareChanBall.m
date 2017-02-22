for n = 1:14
    clear filt_data ext ichan ichan1 ichan2 xcf lags 
    ext = num2str(n);
    if size(ext,2)==1
        ext = ['0' ext];
    end
    
    load(['\\ZSERVER\Data\Cerebus\M101201_BALL\4\spont\M101201_BALL_spont0' ext]);
    
    for ichan = 1:size(data,1)
        [b, a] = ellip(2,0.1,40,[300*2/30000], 'low');
        filt_data(ichan,:) = filtfilt(b,a,double(data(ichan,:)));
    end
    
    for ichan1 = 1:size(data,1)
        ichan1
        for ichan2 = 1:size(data,1)
            [xcf{ichan1, ichan2}, lags{ichan1,ichan2}] = crosscorr(downsample(double(filt_data(ichan1,:)),2), downsample(double(filt_data(ichan2,:)),2),2000);
            disp('.')
        end
    end
    
    clear b a 
    save(['\\ZSERVER\Data\Cerebus\M101201_BALL\4\spont\M101201_BALL_spont0' ext]);
    
end