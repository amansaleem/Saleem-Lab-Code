function out = smthInTime(in, sampFreq, smthWin, type, subset)
%
% Usage: out = smthInTime(in, sampFreq, smthWin, type, subset)
% to smooth in time the signal 'in'
%
% sampFreq = the sampling freq
% smthWin  = the half width of the smoothing window (in ms), 0 gives no
% smoothing
% type     = Shape argument of conv (optional), 'same' is default
% subset   = logical of region to be considered same size as input, everything else comes out
% as NaNs
%
% Aman Saleem
% November 2013

if nargin<5
    subset = true(size(in));
end
if nargin<4
    type = 'same';
end

if smthWin==0
    out = in;
elseif ~strcmp(type,'box')
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(smthWin/samp_int);
    
    sGrid = max([100 win*5]);
    s = (-sGrid:sGrid);
    
    sfilt = (1./(win*(sqrt(2*pi)))).*exp(-(s.^2)./(2*(win^2)));
    
    % padding the inputs
    pad = zeros(size(sfilt));
    pad_length = length(pad);
    if size(in,1)~=1
        flip_again = true;
        in = in';
        subset = subset';
    else
        flip_again = false;
    end
    pad_in = [pad in pad];
    pad_subset = [pad double(subset) pad];
    
    pad_out         = conv(pad_in, sfilt, type);
    norm_subset     = conv(pad_subset, sfilt, type);
    
    pad_out = pad_out./norm_subset;
    out = pad_out((pad_length+1) : (end-pad_length));
    out(~subset) = NaN;
    
    if flip_again
        out = out';
    end
else
    nanentries = isnan(in);
    subset = ~nanentries & subset;
    
    in(nanentries) = 0;
    samp_int = 1000/sampFreq;
    
    win = round(smthWin/samp_int);
    
%     sGrid = max([100 win*5]);
%     s = (-sGrid:sGrid);
%     
%     sfilt = (1./(win*(sqrt(2*pi)))).*exp(-(s.^2)./(2*(win^2)));
    
    % padding the inputs
    pad = zeros(size(win));
    pad_length = length(pad);
    if size(in,1)~=1
        flip_again = true;
        in = in';
        subset = subset';
    else
        flip_again = false;
    end
    pad_in = [pad in pad];
    pad_subset = [pad double(subset) pad];
    
    pad_out         = fastsmooth(pad_in,     win, 1, 1)*win;
    norm_subset     = fastsmooth(pad_subset, win, 1, 1)*win;
    
    pad_out = pad_out./norm_subset;
    out = pad_out((pad_length+1) : (end-pad_length));
    out(~subset) = NaN;
    
    if flip_again
        out = out';
    end
end