function out = smthInTime(in, sampFreq, smthWin, type, subset)
%
% Usage: out = smthInTime(in, sampFreq, smthWin, type subset)
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

nanentries = isnan(in);

subset = ~nanentries & subset;

if smthWin==0
    out = in;
else
    samp_int = 1000/sampFreq;
    
    win = round(smthWin/samp_int);
    
    sGrid = max([100 win*5]);
    s = (-sGrid:sGrid);
    
    sfilt = (1./(win*(sqrt(2*pi)))).*exp(-(s.^2)./(2*(win^2)));
    
    out = conv(in, sfilt, type);
end