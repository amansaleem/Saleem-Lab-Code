function out = autocorr_as(input, numBins);

if nargin<2
    numBins = length(input)-1;
end

pad = nan*ones(1,numBins+100);
padded_input = [input pad];

for delays = 1:numBins
    out(delays) = ...nanmean([...
        corr(padded_input',circshift(padded_input',delays), 'rows','pairwise')...
        ...,corr(input',circshift(input',-delays), 'rows','complete')...
        ...])...
        ;
end
if nargout<1
    stem(1:numBins, out);
end