function [out]              = setSpdLims(in, low, high)
% low and high and in should be in linear units... 
% This converts to log.
% **** Make the conversion to log scale optional ****
        out = zeros(size(in));
        k = exp(low);
        out((in+k)<=0) = NaN;
        out((in+k)>0)  = log((in((in+k)>0)+k)*0.06);
        out( out < low ) = low;
        out( out > high ) = high;
end