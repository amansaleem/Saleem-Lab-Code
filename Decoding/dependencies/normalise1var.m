function [B bins] = normalise1var(A, aGrid, c);

% Usage: [B bins] = normalise1var(A, aGrid, c);
% 
% [B] = normalise1var(A) or [B] = normalise1var(A,[],c);
% Normalize the variable so it goes from 0 to 1
% 
% [B bins] = normalise1var(A, aGrid, c);
% Normalize the variable so it is discretized into 'aGrid' bins with centres 'bins'
% 'c' is optional logical array, given this the discretizasion is limited to the specified subset 

if nargin<3
    c = true(size(A));
else
    c = logical(c);
end

maxA = max(A(c));
minA = min(A(c));

if nargin<2 | isempty(aGrid)
    B = (A-minA)/(maxA-minA);
    bins = [];
else
    bins = minA:((maxA-minA)/(aGrid -1+eps)):maxA;
    A = (A-minA)/(maxA-minA);
    B = 1+floor(aGrid*A/(1+eps));
end