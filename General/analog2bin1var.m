function [C bins] = analog2bin1var(variable, aGrid, t)
% Usage: [C bins] = analog2bin1var(variable, aGrid, t)
% 'variable' is the variable to be discretised
% aGrid is the number of bins it should be discretised into (deft: 5)
% t is optional subset of logical operators
% 
% C: is the output matrix which is length(variable) x aGrid
% bisn are the values of each of the bins at the input
%
% % Aman Saleem
% % March 2014

if nargin< 3 | isempty(c)
    t = true(size(variable));
end
if nargin< 2 | isempty(aGrid)
    aGrid = 5;
end

[B bins] = normalise1var(variable, aGrid, t);

C = zeros(length(B),aGrid);

idx = [(1:length(B))' B];

goodPoints = ~isnan(B) & t;

for n = find(goodPoints)'
    C(idx(n,1),idx(n,2)) = 1;
end