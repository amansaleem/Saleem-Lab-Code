function [starts, stops] = getStartStop(t, minWidth)
% function to find the starts and stops of a signal. 
% Should be useful to split the data to remove nan's 0's etc.
% [starts, stops] = getStartStop(t)
% starts: start bins
% stops: stop bins
% t: logical array
% [starts, stops] = getStartStop(t, minWidth)
% minWidth: removes period that have a length less than minWidth

starts = find(t);
stops = find(t);
starts(find(diff(starts)==1)+1) = [];
stops(find(diff(stops)==1)) = [];

if nargin>1
    toRid = (stops-starts)<minWidth;
    starts(toRid) = [];
    stops(toRid) = [];
end