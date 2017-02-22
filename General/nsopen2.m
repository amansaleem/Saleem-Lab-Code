function [r,SamplingRateInKHZ,nChans] = nsopen2(filename)
% Created DLR who knows when
% 2010-03-16 AZ Tweaked some code, added memmapfile at end
global pepNEV;

% if ~isempty(pepNEV)
%     nevclose;
% end
% filename = sprintf('%s.ns%u',filename,d);
fid = fopen(filename, 'r', 'ieee-le');
if fid == -1, 
    r = -1;
    return;
end
% AZ 2009-03-04
pepNEV.index.fid = fid;

%%%%%%%%%%%  READ HEADER (see Cerebus manual for info/naming)  %%%%%%%%%%%

if fseek(fid,0,'bof')
   r = -1;
   return;
end

fileid = fread(fid,8,'*char')';
if ~strcmp(fileid,'NEURALSG')
   r = -1;
   return;
end

label = fread(fid,16,'*char')';
a = findstr(label,' kS/s');
label = str2double(label(1:a-1));
% if label == 0 or isempty(a), "LFP Low" sampling rate
% if not LFP, label and SamplingRateInKHZ should be equal
SamplingRateInKHZ = 30/fread(fid,1,'uint32=>double');
disp([filename,sprintf(' sampled at %ukHz',SamplingRateInKHZ)]);

nChans = fread(fid,1,'uint32=>double');
pepNEV.nchan = nChans;

pepNEV.ids = fread(fid,nChans,'uint32=>double');
% This is the end of the header / beginning of the data section

nbytesInHeader = 32 + 4*nChans;
% Specify DATA file format
fseek(fid,0,'eof');
nbytes = ftell(fid);
fclose(fid);
% fseek(fid,nbytesInHeader,'bof');

nsformat = { 'int16' [nChans (nbytes-nbytesInHeader)/(nChans*2)] 'data';};
% Open file
pepNEV.ns = memmapfile(filename,'Format',nsformat,'Offset',nbytesInHeader);

r = 1;