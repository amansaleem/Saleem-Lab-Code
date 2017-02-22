global pepNEV

for n = 1:14
    ext = num2str(n);
    if size(ext,2)==1
        ext = ['0' ext];
    end
%     [nevopen_outcome,foo,foo] = nevopen(['M101201_BALL_spont0' ext '.nev']);
    [nsopen_outcome,SamplingRateInKHZ,nchan] = nsopen(['M101201_BALL_spont0' ext '.ns5']);
    data = pepNEV.ns.Data.data;
    save(['M101201_BALL_spont0' ext]);
end
