function CSD_Ch=LFP2iCSD(LFP_Ch)
RepDIM=size(LFP_Ch,1);
Espacing=50e-6;
el_pos=Espacing:Espacing:RepDIM*Espacing;
CSDgaussStd=2*Espacing;

Fcs = F_cubic_spline(el_pos,500e-6,0.3,0.3);
[zs,CSD_cs] = make_cubic_splines(el_pos,LFP_Ch,Fcs);
[zs,CSD_cs] = gaussian_filtering(zs,CSD_cs,CSDgaussStd,5*CSDgaussStd);

% newpos=zeros(1,RepDIM);
% for ch=1:RepDIM
%     newpos(ch)=find(zs<el_pos(ch),1,'last');
% end
% 
% CSD_Ch=round(CSD_cs(newpos,:));
CSD_Ch = CSD_cs;
end