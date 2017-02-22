% load Posterior_all_150327_trainedCV_250_noInter
% CA = Posterior_all;
% CA = CA([1:6 8]);
% 
% load Posterior_all_V1_150601_trainedCV_250.mat
% V1 = Posterior_all;

for mid = 5:10:45
    for iseries = 1:7
        start = mid-5;
        stop  = mid+5;
        X_ML_l = CA(iseries).X_low_orig';
        V_ML_l = V1(iseries).MAP.low(1:length(X_ML_l));
        C_ML_l = CA(iseries).MAP.low(1:length(X_ML_l));
        
        temp = (X_ML_l>start & X_ML_l<stop);
        temp_s = temp & ((C_ML_l>start & C_ML_l<stop) | (V_ML_l>start & V_ML_l<stop));
        subplot(231)
        plot(C_ML_l(temp)+randn(1,sum(temp)), V_ML_l(temp)+randn(1,sum(temp)),'b.')
        hold on;
%         circle([mid mid],5,'r',3);
        circle([mid mid],10,'r--',3);
        line(xlim, [mid+10 mid+10],'color','r','linewidth',1.5)
        line(xlim, [mid-10 mid-10],'color','r','linewidth',1.5)
        line([mid-10 mid-10],ylim,'color','r','linewidth',1.5)
        line([mid+10 mid+10],ylim,'color','r','linewidth',1.5)
        axis equal; axis square; axis([0 50 0 50]); line(xlim,ylim);
        
        X_ML_n = CA(iseries).X_norm';
        V_ML_n = V1(iseries).MAP.norm(1:length(X_ML_n));
        C_ML_n = CA(iseries).MAP.norm(1:length(X_ML_n));
        
        temp = (X_ML_n>start & X_ML_n<stop);
        subplot(232)
        plot(C_ML_n(temp)+randn(1,sum(temp)), V_ML_n(temp)+randn(1,sum(temp)),'k.')
        hold on;
%         circle([mid mid],5,'r',3);
        circle([mid mid],10,'r--',3);
        line(xlim, [mid+10 mid+10],'color','r','linewidth',1.5)
        line(xlim, [mid-10 mid-10],'color','r','linewidth',1.5)
        line([mid-10 mid-10],ylim,'color','r','linewidth',1.5)
        line([mid+10 mid+10],ylim,'color','r','linewidth',1.5)
        axis equal; axis square; axis([0 50 0 50]); line(xlim,ylim);
        
        X_ML_h = CA(iseries).X_high_orig';
        V_ML_h = V1(iseries).MAP.high(1:length(X_ML_h));
        C_ML_h = CA(iseries).MAP.high(1:length(X_ML_h));
        
        temp = (X_ML_h>start & X_ML_h<stop);
        subplot(233)
        plot(C_ML_h(temp)+randn(1,sum(temp)), V_ML_h(temp)+randn(1,sum(temp)),'r.')
        hold on;
%         circle([mid mid],5,'b',3);
        circle([mid mid],10,'b--',3);
        line(xlim, [mid+10 mid+10],'color','b','linewidth',1.5)
        line(xlim, [mid-10 mid-10],'color','b','linewidth',1.5)
        line([mid-10 mid-10], ylim,'color','b','linewidth',1.5)
        line([mid+10 mid+10], ylim,'color','b','linewidth',1.5)
        axis equal; axis square; axis([0 50 0 50]); line(xlim,ylim);
%         pause;
    end
    
    for n = 1:6; subplot(2,3,n); hold off; end
    
    pause
end