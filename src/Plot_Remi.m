function Plot_Remi(Y_vec, X_vec, theta_vec, vr_vec,theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)

% theta_set_disp = {'180n' '135n' '90n' '45n' '0' '45' '90' '135'};
% theta_set_disp = {'90n' '45n' '0' '45' '90' '135' '180n' '135n' };

load('proc_Haining_052213.mat')
% surface(Ur,a,cmx_tot)
% surface(Ur,a,cmy_tot)
% surface(Ur,a,clv_tot)
% surface(Ur,a,cdv_tot)

y_lim_CLv = [-0.5 0.5];
y_lim_CDv = [-0.7 0.5];
y_lim_Cmy = [-1 03];
% plot all CLv
for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        figure
        for k = 1:length(a)
            if k<3
                h(k) = plot(1./Ur, clv_tot(k+6,:));
                hold on
                text(1./Ur(16),clv_tot(k+6,15),theta_set_disp(k),'Color',col(k,:))
            else
                h(k) = plot(1./Ur, clv_tot(k-2,:));
                hold on
                text(1./Ur(16),clv_tot(k-2,15),theta_set_disp(k),'Color',col(k,:))
            end
            set(h(k),'Color',col(k,:))
            set(h(k),'Line','-')
            set(h(k),'Marker',shp{k})
            ylim(y_lim_CLv)
            xlim(x_lim)
            hold on
            
        end
        plot(0.125*[1 1],y_lim_CLv,'--r')
        plot(0.25*[1 1],y_lim_CLv,'--r')
        suptitle([marker ' CF CLv @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
        legend(theta_set_disp)
        legend1 = legend(gca,'show');
        uistack(legend1, 'top');
        %         set(legend1,'Position',[10 20 100 200]);
        saveas(gcf,[outputDir filesep 'CLv' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(j)) marker '.jpg'])
%         close
    end
end

% plot all CDv
for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        figure
        for k = 1:length(a)
            if k<3
                h(k) = plot(1./Ur, cdv_tot(k+6,:));
                hold on
                text(1./Ur(16),cdv_tot(k+6,15),theta_set_disp(k),'Color',col(k,:))
            else
                h(k) = plot(1./Ur, cdv_tot(k-2,:));
                hold on
                text(1./Ur(16),cdv_tot(k-2,15),theta_set_disp(k),'Color',col(k,:))
            end
            set(h(k),'Color',col(k,:))
            set(h(k),'Line','-')
            set(h(k),'Marker',shp{k})
            ylim(y_lim_CDv)
            xlim(x_lim)
            hold on
        end
        plot(0.125*[1 1],y_lim_CDv,'--r')
        plot(0.25*[1 1],y_lim_CDv,'--r')
        suptitle([marker ' IL CDv Coefficient @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
        legend(theta_set_disp)
        legend1 = legend(gca,'show');
        uistack(legend1, 'top');
        %         set(legend1,'Position',[10 20 100 200]);
        saveas(gcf,[outputDir filesep 'CDv' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(j)) marker '.jpg'])
%         close
    end
end
% plot all Cmy

for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        figure
        for k = 1:length(a)
            if k<3
                h(k) = plot(1./Ur, cmy_tot(k+6,:));
                hold on
                text(1./Ur(11),cmy_tot(k+6,10),theta_set_disp(k),'Color',col(k,:))
            else
                h(k) = plot(1./Ur, cmy_tot(k-2,:));
                hold on
                text(1./Ur(11),cmy_tot(k-2,10),theta_set_disp(k),'Color',col(k,:))
            end
            set(h(k),'Color',col(k,:))
            set(h(k),'Line','-')
            set(h(k),'Marker',shp{k})
            ylim(y_lim_Cmy)
            xlim(x_lim)
            hold on
        end
        plot(0.125*[1 1],y_lim_Cmy,'--r')
        plot(0.25*[1 1],y_lim_Cmy,'--r')
        hold on
        plot(bb_viva(:,1),bb_viva(:,3),'s--')
        suptitle([marker ' CF Added Mass Coefficient @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
        legend(theta_set_disp)
        legend1 = legend(gca,'show');
        uistack(legend1, 'top');
        %         set(legend1,'Position',[10 20 100 200]);
        saveas(gcf,[outputDir filesep 'Cmy' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(j)) marker '.jpg'])
%         close
    end
end


% plot all Cmx

for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        figure
        for k = 1:length(a)
            if k<3
                h(k) = plot(1./Ur, cmx_tot(k+6,:));
                hold on
                text(1./Ur(16),cmx_tot(k+6,15),theta_set_disp(k),'Color',col(k,:))
            else
                h(k) = plot(1./Ur, cmx_tot(k-2,:));
                hold on
                text(1./Ur(16),cmx_tot(k-2,15),theta_set_disp(k),'Color',col(k,:))
            end
            set(h(k),'Color',col(k,:))
            set(h(k),'Line','-')
            set(h(k),'Marker',shp{k})
            ylim(y_lim_Cmx)
            xlim(x_lim)
            hold on
        end
        hold on
        plot(0.125*[1 1],y_lim_Cmx,'--r')
        plot(0.25*[1 1],y_lim_Cmx,'--r')
        plot(basic_bare_inline_MT(:,1)./2,basic_bare_inline_MT(:,3),'s--')
        suptitle([marker ' IL Added Mass Coefficient @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
        legend(theta_set_disp)
        legend1 = legend(gca,'show');
        uistack(legend1, 'top');
        %         set(legend1,'Position',[10 20 100 200]);
        saveas(gcf,[outputDir filesep 'Cmx' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(j)) marker '.jpg'])
%         close
    end
end
