function plotDatabaseComp(ILCFHydroModelHZobj,ILCFHydroModelJMDobj,Y_vec, X_vec, theta_vec, vr_vec,...
    theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)

% plot all Cmy

for i = 1:length(X_vec)

    figure
    for k = 1:length(theta_vec)
        h = plot1D_f(ILCFHydroModelHZobj , 'Cmy', Y_vec(i), X_vec(i), theta_vec(k), vr_vec);
        set(h,'Color',col(k,:))
        set(h,'Line','-')
        set(h,'Marker',shp{k})
        ylim(y_lim_Cmy)
        xlim(x_lim)
        hold on
        text(1./vr_vec(5),getDataPoint(ILCFHydroModelHZobj,'Cmy', Y_vec(i), X_vec(i), theta_vec(k)/180, vr_vec(5)),theta_set_disp(k),'Color',col(k,:))
    end
    legend(theta_set_disp)
    for k = 1:length(theta_vec)
        h = plot1D_f(ILCFHydroModelJMDobj , 'Cmy', Y_vec(i), X_vec(i), theta_vec(k), vr_vec);
        set(h,'Color',col(k,:))
        set(h,'Line','--')
        set(h,'Marker',shp{k})
    end
    plot(0.125*[1 1],y_lim_Cmy,'--r')
    plot(0.25*[1 1],y_lim_Cmy,'--r')
    hold on
    plot(bb_viva(:,1),bb_viva(:,3),'s--')
    suptitle([marker ' CF Added Mass Coefficient @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(i))])
    legend1 = legend(gca,'show');
    uistack(legend1, 'top');
    %         set(legend1,'Position',[10 20 100 200]);
    saveas(gcf,[outputDir filesep 'Cmy' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) marker 'Comp.jpg'])
    saveas(gcf,[outputDir filesep 'Cmy' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) marker 'Comp.fig'])
    close
end


% plot all Cmx

for i = 1:length(X_vec)

    figure
    for k = 1:length(theta_vec)
        h = plot1D_f(ILCFHydroModelHZobj , 'Cmx', Y_vec(i), X_vec(i), theta_vec(k), vr_vec);
        set(h,'Color',col(k,:))
        set(h,'Line','-')
        set(h,'Marker',shp{k})
        ylim(y_lim_Cmx)
        xlim(x_lim)
        hold on
        text(1./vr_vec(5),getDataPoint(ILCFHydroModelHZobj,'Cmx', Y_vec(i), X_vec(i), theta_vec(k)/180, vr_vec(5)),theta_set_disp(k),'Color',col(k,:))
    end
    legend(theta_set_disp)
    for k = 1:length(theta_vec)    
        h = plot1D_f(ILCFHydroModelJMDobj , 'Cmx', Y_vec(i), X_vec(i), theta_vec(k), vr_vec);
        set(h,'Color',col(k,:))
        set(h,'Line','--')
        set(h,'Marker',shp{k})
    end
    hold on
    plot(0.125*[1 1],y_lim_Cmx,'--r')
    plot(0.25*[1 1],y_lim_Cmx,'--r')
    plot(basic_bare_inline_MT(:,1)./2,basic_bare_inline_MT(:,3),'s--')
    suptitle([marker ' IL Added Mass Coefficient @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(i))])
    legend1 = legend(gca,'show');
    uistack(legend1, 'top');
    %         set(legend1,'Position',[10 20 100 200]);
    saveas(gcf,[outputDir filesep 'Cmx' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) marker 'Comp.jpg'])
    saveas(gcf,[outputDir filesep 'Cmx' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) marker 'Comp.fig'])
    close

end



% plot all CLv
for i = 1:length(X_vec)
    figure
    for k = 1:length(theta_vec)
        h = plot1D_f(ILCFHydroModelHZobj , 'CLv', Y_vec(i), X_vec(i), theta_vec(k), vr_vec);
        hold on
        set(h,'Color',col(k,:))
        set(h,'Line','-')
        set(h,'Marker',shp{k})
        text(1./vr_vec(5),getDataPoint(ILCFHydroModelHZobj,'CLv', Y_vec(i), X_vec(i), theta_vec(k)/180, vr_vec(5)),theta_set_disp(k),'Color',col(k,:))
    end
    legend(theta_set_disp)
    for k = 1:length(theta_vec)
        h = plot1D_f(ILCFHydroModelJMDobj , 'CLv', Y_vec(i), X_vec(i), theta_vec(k), vr_vec);
        set(h,'Color',col(k,:))
        set(h,'Line','--')
        set(h,'Marker',shp{k})
        ylim(y_lim_CLv)
        xlim(x_lim)
        hold on
    end
    plot(0.125*[1 1],y_lim_CLv,'--r')
    plot(0.25*[1 1],y_lim_CLv,'--r')
    suptitle([marker ' CF CLv @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(i))])
    legend1 = legend(gca,'show');
    uistack(legend1, 'top');
    %         set(legend1,'Position',[10 20 100 200]);
    saveas(gcf,[outputDir filesep 'CLv' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) marker 'Comp.jpg'])
    saveas(gcf,[outputDir filesep 'CLv' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) marker 'Comp.fig'])
    close
end
% plot all CDv
for i = 1:length(X_vec)

    figure
    for k = 1:length(theta_vec)
        h = plot1D_f(ILCFHydroModelHZobj , 'CDv', Y_vec(i), X_vec(i), theta_vec(k), vr_vec);
        set(h,'Color',col(k,:))
        set(h,'Line','-')
        set(h,'Marker',shp{k})
        ylim(y_lim_CDv)
        xlim(x_lim)
        hold on
        text(1./vr_vec(5),getDataPoint(ILCFHydroModelHZobj,'CDv', Y_vec(i), X_vec(i), theta_vec(k)/180, vr_vec(5)),theta_set_disp(k),'Color',col(k,:))
    end
    legend(theta_set_disp)
    for k = 1:length(theta_vec)
        h = plot1D_f(ILCFHydroModelJMDobj , 'CDv', Y_vec(i), X_vec(i), theta_vec(k), vr_vec);
        set(h,'Color',col(k,:))
        set(h,'Line','--')
        set(h,'Marker',shp{k})
    end
    plot(0.125*[1 1],y_lim_CDv,'--r')
    plot(0.25*[1 1],y_lim_CDv,'--r')
    suptitle([marker ' IL CDv Coefficient @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(i))])
    legend1 = legend(gca,'show');
    uistack(legend1, 'top');
    %         set(legend1,'Position',[10 20 100 200]);
    saveas(gcf,[outputDir filesep 'CDv' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) marker 'Comp.jpg'])
    saveas(gcf,[outputDir filesep 'CDv' 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) marker 'Comp.fig'])
    close
end