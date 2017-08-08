function plotDatabaseCompAll(ILCFHydroModelZDMiTobj,ILCFHydroModelHZobj,ILCFHydroModelJMDobj,targetVar, ...
    Y_vec, X_vec, theta_vec, vr_vec,theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)

if strcmp(targetVar,'CLv')
    ylim_comp = [-1 2.5];
elseif  strcmp(targetVar,'CDv')
    ylim_comp = [-1 2.5];
elseif  strcmp(targetVar,'Cmy')
    ylim_comp = [-2 4];
elseif  strcmp(targetVar,'Cmx')
    ylim_comp = [-2 4];
end

% ylim_comp = [-5 5];

for i = 1:length(X_vec)   
    y_val = Y_vec(i);
    x_val = X_vec(i);
    for k = 1:length(theta_vec)
        figure
        subplot(221),
        h(1) = plot1D_f_ori(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec(k), vr_vec);hold on;
        h(2) = plot1D_f_int(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec(k), vr_vec);
        h(3) = plot1D_f(ILCFHydroModelZDMiTobj     , targetVar, y_val, x_val, theta_vec(k), vr_vec);
        plot(0.125*[1 1],[-10 10],'--r')
        plot(0.25*[1 1],[-10 10],'--r')
        ylim(ylim_comp)
        set(h(1),'Color','g');set(h(1),'Marker','s');
        set(h(2),'Color','r');set(h(2),'Marker','*');
        set(h(3),'Color','b');set(h(3),'Marker','d');
        legend('original', 'inpaint', 'smooth','Location','SouthEast')
        title(['Cmy @ ' 'x ' num2str(x_val) ' y ' num2str(y_val) 'ph ' num2str(theta_vec(k))])

        subplot(222),
        h(1) = plot1D_f_ori(ILCFHydroModelZDMiTobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        hold on
        h(2) = plot1D_f_ori(ILCFHydroModelHZobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        h(3) = plot1D_f_ori(ILCFHydroModelJMDobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        plot(0.125*[1 1],[-10 10],'--r')
        plot(0.25*[1 1],[-10 10],'--r')
        ylim(ylim_comp)
        set(h(1),'Color','g');set(h(1),'Marker','s');
        set(h(2),'Color','r');set(h(2),'Marker','*');
        set(h(3),'Color','b');set(h(3),'Marker','d');
        title('Original')
        legend('Combo', 'HZ', 'JMD','Location','SouthEast')

        subplot(223),
        h(1) = plot1D_f(ILCFHydroModelZDMiTobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        hold on
        h(2) = plot1D_f(ILCFHydroModelHZobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        h(3) = plot1D_f(ILCFHydroModelJMDobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        plot(0.125*[1 1],[-10 10],'--r')
        plot(0.25*[1 1],[-10 10],'--r')
        ylim(ylim_comp)
        set(h(1),'Color','g');set(h(1),'Marker','s');
        set(h(2),'Color','r');set(h(2),'Marker','*');
        set(h(3),'Color','b');set(h(3),'Marker','d');
        legend('Combo', 'HZ', 'JMD','Location','SouthEast')
        title('smooth')

        subplot(224)
        h(1) = plot1D_f_ori(ILCFHydroModelHZobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        hold on
        h(2) = plot1D_f(ILCFHydroModelHZobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        h(3) = plot1D_f_ori(ILCFHydroModelJMDobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        h(4) = plot1D_f(ILCFHydroModelJMDobj, targetVar, y_val, x_val, theta_vec(k), vr_vec);
        plot(0.125*[1 1],[-10 10],'--r')
        plot(0.25*[1 1],[-10 10],'--r')
        ylim(ylim_comp)
        set(h(1),'Color','r');set(h(1),'Marker','s');
        set(h(2),'Color','m');set(h(2),'Marker','*');
        set(h(3),'Color','b');set(h(3),'Marker','d');
        set(h(4),'Color','g');set(h(1),'Marker','o'); % 'v','x'
        legend('HZ Original', 'HZ Smooth', 'JMD Original', 'JMD Smooth','Location','SouthEast')

        saveas(gcf,[outputDir filesep targetVar 'X' num2str(X_vec(i)) 'Y' num2str(Y_vec(i)) 'ph ' num2str(theta_vec(k)) marker 'CompAll.jpg'])
        close
    end

end