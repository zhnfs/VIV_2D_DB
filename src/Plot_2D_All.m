function Plot_2D_All(ILCFHydroModelZDMiTobj,ILCFHydroModelHZobj,ILCFHydroModelJMDobj,theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
   outputDir, marker)
     
Y_vec = [0.15, 0.25, 0.5, 0.75, 1, 1.25, 1.5 ];
X_vec = [0.05,0.1,0.15,0.2,0.25,0.3,0.45,0.6,0.75 ];

targetVar = 'Cmy';
figure
for i = 1:length(Y_vec)
    for j = 1:length(X_vec)        
    y_val = Y_vec(i);
    x_val = X_vec(j);
    subplot(length(Y_vec),length(X_vec),length(X_vec)*(i-1)+j,'CLim',y_lim_Cmy)
    plot2D_theta_vr(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);xlim(x_lim)
%     title(['Y: ' num2str(y_val) ' X: ' num2str(x_val) ])
%     set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
end
h = colorbar;
set(h,'Position',[0.92 0.1 0.03 0.8]);
suptitle(targetVar)   
saveas(gcf,[outputDir filesep targetVar marker 'All2DColorbar.jpg'])
saveas(gcf,[outputDir filesep targetVar marker 'All2DColorbar.fig'])
close

targetVar = 'Cmx';
figure
for i = 1:length(Y_vec)
    for j = 1:length(X_vec)        
    y_val = Y_vec(i);
    x_val = X_vec(j);
    subplot(length(Y_vec),length(X_vec),length(X_vec)*(i-1)+j,'CLim',y_lim_Cmx)
    plot2D_theta_vr(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);xlim(x_lim)
%     title(['Y: ' num2str(y_val) ' X: ' num2str(x_val) ])
%     set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
end
h = colorbar;
set(h,'Position',[0.92 0.1 0.03 0.8]);
suptitle(targetVar)   
saveas(gcf,[outputDir filesep targetVar marker 'All2DColorbar.jpg'])
saveas(gcf,[outputDir filesep targetVar marker 'All2DColorbar.fig'])
close

targetVar = 'CLv';
figure
for i = 1:length(Y_vec)
    for j = 1:length(X_vec)        
    y_val = Y_vec(i);
    x_val = X_vec(j);
    subplot(length(Y_vec),length(X_vec),length(X_vec)*(i-1)+j,'CLim',y_lim_CLv)
    plot2D_theta_vr(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);xlim(x_lim)
%     title(['Y: ' num2str(y_val) ' X: ' num2str(x_val) ])
%     set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
end
h = colorbar;
set(h,'Position',[0.92 0.1 0.03 0.8]);
suptitle(targetVar)   
saveas(gcf,[outputDir filesep targetVar marker 'All2DColorbar.jpg'])
saveas(gcf,[outputDir filesep targetVar marker 'All2DColorbar.fig'])
close

targetVar = 'CDv';
figure
for i = 1:length(Y_vec)
    for j = 1:length(X_vec)        
    y_val = Y_vec(i);
    x_val = X_vec(j);
    subplot(length(Y_vec),length(X_vec),length(X_vec)*(i-1)+j,'CLim',y_lim_CDv)
    plot2D_theta_vr(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);xlim(x_lim)
%     title(['Y: ' num2str(y_val) ' X: ' num2str(x_val) ])
%     set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
    set(gca,'ytick',[]);set(gca,'xtick',[]);
    end
end
h = colorbar;
set(h,'Position',[0.92 0.1 0.03 0.8]);
suptitle(targetVar)   
saveas(gcf,[outputDir filesep targetVar marker 'All2DColorbar.jpg'])
saveas(gcf,[outputDir filesep targetVar marker 'All2DColorbar.fig'])
close
