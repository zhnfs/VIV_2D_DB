function Plot_Remi_2D_compori(ILCFHydroModelZDMiTobj,ILCFHydroModelHZobj,ILCFHydroModelJMDobj,theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)
            
y_val = 0.5;
x_val = 0.05;

load('proc_Haining_052213.mat')

a_remi = 0:45:315;
a_real = -90:45:225;
a_rot = -180:45:180;

clv_tot_rot =[ clv_tot(7:8,:); clv_tot(1:6,:); clv_tot(7,:)];
cdv_tot_rot =[ cdv_tot(7:8,:); cdv_tot(1:6,:); cdv_tot(7,:)];
cmy_tot_rot =[ cmy_tot(7:8,:); cmy_tot(1:6,:); cmy_tot(7,:)];
cmx_tot_rot =[ cmx_tot(7:8,:); cmx_tot(1:6,:); cmx_tot(7,:)];

targetVar = 'Cmy';
figure
subplot(221)
plot2D_theta_vr_ori(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank Test All Re 7600')
subplot(222)
plot2D_theta_vr_ori(ILCFHydroModelHZobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank 2011 Re 7600')
subplot(223)
plot2D_theta_vr_ori(ILCFHydroModelJMDobj , targetVar, y_val, 0.15, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank 2006 Re 7600')
subplot(224)
surface(Ur,a_rot,cmy_tot_rot);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
% xlim([vr_vec_disp(1) vr_vec_disp(end)])
title('DNS 2012 Re 100')
suptitle([targetVar ' Original'])   
saveas(gcf,[outputDir filesep targetVar 'X' num2str(x_val) 'Y' num2str(y_val) marker 'RemiCompOri.jpg'])
saveas(gcf,[outputDir filesep targetVar 'X' num2str(x_val) 'Y' num2str(y_val) marker 'RemiCompOri.fig'])
close

targetVar = 'Cmx';
figure
subplot(221)
plot2D_theta_vr_ori(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank Test All Re 7600')
subplot(222)
plot2D_theta_vr_ori(ILCFHydroModelHZobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank 2011 Re 7600')
subplot(223)
plot2D_theta_vr_ori(ILCFHydroModelJMDobj , targetVar, y_val, 0.15, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank 2006 Re 7600')
subplot(224)
surface(Ur,a_rot,cmx_tot_rot);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
% xlim([vr_vec_disp(1) vr_vec_disp(end)])
title('DNS 2012 Re 100')
suptitle([targetVar ' Original'])   
saveas(gcf,[outputDir filesep targetVar 'X' num2str(x_val) 'Y' num2str(y_val) marker 'RemiCompOri.jpg'])
saveas(gcf,[outputDir filesep targetVar 'X' num2str(x_val) 'Y' num2str(y_val) marker 'RemiCompOri.fig'])
close


targetVar = 'CLv';
figure
subplot(221)
plot2D_theta_vr_ori(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank Test All Re 7600')
subplot(222)
plot2D_theta_vr_ori(ILCFHydroModelHZobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank 2011 Re 7600')
subplot(223)
plot2D_theta_vr_ori(ILCFHydroModelJMDobj , targetVar, y_val, 0.15, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank 2006 Re 7600')
subplot(224)
surface(Ur,a_rot,clv_tot_rot);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
% xlim([vr_vec_disp(1) vr_vec_disp(end)])
title('DNS 2012 Re 100')
suptitle([targetVar ' Original'])   
saveas(gcf,[outputDir filesep targetVar 'X' num2str(x_val) 'Y' num2str(y_val) marker 'RemiCompOri.jpg'])
saveas(gcf,[outputDir filesep targetVar 'X' num2str(x_val) 'Y' num2str(y_val) marker 'RemiCompOri.fig'])
close

targetVar = 'CDv';
figure
subplot(221)
plot2D_theta_vr_ori(ILCFHydroModelZDMiTobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank Test All Re 7600')
subplot(222)
plot2D_theta_vr_ori(ILCFHydroModelHZobj , targetVar, y_val, x_val, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank 2011 Re 7600')
subplot(223)
plot2D_theta_vr_ori(ILCFHydroModelJMDobj , targetVar, y_val, 0.15, theta_vec_disp, vr_vec_disp);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
title('Towtank 2006 Re 7600')
subplot(224)
surface(Ur,a_rot,cdv_tot_rot);colorbar; xlim(x_lim);set(gca,'ytick',[-180 -90 0 90 180]);set(gca,'xtick',[4 6 8]);
% xlim([vr_vec_disp(1) vr_vec_disp(end)])
title('DNS 2012 Re 100')
suptitle([targetVar ' Original'])   
saveas(gcf,[outputDir filesep targetVar 'X' num2str(x_val) 'Y' num2str(y_val) marker 'RemiCompOri.jpg'])
saveas(gcf,[outputDir filesep targetVar 'X' num2str(x_val) 'Y' num2str(y_val) marker 'RemiCompOri.fig'])
close