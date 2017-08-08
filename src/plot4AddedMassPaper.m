function plot4AddedMassPaper
close all
clear all

outputDir =['/Users/haining/Dropbox/Writting/AddMass/figFiles/'];

color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
targetpar_set = {'CL1','CL3','CD2','Cmy', 'Cmx','pow','CDv'}; 


%% fig 1

bb_viva = load('basic_bare');
basic_bare_inline_MT = load('basic_bare_inline_MT');

figure
plot(bb_viva(:,1),bb_viva(:,3),'-')
xlabel('Non-dimensional Frequency')
ylabel('Added Mass Coefficient in Crossflow Direction')
saveas(gcf, [outputDir 'Fig1.fig'])    
saveas(gcf, [outputDir 'Fig1.png'])
saveas(gcf, [outputDir 'Fig1.eps'])


%% fig 2
freqRatio = 2;
energyRatio = 0.5;

ts1 = 0:0.05:1;
freq1=1;
freq2=freqRatio*freq1;
omega1 = freq1*2*pi;
omega2 = freq2*2*pi;

A1 = 0.7;
A2 = sqrt(energyRatio) * A1;

plotjpg =0;
plotfig =1; 

motion1 = A1*sin(freq1*2*pi*ts1);
motion1_cos = A1*cos(freq1*2*pi*ts1);


titleString = {'-\pi', '-3/4\pi','-1/2\pi','-1/4\pi','0','1/4\pi','1/2\pi','3/4\pi','\pi'};
figure('Position', [100,100,900,300])
for i =1:9
subplot(1,9,i)
motion2 = A2*sin(freq2*2*pi*ts1 + (i-5)*1/4*pi);
quiver(motion2(1:end-1),  motion1(1:end-1), motion2(2:end) - ...
    motion2(1:end-1), motion1(2:end) - motion1(1:end-1),'b')
xlim([-1, 1])
ylim([-1, 1])
title(titleString{i})
% set(gca,'xtick',[])
% set(gca,'ytick',[])
axis off
end
saveas(gcf, [outputDir 'Fig2.fig'])    
saveas(gcf, [outputDir 'Fig2.png'])
saveas(gcf, [outputDir 'Fig2.eps'])


% figure 3

load('AddMassPaper/ZDMITDatabase.mat')

Y_vec = [0.15, 0.25, 0.5, 0.75, 1, 1.25, 1.5 ];
X_vec = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.45,0.6,0.75 ];
% theta_vec = [-180:22.5:-45, -30:15:30, 45:22.5:180];
theta_vec =  [-180:22.5:-135, -120:15:-45, -22.5:22.5:180];
vr_vec_ori = [4:0.5:8];


% plot3D_theta(obj, targetpar, y, x, theta, vr)

figure
plot3D_theta(ILCFHydroModelZDMiTobj , 'CLv', Y_vec, X_vec, theta_vec(3), vr_vec_ori)

for k = 1:length(theta_vec)
    
    h(k) = plot1D_f(ILCFHydroModelobj , 'CLv', Y_vec(j), X_vec(i), theta_vec(k), vr_vec);
    set(h(k),'Color',col(k,:))
    set(h(k),'Line','-')
    set(h(k),'Marker',shp{k})
%             ylim(y_lim_CLv)
    xlim(x_lim)
    hold on
%             text(1./vr_vec(5),getDataPoint(ILCFHydroModelobj,'CLv', Y_vec(j), X_vec(i), theta_vec(k)/180, vr_vec(5)),theta_set_disp(k),'Color',col(k,:))
end
%         plot(0.125*[1 1],y_lim_CLv,'--r')
%         plot(0.25*[1 1],y_lim_CLv,'--r')
%         suptitle([marker ' CF CLv @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
legend(theta_set_disp)
legend1 = legend(gca,'show');
uistack(legend1, 'top');
%         set(legend1,'Position',[10 20 100 200]);
saveas(gcf, [outputDir 'Fig3.fig'])
saveas(gcf, [outputDir 'Fig3.png'])
saveas(gcf, [outputDir 'Fig3.eps'])


figure
for k = 1:length(theta_vec)
    h(k) = plot1D_f(ILCFHydroModelobj , 'CLv', Y_vec(j), X_vec(i), theta_vec(k), vr_vec);
    set(h(k),'Color',col(k,:))
    set(h(k),'Line','-')
    set(h(k),'Marker',shp{k})
%             ylim(y_lim_CLv)
    xlim(x_lim)
    hold on
%             text(1./vr_vec(5),getDataPoint(ILCFHydroModelobj,'CLv', Y_vec(j), X_vec(i), theta_vec(k)/180, vr_vec(5)),theta_set_disp(k),'Color',col(k,:))
end
%         plot(0.125*[1 1],y_lim_CLv,'--r')
%         plot(0.25*[1 1],y_lim_CLv,'--r')
%         suptitle([marker ' CF CLv @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
legend(theta_set_disp)
legend1 = legend(gca,'show');
uistack(legend1, 'top');
%         set(legend1,'Position',[10 20 100 200]);
saveas(gcf, [outputDir 'Fig4.fig'])
saveas(gcf, [outputDir 'Fig4.png'])
saveas(gcf, [outputDir 'Fig4.eps'])




