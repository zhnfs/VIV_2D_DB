% function vis_haining_masterscripts()

%% it works with haining's 2D database
% global vr_set theta_set X_set Y_set color_set targetpar_set Y_vec X_vec all_phases
clear all
utility = Utility.getInstance();
filesep = utility.getFilesep();

col = distinguishable_colors(100);

vr_set = {'4' '4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'}; 

X_set = {'0.05','0.1','0.15','0.2','0.25','0.3'};
Y_set = {'0.15','0.25','0.5','0.75', '1'};

theta_set = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
                 '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
     
            
phase_vec = [-180:22.5:-45, -30:15:30, 45:22.5:180];
X_vec = [0.05:0.05:0.3];
Y_vec = [0.15,0.25,0.5,0.75,1];
vr_vec = [4:0.5:8];

Y_vector_2 = [0.15 0.15 0.25 0.25 0.25 0.5 0.5 0.5 0.5 0.75 0.75 0.75 0.75 0.75 0.75 1 1 1 1 1 1];
X_vector_2 = [0.05 0.1 0.05 0.1 0.15 0.05 0.1 0.15 0.2 0.05 0.1 0.15 0.2 0.25 0.3 0.05 0.1 0.15 0.2 0.25 0.3 ];


color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
targetpar_set = {'CLv','CDv','Cmy','Cmx','pow'};   

% dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseJMD/';
% outputPath =  '/Users/haining/VIV/src/ILCFprediction/outputJMD/';

dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseHZ/';
outputPath =  '/Users/haining/VIV/src/ILCFprediction/outputHZ/';
utility.conditionalMkdir(outputPath);

bb_viva = load('basic_bare');

ILCFHydroModelHZobj = ILCFHydroModelHZ();

% Added mass comparion with VIVA
figure
hold on
for k = 1:length(phase_vec)-1
    h(1) = plot1D_f(ILCFHydroModelHZobj , 'Cmy', 0.3, 0.05, phase_vec(k), [4:0.5:8]);
    set(h(1),'Color',col(k,:))
    set(h(1),'Line','-')
    hold on
end
hold on
plot(bb_viva(:,1),bb_viva(:,3),'s--')
xlim([0.11 0.27])
suptitle(['Added mass CF @ Y:0.3 X:0.05'])
legend([theta_set,'viva'],'Location','BestOutside')
legend1 = legend(gca,'show');
uistack(legend1, 'top');

figure
hold on
for k = 1:length(phase_vec)-1
    h(1) = plot1D_f(ILCFHydroModelHZobj , 'Cmy', 0.5, 0.05, phase_vec(k), [4:0.5:8]);
    set(h(1),'Color',col(k,:))
    set(h(1),'Line','-')
    hold on
end
hold on
plot(bb_viva(:,1),bb_viva(:,3),'s--')
xlim([0.11 0.27])
suptitle(['Added mass CF @ Y:0.5 X:0.05'])
legend([theta_set,'viva'],'Location','BestOutside')
legend1 = legend(gca,'show');
uistack(legend1, 'top');


% plot all Clv to see if there is 2 peaks
for i = 1:length(Y_vector_2)
    Y_val = Y_vector_2(i);
    X_val = X_vector_2(i);
    figure
    for k = 1:length(phase_vec)-1
%         subplot(4,5,k)
        h(1) = plot1D_f_ori(ILCFHydroModelHZobj , 'CLv', Y_val, X_val, phase_vec(k), [4:0.5:8]);
        hold on
        set(h(1),'Color',col(k,:))
        set(h(1),'Line','-')
        ylim([0 2.5])
    end
    suptitle(['CF Lift Coefficient @ X:  ' num2str(X_val) ' Y:  ' num2str(Y_val)])
    legend(theta_set)
    legend1 = legend(gca,'show');
    uistack(legend1, 'top');
    %         set(legend1,'Position',[10 20 100 200]);
    saveas(gcf,[outputPath filesep 'CLv' ' X' num2str(X_val) 'Y:  ' num2str(Y_val) '.jpg'])
    close
end

plot1D_f_ori(ILCFHydroModelHZobj , 'CLv', 0.75, 0.15, 0, [4.5:0.5:8]);


figure
plot1D_f(ILCFHydroModelHZobj , 'Cmy', 0.25, 0.15, 90, [4.5:0.5:8]);
hold on
plot1D_f_ori(ILCFHydroModelHZobj , 'Cmy', 0.25, 0.15, 90, [4.5:0.5:8]);
legend('smooth','original')

% compare withJMD
ILCFHydroModelJMDobj = ILCFHydroModel();

figure
plot1D_f(ILCFHydroModelJMDobj , 'Cmy', 0.25, 0.15, 90, [4.5:0.5:8]);
hold on
plot1D_f(ILCFHydroModelHZobj , 'Cmy', 0.25, 0.15, 90, [4.5:0.5:8]);
legend('JMD Smooth','HZ Smooth')

figure
plot1D_f_ori(ILCFHydroModelJMDobj , 'Cmy', 0.25, 0.15, 90, [4.5:0.5:8]);
hold on
plot1D_f_ori(ILCFHydroModelHZobj , 'Cmy', 0.25, 0.15, 90, [4.5:0.5:8]);
legend('JMD Ori','HZ Ori')

figure
plot1D_f(ILCFHydroModelJMDobj , 'CLv', 0.25, 0.15, 90, [4.5:0.5:8]);
hold on
plot1D_f(ILCFHydroModelHZobj , 'CLv', 0.25, 0.15, 90, [4.5:0.5:8]);
legend('JMD Smooth','HZ Smooth')

figure
plot1D_f_ori(ILCFHydroModelJMDobj , 'CLv', 0.25, 0.15, 90, [4.5:0.5:8]);
hold on
plot1D_f_ori(ILCFHydroModelHZobj , 'CLv', 0.25, 0.15, 90, [4.5:0.5:8]);
legend('JMD Ori','HZ Ori')



% [rawDatabaseStrucHZ] = loadDatabaseHZ();
% 
% [rawDatabaseStrucJMD] = loadDatabaseJMD();
% 
% rawDatabaseStruc = combineDatabase(rawDatabaseStrucHZ, rawDatabaseStrucJMD);           
% 
% 
% %-- inpaint using INPAINTN 
% bigCmy_int = inpaintn(ILCFHydroModelHZobj.bigCmy); 
% 
% x = ILCFHydroModelHZobj.bigX(:,:,:,1);
% y = ILCFHydroModelHZobj.bigY(:,:,:,1);
% z = ILCFHydroModelHZobj.bigtheta(:,:,:,1).*180;
% Cmy_temp = ILCFHydroModelHZobj.bigCmy(:,:,:,1);
% Cmy_int = bigCmy_int(:,:,:,1);
% 
% xmin = min(x(:)); xmax = max(x(:)); 
% zmin = min(z(:)); ymax = max(y(:)); 
% 
% figure
% subplot(221), imagesc(Cmy_temp(:,:,1) ), axis equal off 
% hold on; contour(Cmy_temp(:,:,1) )
% title('Original plane, Theta -180, Vr 4.5') 
% % subplot(222), imagesc(ILCFHydroModelJMDobj.bigCmy_smooth(:,:,1,1)), axis equal off 
% subplot(222), imagesc(Cmy_int(:,:,1)), axis equal off 
% hold on; contour(Cmy_int(:,:,1) )
% title('Smoothed plane, Theta -180, Vr 4.5') 
% subplot(223) 
% hsurfaces = slice(x,y,z,Cmy_temp,[xmin,(xmin+xmax)/2,xmax],ymax,zmin); 
% set(hsurfaces,'FaceColor','interp','EdgeColor','none') 
% hcont = contourslice(x,y,z,Cmy_temp ,[xmin,(xmin+xmax)/2,xmax],ymax,zmin); 
% set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5) 
% view(3),  axis tight %,daspect([2,2,1]),
% title('Original data compared with...') 
% subplot(224) 
% hsurfaces = slice(x,y,z,Cmy_int,[xmin,(xmin+xmax)/2,xmax],ymax,zmin); 
% set(hsurfaces,'FaceColor','interp','EdgeColor','none') 
% hcont = contourslice(x,y,z,Cmy_int,[xmin,(xmin+xmax)/2,xmax],ymax,zmin); 
% set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5) 
% view(3),  axis tight %,daspect([2,2,1]),
% title('... reconstructed data')


phase_arosen = [77 84 264]-180;
figure
subplot(4,1,1)
for i =1:3
    h(i) = plot1D_f(ILCFHydroModelHZobj , 'CLv', 0.3, 0.15, phase_arosen(i), [4:0.5:10]);
    hold on
end
set(h(1),'Color','g')
set(h(1),'Marker','o')
set(h(2),'Marker','*')
set(h(3),'Marker','p')
legend('77','84','264')
subplot(4,1,2)
for i =1:3
    h(i) = plot1D_f(ILCFHydroModelHZobj , 'Cmy', 0.3, 0.15, phase_arosen(i), [4:0.5:10]);
    hold on
end
set(h(1),'Color','g')
set(h(1),'Marker','o')
set(h(2),'Marker','*')
set(h(3),'Marker','p')
legend('77','84','264')
subplot(4,1,3)
for i =1:3
    h(i) = plot1D_f(ILCFHydroModelHZobj , 'CDv', 0.3, 0.15, phase_arosen(i), [4:0.5:10]);
    hold on
end
set(h(1),'Color','g')
set(h(1),'Marker','o')
set(h(2),'Marker','*')
set(h(3),'Marker','p')
legend('77','84','264')
subplot(4,1,4)
for i =1:3
    h(i) = plot1D_f(ILCFHydroModelHZobj , 'Cmx', 0.3, 0.15, phase_arosen(i), [4:0.5:10]);
    hold on
end
set(h(1),'Color','g')
set(h(1),'Marker','o')
set(h(2),'Marker','*')
set(h(3),'Marker','p')
legend('77','84','264')



plot2D_vr(ILCFHydroModelHZobj , 'Cmy', 0.25, 0.15, phase_vec(1:9), [3:0.5:12]);
figure
for i = 1:length(phase_vec)-1
    subplot(2,4,i)
    plot1D_vr(ILCFHydroModelHZobj , 'Cmy', 0.25, 0.15, phase_vec(i), [3:0.2:15]);
    ylim([-2.5 5])
end

plot2D_vr(ILCFHydroModelHZobj , 'CLv', 0.25, 0.15, phase_vec(1:9), [3:0.5:12]);
figure
for i = 1:length(phase_vec)-1
    subplot(2,4,i)
    plot1D_vr(ILCFHydroModelHZobj , 'CLv', 0.25, 0.15, phase_vec(i), [3:0.2:15]);
%     ylim([-2.5 5])
end

plot2D_vr(ILCFHydroModelHZobj , 'CLv', 0.75, 0.25, phase_vec(1:9), [3:0.5:12]);
figure
for i = 1:length(phase_vec)-1
    subplot(2,4,i)
    plot1D_vr(ILCFHydroModelHZobj , 'CLv', 0.75, 0.25, phase_vec(i), [3:0.2:15]);
%     ylim([-2.5 5])
end



plot1D_vr(ILCFHydroModelHZobj , 'Cmy', 0.25, 0.15, 90, [3:0.5:12]);
plot3D_theta(ILCFHydroModelHZobj , 'Cmy', [0.25:0.25:1.5], [0:0.15:0.75], phase_vec(2), [3:0.5:12]);


plot_para_2D(dataPath, outputPath, 'X',[2],'Y',[1],'Cmy',1, vr_set, theta_set, X_set, Y_set, Y_vec);
plot_para(dataPath, outputPath, 'theta',[2],'Cmy',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);

% plot1D_vr(ILCFHydroModelobj , 'pow', 1, 0.5, 90/180, [6:0.2:8]);
% plot1D_vr(ILCFHydroModelobj , 'Cmy', 1, 0.5, 90/180, [6:0.2:8]);
% plot1D_theta(ILCFHydroModelobj , 'Cmy', 1, 0.5, [-180:45:180]./180, 6);
% plot1D_y(ILCFHydroModelobj , 'Cmy', [0:0.1:1.5], 0.5, 90/180, 6);
% plot1D_x(ILCFHydroModelobj , 'Cmy', 1, [0:0.05:0.75], 90/180, 6);
% 

%  
% % pca_para(dataPath, outputPath, vr_set, theta_set, X_set, Y_set, Y_vec);
% 
% 3D
plot_para(dataPath, outputPath, 'theta',[2],'CLv',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);
plot_para(dataPath, outputPath, 'theta',[2],'Cmy',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);

plot_para(dataPath, outputPath, 'theta',[2],'Cmx',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);
plot_para(dataPath, outputPath, 'theta',[2],'Pow',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);

% % 2D
plot_para_2D(dataPath, outputPath, 'X',[1],'Y',[2],'Cmy',1, vr_set, theta_set, X_set, Y_set, Y_vec);
plot_para_2D(dataPath, outputPath, 'X',[1],'Y',[2],'Cmx',1, vr_set, theta_set, X_set, Y_set, Y_vec);

plot_para_2D(dataPath, outputPath, 'vr',[4],'theta',[5],'Cmy',1, vr_set, theta_set, X_set, Y_set, Y_vec);
plot_para_2D(dataPath, outputPath, 'vr',[4],'theta',[5],'Cmx',1, vr_set, theta_set, X_set, Y_set, Y_vec);

plot_para_2D(dataPath, outputPath, 'X',[1],'Y',[2],'CLv',1, vr_set, theta_set, X_set, Y_set, Y_vec);
plot_para_2D(dataPath, outputPath, 'X',[2],'Y',[2],'CLv',1, vr_set, theta_set, X_set, Y_set, Y_vec);

% % 1D
% plot_para_1D(dataPath, outputPath, 'vr',[4],'theta',[5], 'X',[1],'Cmy', vr_set, theta_set, X_set, Y_set, Y_vec);


% vrIndex = [4];
% thetaIndex = [1:8];
% AmpXIndex = [1:6];
% AmpYIndex = [1:6];
% 
% thetaIndex = [1];
% targetpar = 'Clv';
% 
% [value, vr, theta, X,Y] = extract_para_2D(dataPath, vrIndex, thetaIndex, AmpXIndex, AmpYIndex, targetpar, vr_set, theta_set, X_set, Y_set, Y_vec)


% %% Jason's cases
% for i = 1:length(vr_set)
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'CL1',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'Cmx',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'Cmy',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'CL3',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'pow',1, expName);
% end

%% all cases
% for i = 1:length(vr_set)
%     
%     plot_para(dataPath, outputPath, 'vr',[i],'CL1',2, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:6],'CL1',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[7:12],'CL1',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[13:18],'CL1',1, expName);
%     
%     plot_para(dataPath, outputPath, 'vr',[i],'Cmy',2, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:6],'Cmy',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[7:12],'Cmy',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[13:18],'Cmy',1, expName);
% 
%     plot_para(dataPath, outputPath, 'vr',[i],'Cmx',2, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:6],'Cmx',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[7:12],'Cmx',1, expName);    
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[13:18],'Cmx',1, expName);
%     
%     plot_para(dataPath, outputPath, 'vr',[i],'CL3',2, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:6],'CL3',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[7:12],'CL3',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[13:18],'CL3',1, expName);
%     
%     plot_para(dataPath, outputPath, 'vr',[i],'pow',2, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:6],'pow',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[7:12],'pow',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[13:18],'pow',1, expName);
% 
% end

% vrNumber = [4];
% phaseNumber = [1:8];
% for i=1
%     plot_para(dataPath, outputPath, 'vr',[i],'CL1',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,Y_vec);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'CL1',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, Y_vec);
%     plot_para(dataPath, outputPath, 'vr',[i],'Cmy',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,Y_vec);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'Cmy',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, Y_vec);
%     plot_para(dataPath, outputPath, 'vr',[i],'Cmx',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,Y_vec);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'Cmx',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, Y_vec);
%     plot_para(dataPath, outputPath, 'vr',[i],'CL3',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,Y_vec);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'CL3',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, Y_vec);
%     plot_para(dataPath, outputPath, 'vr',[i],'pow',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,Y_vec);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:8],'pow',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, Y_vec);
% end


% close all
