% function vis_jason_masterscripts()

%% it works with haining's 2D database
% global vr_set theta_set X_set Y_set color_set targetpar_set Y_vec X_vec all_phases
clear all
col = distinguishable_colors(100);

utility = Utility.getInstance();
filesep = utility.getFilesep();

vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'}; 

X_set = {'0','0.15','0.3','0.45','0.6','0.75'};
Y_set = {'0.25','0.5','0.75','1', '1.25', '1.5'};

theta_set = {'180n' '135n' '90n' '45n' '0' '45' '90' '135' '180'};

phase_vec = [-180 -135 -90 -45 0 45 90 135 180];
X_vec = [0:0.15:0.75];
Y_vec = [0.25:0.25:1.5];
vr_vec = [4.5:0.5:8];

color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
targetpar_set = {'CLv','CDv','Cmy','Cmx','pow'};   

dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseJMD/';
outputPath =  '/Users/haining/VIV/src/ILCFprediction/outputJMD/';
utility.conditionalMkdir(outputPath);
%  
bb_viva = load('basic_bare');

ILCFHydroModelobj = ILCFHydroModel();


for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        figure
        for k = 1:length(phase_vec)-1
            h(1) = plot1D_f_ori(ILCFHydroModelobj , 'CLv', Y_vec(j), X_vec(i), phase_vec(k), [4.5:0.5:8]);
            set(h(1),'Color',col(k,:))
            set(h(1),'Line','-')
            ylim([0 2.5])
            hold on
        end
        suptitle(['CF Lift Coefficient @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
        legend(theta_set(1:end-1))
        legend1 = legend(gca,'show');
        uistack(legend1, 'top');
        %         set(legend1,'Position',[10 20 100 200]);
        saveas(gcf,[outputPath filesep 'CLv' ' X' num2str(X_vec(i)) 'Y:  ' num2str(Y_vec(j)) '.jpg'])
        close
    end
end

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


figure
h1 = plot2D_vr(ILCFHydroModelobj , 'Cmy', 0.5, 0, phase_vec(1:9), [4.5:0.5:8]);
hold on
h2 = plot2D_vr_ori(ILCFHydroModelobj , 'Cmy', 0.5, 0, phase_vec(1:9), [4.5:0.5:8]);
alpha(h2, 0.5)

% figure
plot2D_vr_comp(ILCFHydroModelobj , 'Cmy', 0.5, 0, phase_vec(1:9), [4.5:0.5:8]);

figure
for k = 1:length(phase_vec)-1     
    h(1) = plot1D_y_ori(ILCFHydroModelobj , 'Cmy', Y_vec, 0, phase_vec(k), 1/0.13);
    set(h(1),'Color',col(k,:))
    set(h(1),'Line','-')
    hold on
end
suptitle(['Cmy at Vr ' num2str(1/0.13,2) ' x:0'])
legend(theta_set(1:8))


figure
for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        subplot(6,6,(i-1)*6+j)
        h1 = plot1D_theta(ILCFHydroModelobj , 'Cmy', Y_vec(j), X_vec(i), phase_vec(1:9), 1/0.13);
        hold on
        h2 = plot1D_theta_ori(ILCFHydroModelobj , 'Cmy', Y_vec(j),X_vec(i), phase_vec(1:9), 1/0.13);
        title(['X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
    end
end
legend('Smooth','Original')


phase_arosen = [77 84 264]-180;
figure
subplot(4,1,1)
for i =1:3
    h(i) = plot1D_f(ILCFHydroModelobj , 'CLv', 0.3, 0.15, phase_arosen(i), [4:0.5:10]);
    hold on
end
set(h(1),'Color','g')
set(h(1),'Marker','o')
set(h(2),'Marker','*')
set(h(3),'Marker','p')
legend('77','84','264')
subplot(4,1,2)
for i =1:3
    h(i) = plot1D_f(ILCFHydroModelobj , 'Cmy', 0.3, 0.15, phase_arosen(i), [4:0.5:10]);
    hold on
end
set(h(1),'Color','g')
set(h(1),'Marker','o')
set(h(2),'Marker','*')
set(h(3),'Marker','p')
legend('77','84','264')
subplot(4,1,3)
for i =1:3
    h(i) = plot1D_f(ILCFHydroModelobj , 'CDv', 0.3, 0.15, phase_arosen(i), [4:0.5:10]);
    hold on
end
set(h(1),'Color','g')
set(h(1),'Marker','o')
set(h(2),'Marker','*')
set(h(3),'Marker','p')
legend('77','84','264')
subplot(4,1,4)
for i =1:3
    h(i) = plot1D_f(ILCFHydroModelobj , 'Cmx', 0.3, 0.15, phase_arosen(i), [4:0.5:10]);
    hold on
end
set(h(1),'Color','g')
set(h(1),'Marker','o')
set(h(2),'Marker','*')
set(h(3),'Marker','p')
legend('77','84','264')


figure
hold on
for k = 1:length(phase_vec)-1
%     subplot(2,4,k)
    h(1) = plot1D_f(ILCFHydroModelobj , 'Cmy', 0.7, 0.1, phase_vec(k), [4.5:0.5:8]);
    set(h(1),'Color',col(k,:))
    set(h(1),'Line','-')
    hold on
end
hold on
plot(bb_viva(:,1),bb_viva(:,3),'s--')
xlim([0.12 0.24])
suptitle(['Added mass CF @ Y:0.7 X:0.1'])
legend(theta_set(1:8),'viva')
legend1 = legend(gca,'show');
uistack(legend1, 'top');

figure
hold on
for k = 1:length(phase_vec)-1
%     subplot(2,4,k)
    h(1) = plot1D_f(ILCFHydroModelobj , 'Cmy', 0.3, 0.15, phase_vec(k), [4.5:0.5:8]);
    set(h(1),'Color',col(k,:))
    set(h(1),'Line','-')
    hold on
end
hold on
plot(bb_viva(:,1),bb_viva(:,3),'s--')
xlim([0.12 0.24])
suptitle(['Added mass CF @ Y:0.3 X:0.15'])
legend(theta_set(1:8),'viva')
legend1 = legend(gca,'show');
uistack(legend1, 'top');

figure
hold on
for k = 1:length(phase_vec)-1
%     subplot(2,4,k)
    h(1) = plot1D_f(ILCFHydroModelobj , 'Cmx', 0.3, 0.15, phase_vec(k), [4.5:0.5:8]);
    set(h(1),'Color',col(k,:))
%     set(h(1),'Line','--')
    hold on
end
suptitle(['Added mass IL'])
legend(theta_set)
legend1 = legend(gca,'show');
uistack(legend1, 'top');



for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        figure
        for k = 1:length(phase_vec)-1
            subplot(2,4,k)
            h(1) = plot1D_vr(ILCFHydroModelobj , 'Cmy', Y_vec(j), X_vec(i), phase_vec(k), [4.5:0.5:8]);
            set(h(1),'Color','r')
            set(h(1),'Line','--')
            hold on
            h(2) = plot1D_vr(ILCFHydroModelobj , 'Cmx', Y_vec(j), X_vec(i), phase_vec(k), [4.5:0.5:8]);
            ylim([-2.5 5])
        end
        suptitle(['Added mass @ X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
        legend('Cmy', 'cmx')
        legend1 = legend(gca,'show');
        uistack(legend1, 'top');
        %         set(legend1,'Position',[10 20 100 200]);
        saveas(gcf,[outputPath filesep 'Cm' ' X' num2str(X_vec(i)) 'Y:  ' num2str(Y_vec(j)) '.jpg'])
        close
    end
end

% X_vec = [0:0.15:0.75];
% Y_vec = [0.25:0.25:1.5];

figure
for k = 1:length(phase_vec)-1     
    h(1) = plot1D_y(ILCFHydroModelobj , 'Cmy', Y_vec, 0.1, phase_vec(k), 1/0.13);
    set(h(1),'Color',col(k,:))
    set(h(1),'Line','-')
    hold on
end
suptitle(['Cmy at Vr ' num2str(1/0.13,2) 'x:0.1'])
legend(theta_set(1:8))

figure
for k = 1:length(phase_vec)-1     
    h(1) = plot1D_y(ILCFHydroModelobj , 'Cmy', Y_vec, 0, phase_vec(k), 1/0.13);
    set(h(1),'Color',col(k,:))
    set(h(1),'Line','-')
    hold on
end
suptitle(['Cmy at Vr ' num2str(1/0.13,2) ' x:0'])
legend(theta_set(1:8))



figure
for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        subplot(length(X_vec), length(Y_vec), (i-1)*length(Y_vec)+j)
        plot1D_theta(ILCFHydroModelobj , 'Cmy', Y_vec(j), X_vec(i), [-180:45:180], 1/0.13);
        title(['X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
    end
end
suptitle('Cmy at Vr 1/0.13')

figure
for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        subplot(length(X_vec), length(Y_vec), (i-1)*length(Y_vec)+j)
        plot1D_theta(ILCFHydroModelobj , 'Cmy', Y_vec(j), X_vec(i), [-180:45:180], 4.5);
        title(['X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
    end
end
suptitle('Cmy at Vr 4.5')

figure
for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        subplot(length(X_vec), length(Y_vec), (i-1)*length(Y_vec)+j)
        plot1D_theta(ILCFHydroModelobj , 'Cmx', Y_vec(j), X_vec(i), [-180:45:180], 8);
        title(['X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
    end
end
suptitle('Cmx at Vr 8')

figure
for i = 1:length(X_vec)
    for j = 1:length(Y_vec)
        subplot(length(X_vec), length(Y_vec), (i-1)*length(Y_vec)+j)
        plot1D_theta(ILCFHydroModelobj , 'Cmx', Y_vec(j), X_vec(i), [-180:45:180], 4.5);
        title(['X:  ' num2str(X_vec(i)) ' Y:  ' num2str(Y_vec(j))])
    end
end
suptitle('Cmx at Vr 4.5')


plot2D_vr(ILCFHydroModelobj , 'Cmy', 0.25, 0.15, phase_vec(1:9), [3:0.5:12]);
figure
for i = 1:length(phase_vec)-1
    subplot(2,4,i)
    plot1D_vr(ILCFHydroModelobj , 'Cmy', 0.25, 0.15, phase_vec(i), [3:0.2:15]);
    ylim([-2.5 5])
end

plot2D_vr(ILCFHydroModelobj , 'CLv', 0.25, 0.15, phase_vec(1:9), [3:0.5:12]);
figure
for i = 1:length(phase_vec)-1
    subplot(2,4,i)
    plot1D_vr(ILCFHydroModelobj , 'CLv', 0.25, 0.15, phase_vec(i), [3:0.2:15]);
%     ylim([-2.5 5])
end

plot2D_vr(ILCFHydroModelobj , 'CLv', 0.75, 0.25, phase_vec(1:9), [3:0.5:12]);
figure
for i = 1:length(phase_vec)-1
    subplot(2,4,i)
    plot1D_vr(ILCFHydroModelobj , 'CLv', 0.75, 0.25, phase_vec(i), [3:0.2:15]);
%     ylim([-2.5 5])
end


plot1D_vr(ILCFHydroModelobj , 'Cmy', 0.25, 0.5, 90, [3:0.5:12]);
plot3D_theta(ILCFHydroModelobj , 'Cmy', [0.25:0.25:1.5], [0:0.15:0.75], phase_vec(2), [3:0.5:12]);


plot_para_2D(dataPath, outputPath, 'X',[2],'Y',[1],'Cmy',1, vr_set, theta_set, X_set, Y_set, Y_vec);
plot_para(dataPath, outputPath, 'theta',[2],'Cmy',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);

% plot1D_vr(ILCFHydroModelobj , 'pow', 1, 0.5, 90/180, [6:0.2:8]);
% plot1D_vr(ILCFHydroModelobj , 'Cmy', 1, 0.5, 90/180, [6:0.2:8]);
% plot1D_theta(ILCFHydroModelobj , 'Cmy', 1, 0.5, [-180:45:180]./180, 6);
% plot1D_y(ILCFHydroModelobj , 'Cmy', [0:0.1:1.5], 0.5, 90/180, 6);
% plot1D_x(ILCFHydroModelobj , 'Cmy', 1, [0:0.05:0.75], 90/180, 6);
% 



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
% % pca_para(dataPath, outputPath, vr_set, theta_set, X_set, Y_set, Y_vec);
% 
% 3D
plot_para(dataPath, outputPath, 'theta',[2],'CLv',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);
plot_para(dataPath, outputPath, 'theta',[2],'Cmy',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);

plot_para(dataPath, outputPath, 'theta',[2],'Cmx',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);
plot_para(dataPath, outputPath, 'theta',[2],'Pow',2, vr_set, theta_set, X_set, Y_set, color_set,Y_vec);

% % 2D
plot_para_2D(dataPath, outputPath, 'X',[1],'Y',[3],'Cmy',1, vr_set, theta_set, X_set, Y_set, Y_vec);
plot_para_2D(dataPath, outputPath, 'X',[1],'Y',[2],'Cmx',1, vr_set, theta_set, X_set, Y_set, Y_vec);

plot_para_2D(dataPath, outputPath, 'vr',[4],'theta',[5],'Cmy',1, vr_set, theta_set, X_set, Y_set, Y_vec);
plot_para_2D(dataPath, outputPath, 'vr',[4],'theta',[5],'Cmx',1, vr_set, theta_set, X_set, Y_set, Y_vec);

plot_para_2D(dataPath, outputPath, 'X',[1],'Y',[2],'CLv',1, vr_set, theta_set, X_set, Y_set, Y_vec);
plot_para_2D(dataPath, outputPath, 'X',[2],'Y',[2],'CLv',1, vr_set, theta_set, X_set, Y_set, Y_vec);

% % 1D
plot_para_1D(dataPath, outputPath, 'vr',[4],'theta',[5], 'X',[1],'Cmy', vr_set, theta_set, X_set, Y_set, Y_vec);

plot_para_1D(dataPath, outputPath, 'X',[2],'Y',[1], 'theta',[1],'Cmy', vr_set, theta_set, X_set, Y_set, Y_vec);


% close all
