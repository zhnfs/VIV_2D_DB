function vis_smalltank_masterscripts()

%% it works with haining's 2D database
% global vr_set theta_set X_set Y_set color_set targetpar_set all_Y all_X all_phases

utility = Utility.getInstance();
filesep = utility.getFilesep();

vr_set = {'4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8'}; 

X_set = {'0.05','0.1','0.15','0.2','0.25','0.3'};
Y_set = {'0.15','0.25','0.5','0.75', '1'};

all_X = [0.05,0.1,0.15,0.2, 0.25, 0.3];
all_Y = [0.15,0.25,0.5,0.75, 1 ];

Y_vector = [0.15 0.15 0.25 0.25 0.25 0.5 0.5 0.5 0.5 0.75 0.75 0.75 0.75 0.75 0.75 1 1 1 1 1 1];
X_vector = [0.05 0.1 0.05 0.1 0.15 0.05 0.1 0.15 0.2 0.05 0.1 0.15 0.2 0.25 0.3 0.05 0.1 0.15 0.2 0.25 0.3 ];

color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
targetpar_set = {'CL1','CL3','CD2','Cmy', 'Cmx','pow','CDv'};   

% dataPath =  '/Users/haining/VIV/src/2012SpringSmalltank/RepeatTest/RepeatTest6/database/';
% outputPath =  '/Users/haining/VIV/src/2012SpringSmalltank/output/RepeatTest/RepeatTest6/';

dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseHZ/';
outputPath =  '/Users/haining/VIV/src/ILCFprediction/outputHZ/';

utility.conditionalMkdir(outputPath);
    
plot_testMatrix(X_vector, Y_vector, outputPath);

theta_set = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
     '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
phases = [-180 -157.5 -135 -112.5 -90 -67.5 -45 -30 -15 0 15 30 45 67.5 90 112.5 135 157.5];

plot_para_2D(dataPath, outputPath, 'X',[1],'Y',[3],'CLv',1, vr_set, theta_set, X_set, Y_set, all_Y);
% plot_para(dataPath, outputPath, 'X',[1],'CLv',2, vr_set, theta_set, X_set, Y_set, color_set,all_Y);


% for j = [5]
%     vrNumber = [j];
%     phaseNumber = [10:12];
%     for i=1
%         plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:3],'CL1',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
%         plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:3],'Cmy',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
%         plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:3],'Cmx',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
%         plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:3],'CL3',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
%         plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1:3],'pow',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
% 
%     %     plot_para(dataPath, outputPath, 'vr',[i],'CL1',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     %     plot_para(dataPath, outputPath, 'vr',[i],'Cmy',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     %     plot_para(dataPath, outputPath, 'vr',[i],'Cmx',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     %     plot_para(dataPath, outputPath, 'vr',[i],'CL3',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     %     plot_para(dataPath, outputPath, 'vr',[i],'pow',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     end
% 
% end

% %% Jason's cases
% vrNumber = [1];
% phaseNumber = [1:18];
% for i = 1:length(vr_set)
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'CL1',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'Cmx',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'Cmy',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'CL3',1, expName);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'pow',1, expName);
% end
% for i=1
%     plot_para(dataPath, outputPath, 'vr',[i],'CL1',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'CL1',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
%     plot_para(dataPath, outputPath, 'vr',[i],'Cmy',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'Cmy',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
%     plot_para(dataPath, outputPath, 'vr',[i],'Cmx',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'Cmx',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
%     plot_para(dataPath, outputPath, 'vr',[i],'CL3',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'CL3',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
%     plot_para(dataPath, outputPath, 'vr',[i],'pow',2, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, color_set,all_Y);
%     plot_para_2D(dataPath, outputPath, 'vr',[i],'theta',[1, 3, 5, 7, 10, 13, 15, 17],'pow',1, vr_set(vrNumber), theta_set(phaseNumber), X_set, Y_set, all_Y);
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


close all
