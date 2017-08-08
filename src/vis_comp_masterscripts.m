% function vis_comp_masterscripts()

utility = Utility.getInstance();
filesep = utility.getFilesep();

col = distinguishable_colors(100);
shp = {'+','o','*','s','d','^','v','>','<','p','h','x','.','+','o','*','s','d','^','v','>','<','p','h','x','.'};

[vivaRoot, optRoot, experimentRoot, fatigueRoot, utilityRoot, InlineVIVRoot, CrossflowVIVRoot, mooringRoot, ILCFVIVRoot] = getRoots();
outputDir =[ILCFVIVRoot filesep 'outputComp'];
utility.conditionalMkdir(outputDir);


bb_viva = load('basic_bare');
basic_bare_inline_MT = load('basic_bare_inline_MT');
basic_bare_inline_MT2 = load('basic_bare_inline_MT2');

color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
targetpar_set = {'CL1','CL3','CD2','Cmy', 'Cmx','pow','CDv'}; 


vr_set = {'2' '2.5' '3' '3.5' '4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8' '8.5' '9' '9.5' '10'}; 
X_set = {'0','0.05','0.1','0.15','0.2','0.25','0.3','0.45','0.6','0.75'};
Y_set = {'0.15','0.25','0.5','0.75', '1', '1.25', '1.5'};
% theta_set = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
%      '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
theta_set = {'180n' '157.5n' '135n' '120n' '105n' '90n' '75n' '60n'  '45n' '22.5n' '0'...
                '22.5' '45' '67.5' '90' '112.5' '135' '157.5' '180'};
            
Y_vec = [0.15, 0.25, 0.5, 0.75, 1, 1.25, 1.5 ];
X_vec = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.45,0.6,0.75 ];
% theta_vec = [-180:22.5:-45, -30:15:30, 45:22.5:180];
theta_vec =  [-180:22.5:-135, -120:15:-45, -22.5:22.5:180];
vr_vec_ori = [4:0.5:8];


% JMD Exp

% 6*6*8*8 = 2304  6*6*9*8 = 2592
Y_set_JMD = {'0.25','0.5','0.75','1', '1.25', '1.5'};
X_set_JMD = {'0','0.15','0.3','0.45','0.6','0.75'};
theta_set_JMD = {'180n' '135n' '90n' '45n' '0' '45' '90' '135' '180'};
vr_set_JMD = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'}; 


%  HZ Exp
    
vr_set_HZ = {'4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8'}; 
X_set_HZ = {'0.05','0.1','0.15','0.2','0.25','0.3'};
Y_set_HZ = {'0.15','0.25','0.5','0.75', '1'};
% theta_set_HZ = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
%      '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
theta_set_HZ = {'180n' '157.5n' '135n' '120n' '105n' '90n' '75n' '60n'  '45n' '22.5n' '0'...
                '22.5' '45' '67.5' '90' '112.5' '135' '157.5' '180'};

% combine 2 database

clear ILCFHydroModelHZobj
ILCFHydroModelHZobj = ILCFHydroModelHZ();

clear ILCFHydroModelJMDobj
ILCFHydroModelJMDobj = ILCFHydroModelJMD();

clear ILCFHydroModelZDMiTobj
ILCFHydroModelZDMiTobj = ILCFHydroModelZDMiT();

%% plot all database
vr_vec_disp = 3:0.5:10;
theta_vec_disp = [-180 -135 -90 -45 0 45 90 135];
theta_set_disp = {'180n' '135n' '90n' '45n' '0' '45' '90' '135'};

y_lim_CLv = [-1 2.5];
y_lim_CDv = [-1 1];
y_lim_Cmy = [-2 4];
y_lim_Cmx = [-2 4];
x_lim = [0.1 0.3];

marker = 'Combined';
plotDatabase(ILCFHydroModelZDMiTobj,Y_vec, X_vec, theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,...
    y_lim_CDv,x_lim, bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)

marker = 'JMD';
plotDatabase(ILCFHydroModelJMDobj,Y_vec, X_vec, theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,...
    y_lim_CDv,x_lim, bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)

marker = 'HZ';
plotDatabase(ILCFHydroModelHZobj,Y_vec, X_vec, theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,...
    y_lim_CDv,x_lim, bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)



Y_vec_comp = [0.25,0.5,0.75,1,0.75,1];
X_vec_comp = [0.15,0.15,0.15,0.15,0.3,0.3];

marker = 'HZvsJMD';
plotDatabaseComp(ILCFHydroModelHZobj,ILCFHydroModelJMDobj,Y_vec_comp, X_vec_comp, theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,...
    y_lim_CDv,x_lim, bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)

marker = 'OriHZvsJMD';
plotDatabaseCompOri(ILCFHydroModelHZobj,ILCFHydroModelJMDobj,Y_vec_comp, X_vec_comp, theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,...
    y_lim_CDv,x_lim, bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)


marker = 'Remi';
Plot_Remi(0.5, 0.05,  theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)


%% compare

Y_vec_comp = [0.25,0.75];
X_vec_comp = [0.15,0.3];

marker = 'Cmy';
plotDatabaseCompAll(ILCFHydroModelZDMiTobj,ILCFHydroModelHZobj,ILCFHydroModelJMDobj,marker,Y_vec_comp, X_vec_comp,theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)
marker = 'Cmx';
plotDatabaseCompAll(ILCFHydroModelZDMiTobj,ILCFHydroModelHZobj,ILCFHydroModelJMDobj,marker,Y_vec_comp, X_vec_comp,theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)
marker = 'CLv';
plotDatabaseCompAll(ILCFHydroModelZDMiTobj,ILCFHydroModelHZobj,ILCFHydroModelJMDobj,marker,Y_vec_comp, X_vec_comp,theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)
marker = 'CDv';
plotDatabaseCompAll(ILCFHydroModelZDMiTobj,ILCFHydroModelHZobj,ILCFHydroModelJMDobj,marker,Y_vec_comp, X_vec_comp,theta_vec_disp, vr_vec_disp, theta_set_disp, y_lim_Cmy,y_lim_Cmx, y_lim_CLv,y_lim_CDv,x_lim, ...
    bb_viva, basic_bare_inline_MT, col, shp, outputDir, marker)


%% 2D plots
y_val =0.5;
x_val = 0.05;

load('proc_Haining_052213.mat')


figure
subplot(221)
plot2D_theta_vr(ILCFHydroModelZDMiTobj , 'Cmy', y_val, x_val, theta_vec_disp, vr_vec_disp);
title('Combined')
subplot(222)
plot2D_theta_vr(ILCFHydroModelHZobj , 'Cmy', y_val, x_val, theta_vec_disp, vr_vec_disp);
title('Combined')
subplot(223)
plot2D_theta_vr(ILCFHydroModelJMDobj , 'Cmy', y_val, x_val, theta_vec_disp, vr_vec_disp);
title('Combined')
subplot(224)
surface(Ur,a,cmx_tot);
title('Remi')
                                

% %% Dataset 1
% vr_set_1 = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'}; 
% X_set_1 = {'0','0.15','0.3','0.45','0.6','0.75'};
% Y_set_1 = {'0.25','0.5','0.75','1', '1.25', '1.5'};
% theta_set_1 = {'180n' '135n' '90n' '45n' '0' '45' '90' '135'};
% 
% all_Y_1 = [0.25:0.25:1.5];
% 
% X_vector_1 = [0:0.15:0.75 [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75] [0:0.15:0.75]];
% Y_vector_1 = [0.25*ones(1,6) 0.5*ones(1,6) 0.75*ones(1,6) 1*ones(1,6) 1.25*ones(1,6) 1.5*ones(1,6)];
% 
% dataPath_1 = '/Users/haining/VIV/src/ILCFprediction/databaseJMD/';
% 
% vrIndex_1 = [4];
% AmpXIndex_1 = [1:6];
% AmpYIndex_1 = [1:6];
% targetpar_1 = 'CL1';
% 
% %% Dataset 2
% 
% vr_set_2 = {'4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8'}; 
% X_set_2 = {'0.05','0.1','0.15','0.2','0.25','0.3'};
% Y_set_2 = {'0.15','0.25','0.5','0.75', '1'};
% theta_set_2 = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
%      '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
%  
% all_Y_2 = [0.15,0.25,0.5,0.75, 1 ];
% 
% Y_vector_2 = [0.15 0.15 0.25 0.25 0.25 0.5 0.5 0.5 0.5 0.75 0.75 0.75 0.75 0.75 0.75 1 1 1 1 1 1];
% X_vector_2 = [0.05 0.1 0.05 0.1 0.15 0.05 0.1 0.15 0.2 0.05 0.1 0.15 0.2 0.25 0.3 0.05 0.1 0.15 0.2 0.25 0.3 ];
% 
% dataPath_2 =  '/Users/haining/VIV/src/ILCFprediction/databaseHZ/';
% 
% vrIndex_2 = [5];
% AmpXIndex_2 = [1:6];
% AmpYIndex_2 = [1:5];
% targetpar_2 = 'CL1';
% 
% %% combined plot
% 
% JasonThetaIndex = [1, 3, 5, 7, 10, 13, 15, 17];
% repeatExpAmpXIndex_1 = [2 2 2 2 3 3];
% repeatExpAmpYIndex_1 = [1 2 3 4 3 4];
% repeatExpAmpXIndex_2 = [3 3 3 3 6 6];
% repeatExpAmpYIndex_2 = [2 3 4 5 4 5];
% 
% X_set_3 = {'0','0.05','0.1','0.15','0.2','0.25','0.3','0.45','0.6','0.75'};
% Y_set_3 = {'0.15','0.25','0.5','0.75', '1', '1.25', '1.5'};
% 
% all_Y_3 = [0.15,0.25,0.5,0.75, 1, 1.25, 1.5 ];
% 
% 
% 
% for i = 1:length(JasonThetaIndex)
%     
%     thetaIndex_1 = i;
%     thetaIndex_2 = JasonThetaIndex(i);
%     
%     [value_1, vr_1, theta_1, X_1, Y_1] = extract_para_2D(dataPath_1, vrIndex_1, thetaIndex_1, ...
%     AmpXIndex_1, AmpYIndex_1, targetpar_1, vr_set_1, theta_set_1, X_set_1, Y_set_1, all_Y_1);
% 
% 
%     [value_2, vr_2, theta_2, X_2, Y_2] = extract_para_2D(dataPath_2, vrIndex_2, thetaIndex_2, ...
%         AmpXIndex_2, AmpYIndex_2, targetpar_2, vr_set_2, theta_set_2, X_set_2, Y_set_2, all_Y_2);
% 
% 
%     newValueMatrix = NAN(length(Y_set_3), length(X_set_3));
%     
%     % average value at same position
%     
%     % plot big matrix
%     
%     figure
%     [C h] = contour( X_1, Y_1, value_1);
%     set(h,'ShowText','on')
%     hold on
%     [C h] = contour( X_2, Y_2, value_2,'--');
%     set(h,'ShowText','on')
%     xlabel('Inline Motion')
%     ylabel('Crossflow Motion')
%     colorbar
%     
%     plot(X_vector_1, Y_vector_1, 'ob')
%     plot(X_vector_2, Y_vector_2, '*r')
% %     xlim([min(0,0.9*min(X_vector_1)-0.1) 1.1*max(X_vector_1)])
% %     ylim([min(0,0.9*min(Y_vector_1)-0.1) 1.1*max(Y_vector_1)])  
%     title([targetpar_2 ' \theta: ', theta_set_2{thetaIndex_2}, ' Vr: ', vr_set_2{vrIndex_2}])
%     saveas(gcf, [outputDir 'contourLineWMatrixComparison' 'theta' theta_set_2{thetaIndex_2} 'vr' vr_set_2{vrIndex_2} '.jpg'])
%     close
% end
% 
% 
% %% contour comparison
% JasonThetaIndex = [1, 3, 5, 7, 10, 13, 15, 17];
% 
% for i = 1:length(JasonThetaIndex)
%     
%     thetaIndex_1 = i;
%     thetaIndex_2 = JasonThetaIndex(i);
%     
%     [value_1, vr_1, theta_1, X_1, Y_1] = extract_para_2D(dataPath_1, vrIndex_1, thetaIndex_1, ...
%     AmpXIndex_1, AmpYIndex_1, targetpar_1, vr_set_1, theta_set_1, X_set_1, Y_set_1, all_Y_1);
% 
%     [value_2, vr_2, theta_2, X_2, Y_2] = extract_para_2D(dataPath_2, vrIndex_2, thetaIndex_2, ...
%         AmpXIndex_2, AmpYIndex_2, targetpar_2, vr_set_2, theta_set_2, X_set_2, Y_set_2, all_Y_2);
% 
%     figure
%     [C h] = contour( X_1, Y_1, value_1);
%     set(h,'ShowText','on')
%     hold on
%     [C h] = contour( X_2, Y_2, value_2,'--');
%     set(h,'ShowText','on')
%     xlabel('Inline Motion')
%     ylabel('Crossflow Motion')
%     colorbar
%     
%     plot(X_vector_1, Y_vector_1, 'ob')
%     plot(X_vector_2, Y_vector_2, '*r')
% %     xlim([min(0,0.9*min(X_vector_1)-0.1) 1.1*max(X_vector_1)])
% %     ylim([min(0,0.9*min(Y_vector_1)-0.1) 1.1*max(Y_vector_1)])  
%     title([targetpar_2 ' \theta: ', theta_set_2{thetaIndex_2}, ' Vr: ', vr_set_2{vrIndex_2}])
%     saveas(gcf, [outputDir 'contourLineWMatrixComparison' 'theta' theta_set_2{thetaIndex_2} 'vr' vr_set_2{vrIndex_2} '.jpg'])
%     close
% end
% 
% %% comparison of two curves
% 
% repeatExpAmpXIndex_1 = [2 2 2 2 3 3];
% repeatExpAmpYIndex_1 = [1 2 3 4 3 4];
% 
% repeatExpAmpXIndex_2 = [3 3 3 3 6 6];
% repeatExpAmpYIndex_2 = [2 3 4 5 4 5];
% 
% for i = 1:length(JasonThetaIndex)
%     
%     thetaIndex_1 = i;
%     thetaIndex_2 = JasonThetaIndex(i);
%     
%     for j = 1:length(repeatExpAmpXIndex_1)
% 
%         [value_1(i,j)] = extract_para_2D(dataPath_1, vrIndex_1, thetaIndex_1, ...
%         AmpXIndex_1(repeatExpAmpXIndex_1(j)), AmpYIndex_1(repeatExpAmpYIndex_1(j)),...
%         targetpar_1, vr_set_1, theta_set_1, X_set_1, Y_set_1, all_Y_1,1);
% 
% 
%         [value_2(i,j)] = extract_para_2D(dataPath_2, vrIndex_2, thetaIndex_2, ...
%             AmpXIndex_2(repeatExpAmpXIndex_2(j)), AmpYIndex_2(repeatExpAmpYIndex_2(j)),...
%             targetpar_2, vr_set_2, theta_set_2, X_set_2, Y_set_2, all_Y_2,1);
%     end
% end
% 
% figure
% for i=1:length(JasonThetaIndex)
%     subplot(2,4,i)
%     plot(value_1(i,:), 'o-')
%     hold on
%     plot(value_2(i,:), '*-r')
%     title(['theta' theta_set_2{JasonThetaIndex(i)}])
% end
% suptitle([targetpar_2 ' Comparison of Value at vr ' vr_set_2{vrIndex_2} ' Blue: Dahl Red: New']);
% saveas(gcf, [outputDir targetpar_2 'contourLineComparison' 'theta' theta_set_2{thetaIndex_2} 'vr' vr_set_2{vrIndex_2} '.jpg'])  
% 
% 
% close all
freqRatio =[1	1.22	1.37	1.52	1.67	1.9];
freqRatioNames ={ '1p0', '1p22', '1p37', '1p52', '1p67', '1p9'};
fn_y = [0.715	0.799	0.894	0.977	0.698	0.704];
fn_x = [0.715	0.974   1.226   1.488   1.167   1.336];

springCons_y = [697  902 1118 1386 967 1013];
springCons_x = [609  1301 2130 2850 2580 3235];

displacedMass = fluidDensity * pi*(Diameter/2)^2*Riserlength;

massRatio_y = (springCons_y./(fn_y*2*pi).^2)/displacedMass-1;
massRatio_x = (springCons_x./(fn_x*2*pi).^2)/displacedMass-1;
% 
% % begin added for mass ratio 1
% mRatio = massRatio_y(6);
% massRatio_y = 1 *ones(size(massRatio_y));
% massRatio_x = 1 *ones(size(massRatio_x));
% % end added for mass ratio 1
    
    
% massRatio_y = (springCons_y./(fn_y*2*pi).^2)/displacedMass;
% massRatio_x = (springCons_x./(fn_x*2*pi).^2)/displacedMass;

strucDampRatio_y =0.00001*[2.2	1.3	1.1	1.6	2.6	6.2];
strucDampRatio_x =0.00001*[2.2	1.7	2.5	3.2	2.9	2.5];

% massRatio_x =[3.3 3.8 3.7 3.6 5.3 5]; % this massRatio in paper contains
% normal added mass 1
% massRatio_y =[3.8 3.9 3.9 4   5.5	5.7];

massDensity_x = fluidDensity.*massRatio_x;
massDensity_y = fluidDensity.*massRatio_y;


damperCons_x = strucDampRatio_x*2.*(1+massRatio_x)* displacedMass*2*pi.*fn_x;
damperCons_y = strucDampRatio_y*2.*(1+massRatio_y)* displacedMass*2*pi.*fn_y;

freeVibDataDir = ['..' filesep 'freeVibJMD' filesep];
load([freeVibDataDir 'vel.mat']);

newmatfiles = {'new1p0' 'new1p22' 'new1p37' 'new1p52' 'new1p67' 'new1p9'};
phasematfiles = {'phase1p0' 'phase1p22' 'phase1p37' 'phase1p52' 'phase1p67' 'phase1p9'};
corrmatfiles = {'corr1p0' 'corr1p22' 'corr1p37' 'corr1p52' 'corr1p67' 'corr1p9'};
ampfixmatfiles = {'ampfix1p0' 'ampfix1p22' 'ampfix1p37' 'ampfix1p52' 'ampfix1p67' 'ampfix1p9'};

compIndx = [9:25];


% for freqRatioInd= 4:6 
for freqRatioInd= 6
    
    eval(['FluidVel = vels' freqRatioNames{freqRatioInd} ';'])
    eval(['load ' freeVibDataDir char(newmatfiles(freqRatioInd)) '.mat;']);
    eval(['load ' freeVibDataDir char(phasematfiles(freqRatioInd)) '.mat;']);
    phase_ind = find(phasextoy > 180);
    phasextoy(phase_ind) = phasextoy(phase_ind)-360;
    eval(['load ' freeVibDataDir char(corrmatfiles(freqRatioInd)) '.mat;']);
    eval(['load ' freeVibDataDir char(ampfixmatfiles(freqRatioInd)) '.mat;']); 
    eval(['load ' freeVibDataDir char(corrmatfiles(freqRatioInd)) '.mat;']);
    eval(['load ' freeVibDataDir char(ampfixmatfiles(freqRatioInd)) '.mat;']);
% 
    vr_free = Vrn.*fy./(ypeakfreq./2./pi);

    vr_fr(freqRatioInd,:) = vr_free;
    Yad_fr(freqRatioInd,:) = Yad_fix;
    Xad_fr(freqRatioInd,:) = Xad_fix;
    Theta_fr(freqRatioInd,:) = phasextoy;
    vrn_fr(freqRatioInd,:) = Vrn;
end

% plot the added mass of free vibrations in the forced vibration database

load([freeVibDataDir 'freeVibForceCoef.mat'])


for i = compIndx
    Cmy_fr_JMD(i) = getDataPoint(ILCFHydroModelJMDobj , 'Cmy', Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i), Theta_fr(freqRatioInd,i)/180, vr_fr(freqRatioInd,i));
    Cmy_fr_HZ(i) = getDataPoint(ILCFHydroModelHZobj , 'Cmy', Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i), Theta_fr(freqRatioInd,i)/180, vr_fr(freqRatioInd,i));
    Cmy_fr_ZDMiT(i) = getDataPoint(ILCFHydroModelZDMiTobj , 'Cmy', Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i), Theta_fr(freqRatioInd,i)/180, vr_fr(freqRatioInd,i));
end
figure
plot(Vrn_all(freqRatioInd,compIndx), Cmy_all(freqRatioInd,compIndx),'ko-')
hold on
plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_ZDMiT(compIndx),'g*-')
% plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_JMD(compIndx),'rs-')
% plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_HZ(compIndx),'b+-')
legend('Free Vibration Experiment', 'Forced Database','JMD Database', 'HZ Database')
title('Cmy')


figure
plot(Vrn_all(freqRatioInd,compIndx), Cmy_all(freqRatioInd,compIndx),'ko-')
hold on
for i =compIndx
    plot(Vrn_all(freqRatioInd,i), Cmy_pred{i}(1),'g*-')
    hold on
end
% plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_JMD(compIndx),'rs-')
% plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_HZ(compIndx),'b+-')
legend('Free Vibration Experiment', 'Forced Database')
title('Cmy')

figure
plot(Vrn_all(freqRatioInd,compIndx), Cmy_all(freqRatioInd,compIndx)./Cmx_all(freqRatioInd,compIndx),'ko-')
hold on
for i =compIndx
    plot(Vrn_all(freqRatioInd,i), Cmy_pred{i}(1)/Cmx_pred{i}(1),'g*-')
    hold on
end
% plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_JMD(compIndx),'rs-')
% plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_HZ(compIndx),'b+-')
legend('Free Vibration Experiment', 'Forced Database')
title('Cmy/Cmx')

figure
plot(Vrn_all(freqRatioInd,compIndx), sqrt((massDensity_y(freqRatioInd)+Cmy_all(freqRatioInd,compIndx))./(massDensity_x(freqRatioInd)+Cmx_all(freqRatioInd,compIndx))),'ko-')
hold on
for i =compIndx
    plot(Vrn_all(freqRatioInd,i), sqrt((massDensity_y(freqRatioInd)+Cmy_pred{i}(1))/(massDensity_x(freqRatioInd)+Cmx_pred{i}(1))),'g*-')
    hold on
end
% plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_JMD(compIndx),'rs-')
% plot(vrn_fr(freqRatioInd, compIndx),Cmy_fr_HZ(compIndx),'b+-')
legend('Free Vibration Experiment', 'Forced Database')
title('fx/fy')
