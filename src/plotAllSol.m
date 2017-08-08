close all

freqRatioInd = 6;
compIndx = 9:24;

outputDir = '/Users/haining/Dropbox/viv/src/ILCFprediction/output';

figure; 

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp001_run10.mat','Yad_pred_res_min', 'vrn_fr', 'Yad_fr', 'vrn_pred');
plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-')
hold on

plot(vrn_pred(freqRatioInd,compIndx),Yad_pred_res_min(freqRatioInd,compIndx),'r+','MarkerSize',9);

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp01_run10.mat','Yad_pred_res_min', 'vrn_fr', 'Yad_fr', 'vrn_pred');
plot(vrn_pred(freqRatioInd,compIndx),Yad_pred_res_min(freqRatioInd,compIndx),'co','MarkerSize',9);
% 
% load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp1_run10.mat','Yad_pred_res_min', 'vrn_fr', 'Yad_fr', 'vrn_pred');
% plot(vrn_pred(freqRatioInd,compIndx),Yad_pred_res_min(freqRatioInd,compIndx),'g*','MarkerSize',9);

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp2_run10.mat','Yad_pred_res_min', 'vrn_fr', 'Yad_fr', 'vrn_pred');
plot(vrn_pred(freqRatioInd,compIndx),Yad_pred_res_min(freqRatioInd,compIndx),'m<','MarkerSize',9);

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp5_run10.mat','Yad_pred_res_min', 'vrn_fr', 'Yad_fr', 'vrn_pred');
plot(vrn_pred(freqRatioInd,compIndx),Yad_pred_res_min(freqRatioInd,compIndx),'ks','MarkerSize',9);

legend('free vib','Damping = 10^{-4}','Damping = 10^{-3}','Damping = 0.01','Damping = 0.05')
xlabel('V_{rn}','FontSize',16);
ylabel('A_{y}/D','FontSize',16);
set(gca,'FontSize',14);
axis([4 10 0 1.5]);
set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
grid
% suptitle('AmpY Prediction vs Free Virbation Experiment Data')
saveas(gcf,[outputDir filesep 'PredvsFree'  'AmpYAllDamp.png'] );
saveas(gcf,[outputDir filesep 'PredvsFree'  'AmpYAllDamp.fig'] );

    
figure; 

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp001_run10.mat','Xad_pred_res_min', 'vrn_fr', 'Xad_fr', 'vrn_pred');
plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-')
hold on

plot(vrn_pred(freqRatioInd,compIndx),Xad_pred_res_min(freqRatioInd,compIndx),'r+','MarkerSize',9);

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp01_run10.mat','Xad_pred_res_min', 'vrn_fr', 'Xad_fr', 'vrn_pred');
plot(vrn_pred(freqRatioInd,compIndx),Xad_pred_res_min(freqRatioInd,compIndx),'co','MarkerSize',9);

% load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp1_run10.mat','Xad_pred_res_min', 'vrn_fr', 'Xad_fr', 'vrn_pred');
% plot(vrn_pred(freqRatioInd,compIndx),Xad_pred_res_min(freqRatioInd,compIndx),'g*','MarkerSize',9);
% 
load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp2_run10.mat','Xad_pred_res_min', 'vrn_fr', 'Xad_fr', 'vrn_pred');
plot(vrn_pred(freqRatioInd,compIndx),Xad_pred_res_min(freqRatioInd,compIndx),'m<','MarkerSize',9);

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp5_run10.mat','Xad_pred_res_min', 'vrn_fr', 'Xad_fr', 'vrn_pred');
plot(vrn_pred(freqRatioInd,compIndx),Xad_pred_res_min(freqRatioInd,compIndx),'ks','MarkerSize',9);

legend('free vib','Damping = 10^{-4}','Damping = 10^{-3}','Damping = 0.01','Damping = 0.05')
xlabel('V_{rn}','FontSize',16);
ylabel('A_{x}/D','FontSize',16);
set(gca,'FontSize',14);
axis([4 10 0 0.6]);
set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
grid
% suptitle('AmpX Prediction vs Free Virbation Experiment Data')
saveas(gcf,[outputDir filesep 'PredvsFree'  'AmpXAllDamp.png'] );
saveas(gcf,[outputDir filesep 'PredvsFree'  'AmpXAllDamp.fig'] );


figure
load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp001_run10.mat','CL3_pred_res_min', 'vrn_fr', 'CL3_pred_free', 'CL3_fr', 'vrn_pred');
plot(vrn_fr(freqRatioInd,compIndx),CL3_fr(freqRatioInd,compIndx),'bo-')
hold on
plot(vrn_fr(freqRatioInd,compIndx),CL3_pred_free(freqRatioInd,compIndx),'gx-')
plot(vrn_pred(freqRatioInd,compIndx),CL3_pred_res_min(freqRatioInd,compIndx),'r+','MarkerSize',9);

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp01_run10.mat','CL3_pred_res_min');
plot(vrn_pred(freqRatioInd,compIndx),CL3_pred_res_min(freqRatioInd,compIndx),'co','MarkerSize',9);

% load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp1_run10.mat','CL3_pred_res_min');
% plot(vrn_pred(freqRatioInd,compIndx),CL3_pred_res_min(freqRatioInd,compIndx),'g*','MarkerSize',9);

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp2_run10.mat','CL3_pred_res_min');
plot(vrn_pred(freqRatioInd,compIndx),CL3_pred_res_min(freqRatioInd,compIndx),'m<','MarkerSize',9);

load('ICFLMotion_ZDMiT_RobustSmooth_20130906_dampp5_run10.mat','CL3_pred_res_min');
plot(vrn_pred(freqRatioInd,compIndx),CL3_pred_res_min(freqRatioInd,compIndx),'ks','MarkerSize',9);

xlabel('V_{rn}','FontSize',16);
ylabel('CL3','FontSize',16);
set(gca,'FontSize',14);
legend({'Measurement from Free Vib Exp', 'Predicted from Free Vib Motion','Damping = 10^{-4}','Damping = 10^{-3}','Damping = 0.01','Damping = 0.05'})
axis([4 10 0 2]);
set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]);
grid
title(['CL3 Prediction vs Free Virbation Experiment Data '])    
saveas(gcf,[outputDir filesep 'HHPredvsFree'  'CL3AllDamp.png'] );
saveas(gcf,[outputDir filesep 'HHPredvsFree'  'CL3AllDamp.fig'] );
