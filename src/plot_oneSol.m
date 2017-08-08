    function plot_oneSol(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, CL1_fr, CL3_fr, CL1_pred_free, CL3_pred_free, vrn_pred, ...
        vr_pred_oneSol, Yad_pred_oneSol, Xad_pred_oneSol, Theta_pred_oneSol, CL1_pred_oneSol, CL3_pred_oneSol, ...
        freqRatio, freqRatioInd, compIndx, col, shp, outputDir)
        
    
    figure; 
    subplot(2,2,1)
    plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
        vrn_pred(freqRatioInd,compIndx),vr_pred_oneSol(freqRatioInd,compIndx)./vrn_pred(freqRatioInd,compIndx),'m<','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,2)
    plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Theta_pred_oneSol(freqRatioInd,compIndx),'m<','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
%     legend('Observed Free Vib','Predicted from Forced Vib','Location','Southwest')
    axis([4 10 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    subplot(2,2,3)
    plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Yad_pred_oneSol(freqRatioInd,compIndx),'m<','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,4)
    plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Xad_pred_oneSol(freqRatioInd,compIndx),'m<','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('Prediction vs Free Virbation Experiment Data')
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'oneSol.png'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'oneSol.fig'] );
    
    figure; 
    plot(vrn_fr(freqRatioInd,compIndx),CL3_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),CL3_pred_oneSol(freqRatioInd,compIndx),'m<','MarkerSize',9);
    hold on
%     plot(vrn_fr(freqRatioInd,compIndx),CL3_pred_free(freqRatioInd,compIndx),'rx-')
    xlabel('V_{rn}','FontSize',16);
    ylabel('CL3','FontSize',16);
    set(gca,'FontSize',14);
%       legend({'Observed Free Vib', 'Predicted from Forced Vib', 'Predicted from Free Vib'},'Location','Southwest')
    legend({'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 0 2]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]);
    grid
    title(['CL3 Prediction vs Free Virbation Experiment Data ' num2str(freqRatio(freqRatioInd))])    
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'CL3.png'] );
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'CL3.fig'] );
 
    figure; 
    plot(vrn_fr(freqRatioInd,compIndx),CL3_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),CL3_pred_oneSol(freqRatioInd,compIndx),'m<','MarkerSize',9);
    hold on
    plot(vrn_fr(freqRatioInd,compIndx),CL3_pred_free(freqRatioInd,compIndx),'rx-')
    xlabel('V_{rn}','FontSize',16);
    ylabel('CL3','FontSize',16);
    set(gca,'FontSize',14);
    legend({'Measurement from Free Vib Exp', 'Predicted in Forced Vib', 'Predicted from Free Vib Motion'},'Location','Southwest')
    axis([4 10 0 2]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]);
    grid
    title(['CL3 Prediction vs Free Virbation Experiment Data ' num2str(freqRatio(freqRatioInd))])    
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'CL3_all.png'] );
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'CL3_all.fig'] );
 
        
    figure; 
    plot(vrn_fr(freqRatioInd,compIndx),CL1_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),CL1_pred_oneSol(freqRatioInd,compIndx),'m<','MarkerSize',9);
    hold on
%     plot(vrn_fr(freqRatioInd,compIndx),CL1_pred_free(freqRatioInd,compIndx),'rx-')
    xlabel('V_{rn}','FontSize',16);
    ylabel('CL1','FontSize',16);
    set(gca,'FontSize',14);
%         legend({'Observed Free Vib', 'Predicted from Forced Vib', 'Predicted from Free Vib'},'Location','Southwest')
    legend({'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 0 3]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 1 2 3]);
    grid
    title([' CL1 Prediction vs Free Virbation Experiment Data ' num2str(freqRatio(freqRatioInd))])    
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'CL1.png'] );
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'CL1.fig'] );
    
    figure; 
    plot(vrn_fr(freqRatioInd,compIndx),CL3_fr(freqRatioInd,compIndx).^2./(CL1_fr(freqRatioInd,compIndx).^2+CL3_fr(freqRatioInd,compIndx).^2),'bo-',...
        vrn_pred(freqRatioInd,compIndx),CL3_pred_oneSol(freqRatioInd,compIndx).^2./(CL3_pred_oneSol(freqRatioInd,compIndx).^2+CL1_pred_oneSol(freqRatioInd,compIndx).^2),'m<','MarkerSize',9);
    hold on
%     plot(vrn_fr(freqRatioInd,compIndx),CL3_pred_free(freqRatioInd,compIndx).^2./(CL1_pred_free(freqRatioInd,compIndx).^2+CL3_pred_free(freqRatioInd,compIndx).^2),'rx-')
    xlabel('V_{rn}','FontSize',16);
    ylabel('A3/A','FontSize',16);
    set(gca,'FontSize',14);
%         legend({'Observed Free Vib', 'Predicted from Forced Vib', 'Predicted from Free Vib'},'Location','Southwest')
    legend({'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 0 1]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1]);
    grid
    title([' A3/A Prediction vs Free Virbation Experiment Data ' num2str(freqRatio(freqRatioInd))])    
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'A3.png'] );
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'A3.fig'] );
 
 
    

  