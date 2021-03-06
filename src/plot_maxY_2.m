    function plot_maxY_2(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred, vr_pred_MaxAmpy, Yad_pred_MaxAmpy, Xad_pred_MaxAmpy, Theta_pred_MaxAmpy, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

    %% plot max amp y solution
    
    figure;  
    
    hold on   
    plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
        vrn_pred(freqRatioInd,compIndx),vr_pred_MaxAmpy(freqRatioInd,compIndx)./vrn_pred(freqRatioInd,compIndx),'m<','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid

    
    figure
    hold on
    plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Theta_pred_MaxAmpy(freqRatioInd,compIndx),'m<','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib','Predicted from Forced Vib','Location','Southwest')
    axis([4 10 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    
    figure
    hold on
    plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Yad_pred_MaxAmpy(freqRatioInd,compIndx),'m<','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid

    
    figure
    
    hold on
    plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Xad_pred_MaxAmpy(freqRatioInd,compIndx),'m<','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    
    
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'MaxAmpy.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'MaxAmpy.fig'] );
    
  