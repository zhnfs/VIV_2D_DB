function plot_min_res(vrn_fr, vr_fr, Yad_fr, Xad_fr, Theta_fr, vrn_pred, vr_pred_res_min, Yad_pred_res_min, Xad_pred_res_min, Theta_pred_res_min, x, res, numRun, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   
   % minimum residuals
    for i=compIndx  
        casenums = [1:numRun];
%         casenums = find(res(freqRatioInd,i,casenums)<1e-8);
        figure
        plot(reshape(x(freqRatioInd,i,casenums,1),1,length(casenums)),'-o')
        hold on
        plot(reshape(x(freqRatioInd,i,casenums,2),1,length(casenums)),'-*r')
        plot(reshape(x(freqRatioInd,i,casenums,3),1,length(casenums))./180,'-*g')
        plot(reshape(x(freqRatioInd,i,casenums,4),1,length(casenums))./vrn_pred(i),'-*m')
        legend('Ay','Ax','theta/180','Vr/Vrn')
        saveas(gcf, [outputDir filesep 'ResultRandomICcase' num2str(i) '.jpg'] )
        saveas(gcf, [outputDir filesep 'ResultRandomICcase' num2str(i) '.fig'] )
        close
        
        figure
        plot(reshape(res(freqRatioInd, i,casenums),1, length(casenums)),'-o')
        title(['Residual Case' num2str(i)])
        saveas(gcf, [outputDir filesep 'ResidualRandomIC' num2str(freqRatio(freqRatioInd)) 'case' num2str(i) '.jpg'] )
        saveas(gcf, [outputDir filesep 'ResidualRandomIC' num2str(freqRatio(freqRatioInd)) 'case' num2str(i) '.fig'] )  
        close
    end

  
        
    figure; 
    subplot(2,2,1)
    plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
        vrn_pred(freqRatioInd,compIndx),vr_pred_res_min(freqRatioInd,compIndx)./vrn_pred(freqRatioInd,compIndx),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,2)
    plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Theta_pred_res_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib','Predicted from Forced Vib','Location','Southwest')
    axis([4 10 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    subplot(2,2,3)
    plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Yad_pred_res_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,4)
    plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Xad_pred_res_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('Min Residual Prediction vs Free Virbation Experiment Data')
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'minRes.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'minRes.fig'] );
