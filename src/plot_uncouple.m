function plot_uncouple(vrn_fr, vr_fr, Yad_fr, Xad_fr, Vrn_y, Vr0_y,Vr0_x, NonDimAmp_linear0_y,NonDimAmp_linear0_x, freqRatio, freqRatioInd, compIndx, outputDir)

%     compIndx = 3:37;
    figure; 
    subplot(2,2,1)
    plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
        Vrn_y(freqRatioInd,compIndx),Vr0_y(freqRatioInd,compIndx)./Vrn_y(freqRatioInd,compIndx),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Crossflow Predicted from Forced Vib')
    axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,2)
    plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
        Vrn_y(freqRatioInd,compIndx),2.*Vr0_x(freqRatioInd,compIndx)./Vrn_y(freqRatioInd,compIndx),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', '2*Inline Predicted from Forced Vib')
    axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,3)
    plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',Vrn_y(freqRatioInd,compIndx),NonDimAmp_linear0_y(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,4)
    plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',Vrn_y(freqRatioInd,compIndx),NonDimAmp_linear0_x(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 0.6]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('Uncoupled Prediction vs Free Virbation Experiment Data')
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'unCoupled.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'unCoupled.fig'] );
   