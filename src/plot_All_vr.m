function plot_All_vr(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred,vr2D_pred, Ampy_pred, Ampx_pred, theta2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

%% All solutions No color coding
    figure; 
    subplot(2,2,1)
    for i=compIndx
        for j = 1:length(vr2D_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j)./vrn_pred(freqRatioInd,i),'rs','MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,2)
    for i=compIndx
        for j = 1:length(theta2D_pred{i})
            plot(vr2D_pred{i}(j) ,theta2D_pred{i}(j),'rs','MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vr_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{r}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 8.5 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    
    subplot(2,2,3)
    for i=compIndx
        for j = 1:length(Ampy_pred{i})
            plot(vr2D_pred{i}(j),Ampy_pred{i}(j),'rs','MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vr_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{r}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 8.5 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,4)
    for i=compIndx
        for j = 1:length(Ampx_pred{i})
            plot(vr2D_pred{i}(j),Ampx_pred{i}(j),'rs','MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vr_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{r}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 8.5 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('All Prediction vs Free Virbation Experiment Data')    
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'allVr.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'allVr.fig'] );
    
    
%% All solutions color coding
    figure; 
    subplot(2,2,1)
    for i=compIndx
        for j = 1:length(vr2D_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j)./vrn_pred(freqRatioInd,i),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,2)
    for i=compIndx
        for j = 1:length(theta2D_pred{i})
            plot(vr2D_pred{i}(j) ,theta2D_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vr_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{r}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 8.5 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    
    subplot(2,2,3)
    for i=compIndx
        for j = 1:length(Ampy_pred{i})
            plot(vr2D_pred{i}(j),Ampy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vr_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{r}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 8.5 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,4)
    for i=compIndx
        for j = 1:length(Ampx_pred{i})
            plot(vr2D_pred{i}(j),Ampx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vr_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{r}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 8.5 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('All Prediction vs Free Virbation Experiment Data')    
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'allColorVr.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'allColorVr.fig'] );
   