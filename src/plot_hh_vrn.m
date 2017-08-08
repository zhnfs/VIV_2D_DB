function plot_hh_vrn(vrn_fr, CL3_fr, CL5_fr, vrn_pred, CL3_pred, CL5_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)


              
    figure; 
    subplot(1,2,1)
    for i=compIndx
        for j = 1:length(CL3_pred{i})
            plot(vrn_pred(freqRatioInd,i),CL3_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),CL3_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('CL3','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,2)
    for i=compIndx
        for j = 1:length(CL5_pred{i})
            plot(vrn_pred(freqRatioInd,i),CL5_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),CL5_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('CL5','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    
    suptitle(['All Prediction vs Free Virbation Experiment Data ' num2str(freqRatio(freqRatioInd))])    
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'allColor.jpg'] );
    saveas(gcf,[outputDir filesep 'HHPredvsFree' num2str(freqRatio(freqRatioInd)) 'allColor.fig'] );
    