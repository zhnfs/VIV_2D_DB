function plot_vr


figure; 
    for i=compIndx
        for j = 1:length(vr2D_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 4 10]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10]);
    grid
    suptitle('Vr Prediction vs Free Virbation Experiment Data')
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'Vr.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'Vr.fig'] );
