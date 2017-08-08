function plot_addedMass(CLv_pred, CDv_pred, Cmy_pred, Cmx_pred, vr2D_pred_x, Vrn_y, Vrn_x, ...
    massRatio_y, massRatio_x, vrn_pred,vr2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)


% Added mass verification
    
    figure; 
    subplot(2,3,1)
    for i=compIndx
        for j = 1:length(CLv_pred{i})
%             plot(vrn_pred(freqRatioInd,i),CLv_pred{i},'Color', col(j,:),'Marker',shp(j),'MarkerSize',9);
            plot(vrn_pred(freqRatioInd,i),CLv_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('CLv Y','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
%     axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    
    subplot(2,3,4)
    for i=compIndx
        for j = 1:length(CDv_pred{i})
            plot(vrn_pred(freqRatioInd,i),CDv_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('CLv X','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
%     axis([4 10 -180 180]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    
    subplot(2,3,3)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j)./Vrn_y(freqRatioInd,i)./sqrt((massRatio_y(freqRatioInd)+Cmy_pred{i}(j))...
                ./(massRatio_y(freqRatioInd)+1)),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{r}(addedmass) Y','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 0.90 1.1]);
    set(gca,'Xtick',[4 5 6 7 8 9 10],'Ytick',[0.95 1 1.05]);
    grid
    
    subplot(2,3,6)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred_x{i}(j)./Vrn_x(freqRatioInd,i)./sqrt((massRatio_x(freqRatioInd)+Cmx_pred{i}(j))...
                ./(massRatio_x(freqRatioInd)+1)),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr2D_pred_x{i}./vrn_pred(freqRatioInd,i).*sqrt((massRatio_x(freqRatioInd)+Cmy_pred{i})./(massRatio_x(freqRatioInd)+1)),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{r}(addedmass) X','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 0.90 1.1]);
    set(gca,'Xtick',[4 5 6 7 8 9 10],'Ytick',[0.95 1 1.05]);
    grid
    
     subplot(2,3,2)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('Cmy','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
%     axis([4 10 0.90 1.1]);
    xlim([4 10])
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    
    subplot(2,3,5)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr2D_pred_x{i}./vrn_pred(freqRatioInd,i).*sqrt((massRatio_x(freqRatioInd)+Cmy_pred{i})./(massRatio_x(freqRatioInd)+1)),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('Cmx','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
%     axis([4 10 0.90 1.1]);
    xlim([4 10])
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid   
    
    suptitle('Prediction Verification')    
    saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'all.jpg'] );
    saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'all.fig'] );
        
       %% Added mass verification
    figure
    subplot(2,3,1)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlim([4 10])
    xlabel('V_{rn}','FontSize',16);
    ylabel('Cmy','FontSize',16);
    set(gca,'FontSize',14);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid    
    subplot(2,3,4)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlim([4 10])
    xlabel('V_{rn}','FontSize',16);
    ylabel('Cmx','FontSize',16);
    set(gca,'FontSize',14);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid   
    subplot(2,3,2)
    for i=compIndx
        for j = 1:length(vr2D_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end  
    xlim([4 10])
    xlabel('V_{rn}','FontSize',16);
    ylabel('Vry','FontSize',16);
    set(gca,'FontSize',14);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid    
    subplot(2,3,5)
    for i=compIndx
        for j = 1:length(vr2D_pred_x{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred_x{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('Vrx','FontSize',16);
    xlim([4 10])
    set(gca,'FontSize',14);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    subplot(2,3,3)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j)./Vrn_y(freqRatioInd,i)./sqrt((massRatio_y(freqRatioInd)+Cmy_pred{i}(j))...
                ./(massRatio_y(freqRatioInd)+1)),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{r}(addedmass) Y','FontSize',16);
    set(gca,'FontSize',14);
    axis([4 10 0.90 1.1]);
    set(gca,'Xtick',[4 5 6 7 8 9 10],'Ytick',[0.95 1 1.05]);
    grid
    
    subplot(2,3,6)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred_x{i}(j)./Vrn_x(freqRatioInd,i)./sqrt((massRatio_x(freqRatioInd)+Cmx_pred{i}(j))...
                ./(massRatio_x(freqRatioInd)+1)),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{r}(addedmass) X','FontSize',16);
    set(gca,'FontSize',14);
    axis([4 10 0.90 1.1]);
    set(gca,'Xtick',[4 5 6 7 8 9 10],'Ytick',[0.95 1 1.05]);
    grid
    suptitle('Prediction Verification Cm Vr')    
    saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'CmVr.jpg'] );
    saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'CmVr.fig'] );
    