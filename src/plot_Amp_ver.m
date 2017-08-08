function plot_Amp_ver(CLv_pred, CDv_pred, Cmy_pred, Cmx_pred, vr2D_pred_x, Vrn_y, Vrn_x, massRatio_y, massRatio_x,...
    Damp_pred_y, Ampy_pred, Ampx_pred, Damp_pred_x, vr2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)


%% Individual Amplitude veification    
   for i= compIndx
        figure; 
        subplot(4,3,1) 
        for j = 1:length(CLv_pred{i})
            plot(j, CLv_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
        ylabel('CLv Y','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid
        subplot(4,3,7)
        for j = 1:length(CDv_pred{i})
            plot(j, CDv_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('CLv X','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid

        subplot(4,3,2)
        for j = 1:length(Cmy_pred{i})
            plot(j,Damp_pred_y{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end  
        ylabel('Damp y','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid

        subplot(4,3,8)
        for j = 1:length(Cmx_pred{i})
            plot(j,Damp_pred_x{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Damp x','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid 

        subplot(4,3,3)
        for j = 1:length(Cmy_pred{i})
            plot(j,Cmy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Cmy','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid

        subplot(4,3,9)
        for j = 1:length(Cmx_pred{i})
            plot(j,Cmx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Cmx','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid   
        
        subplot(4,3,4) 
        for j = 1:length(vr2D_pred{i})
            plot(j, vr2D_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
        ylabel('Vr Y','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(vr2D_pred{i})+1])
        grid
        subplot(4,3,10)
        for j = 1:length(vr2D_pred_x{i})
            plot(j, vr2D_pred_x{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Vr X','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(vr2D_pred_x{i})+1])
        grid
        subplot(4,3,5) 
        for j = 1:length(Ampy_pred{i})
            plot(j, Ampy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
        ylabel('Amp Y','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(Ampy_pred{i})+1])
        grid
        subplot(4,3,11)
        for j = 1:length(Ampx_pred{i})
            plot(j, Ampx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Amp X','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(Ampx_pred{i})+1])
        grid        
        
        subplot(4,3,6)
        for j = 1:length(Cmy_pred{i})
            plot(j, Ampy_pred{i}(j)./(CLv_pred{i}(j).*vr2D_pred{i}(j).*Vrn_y(freqRatioInd,i)/...
                (4*pi^3*(massRatio_y(freqRatioInd)+Cmy_pred{i}(j)).*Damp_pred_y{i}(j))),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
        ylabel('Amp Y / Amp Y(lift, added mess)','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(Cmy_pred{i})+1])
        ylim([0.95 1.05])
        grid

        subplot(4,3,12)
        for j = 1:length(Cmx_pred{i})
            plot(j, Ampx_pred{i}(j)./(CDv_pred{i}(j).*vr2D_pred_x{i}(j).*Vrn_x(freqRatioInd,i)/...
                (4*pi^3*(massRatio_x(freqRatioInd)+Cmx_pred{i}(j)).*Damp_pred_x{i}(j))),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end   
        ylabel('Amp X / Amp X(lift, added mess)','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(Cmx_pred{i})+1])
        ylim([0.95 1.05])
        grid
        
        
        suptitle(['Prediction Verification Case ' num2str(i)])    
        saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'Case' num2str(i) '.jpg'] );
        saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'Case' num2str(i) '.fig'] );
        close
   end 