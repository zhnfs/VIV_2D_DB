function plot_Coef(CLv_pred, CDv_pred, Cmy_pred, Cmx_pred,Damp_pred_y, Damp_pred_x, vrn_pred,...
   vr_pred, Yad_pred, Xad_pred, Theta_pred , ...
   CLv_free, CDv_free, Cmy_free, Cmx_free, vr_fr, Yad_fr, Xad_fr, Theta_fr,RigCylModel2Dobj, ILCFHydroModelobj,freqRatio, freqRatioInd, compIndx,  Vrn_y, Vrn_x, outputDir, marker)
% 
%     for i = compIndx
%         CLv_fr_pred(i) = getDataPoint(ILCFHydroModelobj, 'CLv', Yad_fr(freqRatioInd, i), Xad_fr(freqRatioInd, i), Theta_fr(freqRatioInd, i)./180, vr_fr(freqRatioInd, i));
%         CDv_fr_pred(i) = getDataPoint(ILCFHydroModelobj, 'CDv', Yad_fr(freqRatioInd, i), Xad_fr(freqRatioInd, i), Theta_fr(freqRatioInd, i)./180, vr_fr(freqRatioInd, i));
%         Cmy_fr_pred(i) = getDataPoint(ILCFHydroModelobj, 'Cmy', Yad_fr(freqRatioInd, i), Xad_fr(freqRatioInd, i), Theta_fr(freqRatioInd, i)./180, vr_fr(freqRatioInd, i));
%         Cmx_fr_pred(i) = getDataPoint(ILCFHydroModelobj, 'Cmx', Yad_fr(freqRatioInd, i), Xad_fr(freqRatioInd, i), Theta_fr(freqRatioInd, i)./180, vr_fr(freqRatioInd, i));
%     end


    for i = compIndx
        RigCylModel2Dobj(freqRatioInd,i).MassRatio_y = 3.5;
        RigCylModel2Dobj(freqRatioInd,i).MassRatio_x = 3.5;
        Cmy_fr_pred(i) = (vr_fr(freqRatioInd, i)/Vrn_y(freqRatioInd, i))^2*(RigCylModel2Dobj(freqRatioInd,i).MassRatio_y + 1) - RigCylModel2Dobj(freqRatioInd,i).MassRatio_y;
        Cmx_fr_pred(i) = (vr_fr(freqRatioInd, i)/2/Vrn_x(freqRatioInd, i))^2*(RigCylModel2Dobj(freqRatioInd,i).MassRatio_x + 1) - RigCylModel2Dobj(freqRatioInd,i).MassRatio_x;
        
        dampY_fr_pred(i) = (RigCylModel2Dobj(freqRatioInd,i).DamperCons_y/2/((RigCylModel2Dobj(freqRatioInd,i).MassRatio_y+ Cmy_fr_pred(i))*...
            RigCylModel2Dobj(freqRatioInd,i).FluidDensity*pi*RigCylModel2Dobj(freqRatioInd,i).Diameter^2*RigCylModel2Dobj(freqRatioInd,i).Length/4*2*pi*RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_y));
        dampX_fr_pred(i) = (RigCylModel2Dobj(freqRatioInd,i).DamperCons_x/2/((RigCylModel2Dobj(freqRatioInd,i).MassRatio_x+ Cmx_fr_pred(i))*...
            RigCylModel2Dobj(freqRatioInd,i).FluidDensity*pi*RigCylModel2Dobj(freqRatioInd,i).Diameter^2*RigCylModel2Dobj(freqRatioInd,i).Length/4*2*pi*RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_x));
        
        CLv_fr_pred(i) = Yad_fr(freqRatioInd, i)*(4*pi^3*(RigCylModel2Dobj(freqRatioInd,i).MassRatio_y+Cmy_fr_pred(i))*dampY_fr_pred(i))./Vrn_y(freqRatioInd, i)./vr_fr(freqRatioInd, i);
        CDv_fr_pred(i) = Xad_fr(freqRatioInd, i)*(4*pi^3*(RigCylModel2Dobj(freqRatioInd,i).MassRatio_x+Cmx_fr_pred(i))*dampX_fr_pred(i))./Vrn_x(freqRatioInd, i)./(vr_fr(freqRatioInd, i)/2);

    end
    
    
    figure
    subplot(1,2,1)
    for i=compIndx
        for j = 1:length(CLv_pred{i})
%             plot(vrn_pred(freqRatioInd,i),CLv_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),Damp_pred_y{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),dampY_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('CLv','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
%     axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    subplot(1,2,2)
    for i=compIndx
        for j = 1:length(CLv_pred{i})
%             plot(vrn_pred(freqRatioInd,i),CLv_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),Damp_pred_x{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),dampX_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('CLv','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
%     axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    legend({'Predicted from Forced Vib', 'Derivied from Free Vib Motion'},'Location','Southwest')
    suptitle('Damping Coefficients')    
    saveas(gcf,[outputDir filesep 'Damp' num2str(freqRatio(freqRatioInd)) marker '.png'] );
    saveas(gcf,[outputDir filesep 'Damp' num2str(freqRatio(freqRatioInd)) marker '.fig'] );
    
    
    
% Added mass verification
    
    figure; 
    subplot(2,2,2)
    for i=compIndx
        for j = 1:length(CLv_pred{i})
%             plot(vrn_pred(freqRatioInd,i),CLv_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),CLv_pred{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),CLv_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('CLv','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
%     axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    
    subplot(2,2,4)
    for i=compIndx
        for j = 1:length(CDv_pred{i})
%             plot(vrn_pred(freqRatioInd,i),-CDv_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),CDv_pred{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),CDv_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('CDv','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
%     axis([4 10 -180 180]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid

    
     subplot(2,2,1)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
%             plot(vrn_pred(freqRatioInd,i),Cmy_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),Cmy_pred{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),Cmy_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
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
    legend({'Predicted from Forced Vib', 'Observed Free Vib'},'Location','Southwest')
%     legend({'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')

    
    subplot(2,2,3)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
%             plot(vrn_pred(freqRatioInd,i),-Cmx_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),Cmx_pred{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),Cmx_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
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
    
    
    suptitle('Hydrodynamic Coefficients')    
    saveas(gcf,[outputDir filesep 'Coef' num2str(freqRatio(freqRatioInd)) marker '.png'] );
    saveas(gcf,[outputDir filesep 'Coef' num2str(freqRatio(freqRatioInd)) marker '.fig'] );
    
    
    % Added mass verification
    
    figure; 
    subplot(2,2,2)
    for i=compIndx
        for j = 1:length(CLv_pred{i})
            plot(vrn_pred(freqRatioInd,i),CLv_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),CLv_pred{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),CLv_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('CLv','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
%     axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    
    subplot(2,2,4)
    for i=compIndx
        for j = 1:length(CDv_pred{i})
            plot(vrn_pred(freqRatioInd,i),-CDv_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),CDv_pred{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),CDv_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('CDv','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
%     axis([4 10 -180 180]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid

    
     subplot(2,2,1)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmy_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),Cmy_pred{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),Cmy_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
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
%     legend({'Predicted from Forced Vib', 'Observed Free Vib'},'Location','Southwest')
    legend({'Calculated from Free Vib Force Measurement', 'Predicted from Forced Vib', 'Derivied from Free Vib Motion'},'Location','Southwest')

    
    subplot(2,2,3)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),-Cmx_free(i), 'Color','g','Marker','d','MarkerSize',9);
            hold on
            plot(vrn_pred(freqRatioInd,i),Cmx_pred{i}(j),'Color', 'm','Marker','<','MarkerSize',9);          
            plot(vrn_pred(freqRatioInd,i),Cmx_fr_pred(i),'Color', 'b','Marker','o','MarkerSize',9);
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
    
    
    suptitle('Hydrodynamic Coefficients')    
    saveas(gcf,[outputDir filesep 'Coef' num2str(freqRatio(freqRatioInd)) marker '2.png'] );
    saveas(gcf,[outputDir filesep 'Coef' num2str(freqRatio(freqRatioInd)) marker '2.fig'] );
        
        