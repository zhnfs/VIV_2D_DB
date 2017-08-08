function [y_new, x_new, Vr_new, theta_new] = CalResponse(RigCylModel2Dobj, ILCFHydroModelobj, NonDimAmp_linear0_y, NonDimAmp_linear0_x, Vr0_y)


%% solution 1 using iterations

%     Vrn = RigCylModelobj.FluidSpeed / ...
%         (2*RigCylModelobj.NominalNaturalFreq*RigCylModelobj.Diameter);
%     Vr = Vrn;
% %     CaIL = feval(CaIL_sfitobj.surface_sfit,1/Vr, 0 );
%     CaIL_A0point = getCaILPoint(linearCaIlobj, Vr);
%     VrTemp = Vrn * sqrt((RigCylModelobj.MassRatio + CaIL_A0point )/...
%         (RigCylModelobj.MassRatio+1));
%     iter =0;
%     while(abs(Vr - VrTemp)>0.01)
%         Vr = VrTemp;
% %         CaIL = feval(CaIL_sfitobj.surface_sfit,1/Vr, 0 );
%         CaIL_A0point = getCaILPoint(linearCaIlobj, Vr);
%         VrTemp = Vrn * sqrt((RigCylModelobj.MassRatio + CaIL_A0point )/...
%         (RigCylModelobj.MassRatio+1));
%         iter = iter +1;
%         if iter == 10000
%             Vr = Vrn;
%             break;
%         end
%         
%     end
    
    %% solution 2 using partial interation   
    

    Vrn_y = RigCylModel2Dobj.FluidSpeed / ...
        (RigCylModel2Dobj.NominalNaturalFreq_y*RigCylModel2Dobj.Diameter);
    Vrn_x = RigCylModel2Dobj.FluidSpeed / ...
        (RigCylModel2Dobj.NominalNaturalFreq_x*RigCylModel2Dobj.Diameter);
    Vrn=Vrn_y;
    
    thetaList = [-180:5:180];
    yTemp = NonDimAmp_linear0_y;
    xTemp = NonDimAmp_linear0_x;
    VrTemp = Vr0_y;
    
    m = 1;
    while m <100
        
        % step 3
        for i = 1: length(thetaList)
            thetaTemp = thetaList(i);
            Cm_y(i) = getDataPoint(ILCFHydroModelobj, 'Cmy', yTemp, xTemp, thetaTemp, VrTemp);
            Cm_x(i) = getDataPoint(ILCFHydroModelobj, 'Cmx', yTemp, xTemp, thetaTemp, VrTemp);
            if RigCylModel2Dobj.MassRatio_y + Cm_y(i) >0 && RigCylModel2Dobj.MassRatio_x + Cm_x(i)>0
                Vr_y(i) = Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + Cm_y(i))/(RigCylModel2Dobj.MassRatio_y + 1));
                Vr_x(i) = Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + Cm_x(i))/(RigCylModel2Dobj.MassRatio_x + 1));
            else
                Vr_y(i) = 100; % make sure the error is larger
                Vr_x(i) = 10;
            end
            err(i) = abs( Vr_y(i) - 2* Vr_x(i));
        end

%         figure
%         subplot(2,2,1)
%         plot(thetaList,Cm_y,'o-')
%         set(gca,'XTick',-180:45:180)
%         set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
%         title('Cm_y')
%         subplot(2,2,2)
%         plot(thetaList,Cm_x,'o-') 
%         set(gca,'XTick',-180:45:180)
%         set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
%         title('Cm_x')
%         subplot(2,2,3)
%         plot(thetaList,Vr_y,'o-')
%         hold on
%         plot(thetaList,2.*Vr_x,'r*-')
%         legend('Vry','2*Vrx')
%         set(gca,'XTick',-180:45:180)
%         set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
%         title('Vr')
%         subplot(2,2,4)
%         plot(thetaList,err,'o-')
%         set(gca,'XTick',-180:45:180)
%         set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
%         title('Vr Error')
%         close
        
        [errmin errInd] = min(err);

        Vr_new = Vr_y(errInd);
        theta_new = thetaList(errInd);
        
        if abs(Vr_new-VrTemp) < 0.1 && abs(theta_new-thetaTemp) < 45
            disp('Vr and theta Converged')
            break
        end
        
        % update Vr and theta
        VrTemp =Vr_new;
        thetaTemp = theta_new;

        % step 4
        n=1;
        stepSize=2;
        while (n<100)
            Cm_y = getDataPoint(ILCFHydroModelobj, 'Cmy', yTemp, xTemp, thetaTemp, VrTemp);

            DampingRatio_y = RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ Cm_y )*...
                RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y);      

            CLv_y = getDataPoint(ILCFHydroModelobj, 'CLv', yTemp, xTemp, thetaTemp, VrTemp);

            y_new = CLv_y /(4*pi^3*(RigCylModel2Dobj.MassRatio_y+Cm_y)*DampingRatio_y);


            
            Cm_x = getDataPoint(ILCFHydroModelobj, 'Cmx', yTemp, xTemp, thetaTemp, VrTemp);

            DampingRatio_x = RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ Cm_x )*...
                RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x);      

            CLv_x = getDataPoint(ILCFHydroModelobj, 'CLv', yTemp, xTemp, thetaTemp, VrTemp);

            x_new = CLv_x /(4*pi^3*(RigCylModel2Dobj.MassRatio_x+Cm_x)*DampingRatio_x);

%             if y_new <0 || x_new <0  
%                 y_new = yTemp+0.1;
%                 x_new = xTemp+0.1;
%             end
            
            if abs(yTemp - y_new) < 0.1 && abs(xTemp - x_new) < 0.1
                disp('Amp Converged')
                break
            end

            % update Y and X

            if yTemp+y_new/stepSize <0 ||xTemp+x_new/stepSize<0  
%                 yTemp = 0;
%                 xTemp = 0;
                stepSize=stepSize*2;
            else
                yTemp = yTemp+y_new/stepSize;
                xTemp = yTemp+y_new/stepSize;
            end
            n=n+1
        end
        
        m=m+1

    end
        
end