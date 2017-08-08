function [Vr_y_new, Vr_x_new, theta_new, errmin] = CalVrTheta(RigCylModel2Dobj, ILCFHydroModelobj, NonDimAmp_linear0_y, NonDimAmp_linear0_x, Vr0_y) 

    Vrn_y = RigCylModel2Dobj.FluidSpeed / ...
        (RigCylModel2Dobj.NominalNaturalFreq_y*RigCylModel2Dobj.Diameter);
    Vrn_x = RigCylModel2Dobj.FluidSpeed / ...
        (RigCylModel2Dobj.NominalNaturalFreq_x*RigCylModel2Dobj.Diameter);
    Vrn=Vrn_y;
    
    thetaList = [-180:5:180];
    yTemp = NonDimAmp_linear0_y;
    xTemp = NonDimAmp_linear0_x;
    VrTemp = Vr0_y;
    
        
    % step 3
    for i = 1: length(thetaList)
        thetaTemp = thetaList(i);
        Cm_y(i) = getDataPoint(ILCFHydroModelobj, 'Cmy', yTemp, xTemp, thetaTemp, VrTemp); 
        f = @(x) myVryfun(x, Vrn_y, RigCylModel2Dobj.MassRatio_y, ILCFHydroModelobj,  yTemp, xTemp, thetaTemp);                
        [out fval] =fsolve(f,[ VrTemp Cm_y(i)]);
        Vr_y(i) = out(1);

        Cm_x(i) = getDataPoint(ILCFHydroModelobj, 'Cmx', yTemp, xTemp, thetaTemp, VrTemp);                
        f = @(x) myVryfun(x, Vrn_x, RigCylModel2Dobj.MassRatio_x, ILCFHydroModelobj,  yTemp, xTemp, thetaTemp);                
        [out fval] =fsolve(f,[ VrTemp Cm_x(i)]);
        Vr_x(i) = out(1);
                
%         Cm_y(i) = getDataPoint(ILCFHydroModelobj, 'Cmy', yTemp, xTemp, thetaTemp, VrTemp);
%         Cm_x(i) = getDataPoint(ILCFHydroModelobj, 'Cmx', yTemp, xTemp, thetaTemp, VrTemp);
%         if RigCylModel2Dobj.MassRatio_y + Cm_y(i) >0 && RigCylModel2Dobj.MassRatio_x + Cm_x(i)>0
%             Vr_y(i) = Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + Cm_y(i))/(RigCylModel2Dobj.MassRatio_y + 1));
%             Vr_x(i) = Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + Cm_x(i))/(RigCylModel2Dobj.MassRatio_x + 1));
%         else
%             Vr_y(i) = 100; % make sure the error is larger
%             Vr_x(i) = 10;
%         end
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

    Vr_y_new = Vr_y(errInd);
    Vr_x_new = Vr_x(errInd);
    theta_new = thetaList(errInd);

    if errmin <0.1
        disp('Vr Converged')
    else
        disp('Vr Not Converged')
    end