function [y_output, x_output, Vr_output, theta_output] = CalResponse4(RigCylModel2Dobj, ILCFHydroModelobj, NonDimAmp_linear0_y, NonDimAmp_linear0_x, Vr0_y)


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
    
    stepSize = 0.05;
    yAmpList = [0:stepSize:1.5];
    xAmpList = [0:stepSize:0.75];

    for i = 1:length(yAmpList)
        yTemp = yAmpList(i);       
        for j =1:length(xAmpList)
            xTemp = xAmpList(j);                
            % step 3
            for m = 1: length(thetaList)
                thetaTemp = thetaList(m);
                
                Cm_y(m) = getDataPoint(ILCFHydroModelobj, 'Cmy', yTemp, xTemp, thetaTemp, VrTemp); 
                f = @(x) myVryfun(x, Vrn_y, RigCylModel2Dobj.MassRatio_y, ILCFHydroModelobj,  yTemp, xTemp, thetaTemp);                
                [out fval] =fsolve(f,[ VrTemp Cm_y(m)]);
                Vr_y(m) = out(1);
                
                Cm_x(m) = getDataPoint(ILCFHydroModelobj, 'Cmx', yTemp, xTemp, thetaTemp, VrTemp);                
                f = @(x) myVryfun(x, Vrn_x, RigCylModel2Dobj.MassRatio_x, ILCFHydroModelobj,  yTemp, xTemp, thetaTemp);                
                [out fval] =fsolve(f,[ VrTemp Cm_x(m)]);
                Vr_x(m) = out(1);
                
%                 Cm_y(m) = getDataPoint(ILCFHydroModelobj, 'Cmy', yTemp, xTemp, thetaTemp, VrTemp);
%                 Cm_x(m) = getDataPoint(ILCFHydroModelobj, 'Cmx', yTemp, xTemp, thetaTemp, VrTemp);
%                 if RigCylModel2Dobj.MassRatio_y + Cm_y(m) >0 && RigCylModel2Dobj.MassRatio_x + Cm_x(m)>0
%                     Vr_y(m) = Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + Cm_y(m))/(RigCylModel2Dobj.MassRatio_y + 1));
%                     Vr_x(m) = Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + Cm_x(m))/(RigCylModel2Dobj.MassRatio_x + 1));
%                 else
%                     Vr_y(m) = 100; % make sure the error is larger
%                     Vr_x(m) = 10;
%                 end
                errVr(m) = abs( Vr_y(m) - 2* Vr_x(m));
            end
            [errmin errInd] = min(errVr);

            VrTemp= Vr_y(i,j,errInd);
            thetaTemp = thetaList(errInd);
            VrList(i,j) = VrTemp;
            thetaList(i,j) = thetaTemp;
            
            % step 4
            Cm_y = getDataPoint(ILCFHydroModelobj, 'Cmy', yTemp, xTemp, thetaTemp, VrTemp);

            DampingRatio_y = RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ Cm_y )*...
                RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y);      

            CLv_y = getDataPoint(ILCFHydroModelobj, 'CLv', yTemp, xTemp, thetaTemp, VrTemp);

            y_new_temp = CLv_y /(4*pi^3*(RigCylModel2Dobj.MassRatio_y+Cm_y)*DampingRatio_y);

            Cm_x = getDataPoint(ILCFHydroModelobj, 'Cmx', yTemp, xTemp, thetaTemp, VrTemp);

            DampingRatio_x = RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ Cm_x )*...
                RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x);      

            CLv_x = getDataPoint(ILCFHydroModelobj, 'CDv', yTemp, xTemp, thetaTemp, VrTemp);

            x_new_temp = CLv_x /(4*pi^3*(RigCylModel2Dobj.MassRatio_x+Cm_x)*DampingRatio_x);

            errAmp(i,j) = abs(yTemp - y_new_temp) + abs(xTemp - x_new_temp);

        end
    end
    
            
%         [C h ] = contour(xAmpList, yAmpList, errAmp);
%         clabel(C,h,'FontSize',11,'fontweight','b');
    [errAmpmin] = min(min(errAmp));

    [errAmpyInd,errAmpxInd] = find(errAmp==min(errAmp(:)));

    y_output = yAmpList(errAmpyInd);
    x_output = yAmpList(errAmpxInd);

    Vr_output = VrList(errAmpyInd, errAmpxInd);
    theta_output = thetaList(errAmpyInd, errAmpxInd);
    

    if errAmpmin <0.2
        disp('Amp Converged')
    else
        disp('Amp Not Converged')
    end  

