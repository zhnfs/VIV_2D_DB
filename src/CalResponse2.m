function [y_final, x_final, Vr_y_final, Vr_x_final] = CalResponse2(RigCylModel2Dobj, ILCFHydroModelobj, Amp0_y, Amp0_x, Vr0, Theta0)


%% solution 1 using iterations

    Vrn_y = RigCylModel2Dobj.FluidSpeed / ...
        (RigCylModel2Dobj.NominalNaturalFreq_y*RigCylModel2Dobj.Diameter);
    Vrn_x = RigCylModel2Dobj.FluidSpeed / ...
        (RigCylModel2Dobj.NominalNaturalFreq_x*RigCylModel2Dobj.Diameter);
    Vrn=Vrn_y;
    
    yTemp = Amp0_y;
    xTemp = Amp0_x;
    VrTemp = Vr0;
    thetaTemp = Theta0;
    
    n=1;
    while(n<2)

        Cm_y(n) = getDataPoint(ILCFHydroModelobj, 'Cmy', yTemp, xTemp, thetaTemp, VrTemp);
        Cm_x(n) = getDataPoint(ILCFHydroModelobj, 'Cmx', yTemp, xTemp, thetaTemp, VrTemp);
        Vr_y_new(n) = Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + Cm_y(n))/(RigCylModel2Dobj.MassRatio_y + 1));
        Vr_x_new(n) = Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + Cm_x(n))/(RigCylModel2Dobj.MassRatio_x + 1));
        err(n) = abs( Vr_y_new(n) - 2* Vr_x_new(n));

        DampingRatio_y(n) = RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ Cm_y(n) )*...
            RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y);      

        CLv_y(n) = getDataPoint(ILCFHydroModelobj, 'CLv', yTemp, xTemp, thetaTemp, VrTemp);

        y_new(n) = CLv_y /(4*pi^3*(RigCylModel2Dobj.MassRatio_y+Cm_y(n))*DampingRatio_y);

        DampingRatio_x(n) = RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ Cm_x(n) )*...
            RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x);      

        CLv_x(n) = getDataPoint(ILCFHydroModelobj, 'CLv', yTemp, xTemp, thetaTemp, VrTemp);

        x_new(n) = CLv_x(n) /(4*pi^3*(RigCylModel2Dobj.MassRatio_x+Cm_x(n))*DampingRatio_x(n));

        if y_new(n) <0 || x_new(n) <0  
            y_new(n) = yTemp+0.1;
            x_new(n) = xTemp+0.1;
        end
        
        if ~isreal(Vr_y_new(n))||~isreal(Vr_x_new(n))
            Vr_y_new(n) = VrTemp; % make sure the error is larger
            Vr_x_new(n) = VrTemp/2;
            err(n)  = NaN;
        end
        
        if abs(yTemp - y_new(n)) < 0.1 && abs(xTemp - x_new(n)) < 0.1
            disp('Amp Converged')
            break
        end
        
        % update 
        VrTemp =Vr_y_new(n);
%         thetaTemp = theta_new(n);
        yTemp = y_new(n);
        xTemp = x_new(n);
            
        n=n+1
    end
    
    y_final = y_new(end);
    x_final = x_new(end);
    Vr_y_final = Vr_y_new(end);
    Vr_x_final = Vr_x_new(end);
%     theta_final = theta_new(end);