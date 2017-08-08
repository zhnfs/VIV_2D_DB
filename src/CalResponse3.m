function [y_output, x_output, Vr_output, theta_output] = CalResponse3(RigCylModel2Dobj, ILCFHydroModelobj, NonDimAmp_linear0_y, NonDimAmp_linear0_x, Vr0_y)


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
    while m <30
        
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
        
        [errmin errInd] = min(err);

        Vr_new(m) = Vr_y(errInd);
        theta_new(m) = thetaList(errInd);
        
        if errmin <0.1
            disp('Vr Converged')
        else
            disp('Vr Not Converged')
        end
        
        if m>1
            if abs(Vr_new(m)-Vr_new(m-1))<0.5 && abs(theta_new(m)-theta_new(m-1))<15
                disp('Vr, Theta and Amp Converged')
                disp(['Vrn: ' num2str(Vrn)])
                disp(['Vr: ' num2str(Vr_new(m))])
                disp(['Theta: ' num2str(theta_new(m))])
                disp(['AmpY: ' num2str(y_new(m-1))])
                disp(['AmyX: ' num2str(x_new(m-1))])
                break
            end
        end
        
        % update Vr and theta
        VrTemp =Vr_new(m);
        thetaTemp = theta_new(m);

        % step 4
        stepSize = 0.05;
        yAmpList = [0:stepSize:1.5];
        xAmpList = [0:stepSize:0.75];
        
        for i = 1:length(yAmpList)
            yTemp = yAmpList(i);
            for j =1:length(xAmpList)
            
                xTemp = xAmpList(j);
                
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
        
        
        [C h ] = contour(xAmpList, yAmpList, errAmp);
        clabel(C,h,'FontSize',11,'fontweight','b');
        [errAmpmin] = min(min(errAmp));
        
        [errAmpyInd,errAmpxInd] = find(errAmp==min(errAmp(:)));
        
        y_new(m) = yAmpList(errAmpyInd);
        x_new(m) = yAmpList(errAmpxInd);
        
        % update Amp x and y
        
        yTemp =y_new(m);
        xTemp = x_new(m);
        
        if errAmpmin <0.2
            disp('Amp Converged')
        else
            disp('Amp Not Converged')
        end               
        m=m+1
    end
     
    y_output = y_new(end);
    x_output = x_new(end);
    Vr_output = Vr_new(end);
    theta_output = theta_new(end);
