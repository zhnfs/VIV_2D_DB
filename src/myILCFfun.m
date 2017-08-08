function F = myILCFfun(x, RigCylModel2Dobj, ILCFHydroModelobj, Vrn_y, Vrn_x, targetFreqRatio)

%     Diameter = 0.0762; % m  
%     Riserlength = 2; % m
%     fluidDensity = 1000; %kg/m^3 
% 
%     fn_y = [0.715	0.799	0.894	0.977	0.698	0.704];
%     fn_x = [0.715	0.974   1.226   1.488   1.167   1.336];
% 
%     springCons_y = [697  902 1118 1386 967 1013];
%     springCons_x = [609  1301 2130 2850 2580 3235];
% 
%     displacedMass = fluidDensity * pi*(Diameter/2)^2*Riserlength;
% 
%     massRatio_y = (springCons_y./(fn_y*2*pi).^2)/displacedMass-1;
%     massRatio_x = (springCons_x./(fn_x*2*pi).^2)/displacedMass-1;
% 
%     strucDampRatio_y =0.01*[2.2	1.3	1.1	1.6	2.6	6.2];
%     strucDampRatio_x =0.01*[2.2	1.7	2.5	3.2	2.9	2.5];
% 
%     massDensity_x = fluidDensity.*massRatio_x;
%     massDensity_y = fluidDensity.*massRatio_y;
% 
% 
%     damperCons_x = strucDampRatio_x*2.*(1+massRatio_x)* displacedMass*2*pi.*fn_x;
%     damperCons_y = strucDampRatio_y*2.*(1+massRatio_y)* displacedMass*2*pi.*fn_y;
% 
%     FluidVel = [0.161,0.174,0.188,0.201,0.215,0.228,0.241,0.255,0.268,0.282,...
%         0.295,0.308,0.322,0.335,0.349,0.362,0.376,0.389,0.402,0.416,0.429,...
%         0.443,0.456,0.469,0.483,0.496,0.510,0.523,0.536,0.550,0.563,0.577,0.590,0.604,0.617,0.630,0.644];
%     freqRatioInd=6;
%     i=9;
%     
%     fluidSpeed = FluidVel(i);  
%     RigCylModel2Dobj = RigCylModel2D(Diameter, Riserlength, massDensity_x(freqRatioInd), massDensity_y(freqRatioInd), ...
%         fluidSpeed, fluidDensity, springCons_x(freqRatioInd), springCons_y(freqRatioInd), damperCons_x(freqRatioInd), damperCons_y(freqRatioInd));
% %     RigCylModel2Dobj.PrintRigCylModel2D()    
%  
%     Vrn_y = RigCylModel2Dobj.FluidSpeed / ...
%         (RigCylModel2Dobj.NominalNaturalFreq_y*RigCylModel2Dobj.Diameter);
%     Vrn_x = RigCylModel2Dobj.FluidSpeed / ...
%         (RigCylModel2Dobj.NominalNaturalFreq_x*RigCylModel2Dobj.Diameter);
% 
%     
%     ILCFHydroModelobj = ILCFHydroModel();


    F = [   
         x(1) - (x(8)*x(4)*Vrn_y/(4*pi^3*(RigCylModel2Dobj.MassRatio_y+x(6))*x(10)));
         x(2) - (x(9)*x(5)*Vrn_x/(4*pi^3*(RigCylModel2Dobj.MassRatio_x+x(7))*x(11)));
         x(4) - (Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + x(6))/(RigCylModel2Dobj.MassRatio_y + 1)));
         x(5) - (Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + x(7))/(RigCylModel2Dobj.MassRatio_x + 1)));
%          x(4) - (2*x(5));
         x(4) - (targetFreqRatio*x(5));
         x(6) - getDataPoint(ILCFHydroModelobj, 'Cmy', x(1), x(2), x(3), x(4));
         x(7) - getDataPoint(ILCFHydroModelobj, 'Cmx', x(1), x(2), x(3), x(4));
         x(8) - getDataPoint(ILCFHydroModelobj, 'CLv', x(1), x(2), x(3), x(4));
         x(9) - getDataPoint(ILCFHydroModelobj, 'CDv', x(1), x(2), x(3), x(4));
         x(10) - (RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ x(6))*...
            RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y));
         x(11) - (RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ x(7))*...
            RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x));   
            ];
        
%             F = [   
%          x(1) / (x(8)*x(4)*Vrn_y/(4*pi^3*(RigCylModel2Dobj.MassRatio_y+x(6))*x(10)))-1;
%          x(2) / (x(9)*x(5)*Vrn_x/(4*pi^3*(RigCylModel2Dobj.MassRatio_x+x(7))*x(11)))-1;
%          x(4) / (Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + x(6))/(RigCylModel2Dobj.MassRatio_y + 1)))-1;
%          x(5) / (Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + x(7))/(RigCylModel2Dobj.MassRatio_x + 1)))-1;
%          x(4) / (2*x(5))-1;
%          x(6) / getDataPoint(ILCFHydroModelobj, 'Cmy', x(1), x(2), x(3), x(4))-1;
%          x(7) / getDataPoint(ILCFHydroModelobj, 'Cmx', x(1), x(2), x(3), x(4))-1;
%          x(8) / getDataPoint(ILCFHydroModelobj, 'CLv', x(1), x(2), x(3), x(4))-1;
%          x(9) / getDataPoint(ILCFHydroModelobj, 'CDv', x(1), x(2), x(3), x(4))-1;
%          x(10) / (RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ x(6))*...
%             RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y))-1;
%          x(11) / (RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ x(7))*...
%             RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x))-1;   
%             ];
%      
%     Cm_y = getDataPoint(ILCFHydroModelobj, 'Cmy', x(1), x(2), x(3), x(4));
%     Cm_x = getDataPoint(ILCFHydroModelobj, 'Cmx', x(1), x(2), x(3), x(4));
% 
%     CLv_y = getDataPoint(ILCFHydroModelobj, 'CLv', x(1), x(2), x(3), x(4));
%     CLv_x = getDataPoint(ILCFHydroModelobj, 'CDv', x(1), x(2), x(3), x(4));
%     
%     DampingRatio_y = RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ Cm_y )*...
%         RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y);
% 
%     DampingRatio_x = RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ Cm_x )*...
%         RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x);      
%     
%     F = [x(4) - Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + Cm_y)/(RigCylModel2Dobj.MassRatio_y + 1));
%          x(5) - Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + Cm_x)/(RigCylModel2Dobj.MassRatio_x + 1));
%          x(1) - CLv_y /(4*pi^3*(RigCylModel2Dobj.MassRatio_y+Cm_y)*DampingRatio_y);
%          x(2) - CLv_x /(4*pi^3*(RigCylModel2Dobj.MassRatio_x+Cm_x)*DampingRatio_x);
%          x(4) - x(5);];
% 
end