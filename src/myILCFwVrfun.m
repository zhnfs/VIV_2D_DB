function F = myILCFwVrfun(x, RigCylModel2Dobj, ILCFHydroModelobj, Vrn_y, Vrn_x, Vr)




    F = [   
         x(1) - x(6)*Vr*Vrn_y/(4*pi^3*(RigCylModel2Dobj.MassRatio_y+x(4))*x(8));
         x(2) - x(7)*1/2*Vr*Vrn_x/(4*pi^3*(RigCylModel2Dobj.MassRatio_x+x(5))*x(9));
         Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + x(4))/(RigCylModel2Dobj.MassRatio_y + 1)) ...
         - 2*Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + x(5))/(RigCylModel2Dobj.MassRatio_x + 1));
         x(4) - getDataPoint(ILCFHydroModelobj, 'Cmy', x(1), x(2), x(3), Vr);
         x(5) - getDataPoint(ILCFHydroModelobj, 'Cmx', x(1), x(2), x(3), Vr);
         x(6) - getDataPoint(ILCFHydroModelobj, 'CLv', x(1), x(2), x(3), Vr);
         x(7) - getDataPoint(ILCFHydroModelobj, 'CDv', x(1), x(2), x(3), Vr);
         x(8) - RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ x(4))*...
            RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y);
         x(9) - RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ x(5))*...
            RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x);   
            ];

        
%     F = [   
%          x(1) - x(8)*x(4)*Vrn_y/(4*pi^3*(RigCylModel2Dobj.MassRatio_y+x(6))*x(10));
%          x(2) - x(9)*x(5)*Vrn_x/(4*pi^3*(RigCylModel2Dobj.MassRatio_x+x(7))*x(11));
%          x(4) - Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + x(6))/(RigCylModel2Dobj.MassRatio_y + 1));
%          x(5) - Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + x(7))/(RigCylModel2Dobj.MassRatio_x + 1));
%          x(4) - 2*x(5);
%          x(6) - getDataPoint(ILCFHydroModelobj, 'Cmy', x(1), x(2), x(3), x(4));
%          x(7) - getDataPoint(ILCFHydroModelobj, 'Cmx', x(1), x(2), x(3), x(4));
%          x(8) - getDataPoint(ILCFHydroModelobj, 'CLv', x(1), x(2), x(3), x(4));
%          x(9) - getDataPoint(ILCFHydroModelobj, 'CDv', x(1), x(2), x(3), x(4));
%          x(10) - RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ x(6))*...
%             RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y);
%          x(11) - RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ x(7))*...
%             RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x);   
%             ];
 
end