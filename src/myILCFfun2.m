function F = myILCFfun2(x, RigCylModel2Dobj, ILCFHydroModelobj, Vrn_y, Vrn_x, targetFreqRatio)

     F = [   
         x(1) / (x(8)*x(4)*Vrn_y/(4*pi^3*(RigCylModel2Dobj.MassRatio_y+x(6))*x(10)))-1;
         x(2) / (x(9)*x(5)*Vrn_x/(4*pi^3*(RigCylModel2Dobj.MassRatio_x+x(7))*x(11)))-1;
         x(4) / (Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + x(6))/(RigCylModel2Dobj.MassRatio_y + 1)))-1;
         x(5) / (Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + x(7))/(RigCylModel2Dobj.MassRatio_x + 1)))-1;
%          x(4) / (2*x(5))-1;
         x(4) / (targetFreqRatio*x(5))-1;
         x(6) / getDataPoint(ILCFHydroModelobj, 'Cmy', x(1), x(2), x(3), x(4))-1;
         x(7) / getDataPoint(ILCFHydroModelobj, 'Cmx', x(1), x(2), x(3), x(4))-1;
         x(8) / getDataPoint(ILCFHydroModelobj, 'CLv', x(1), x(2), x(3), x(4))-1;
         x(9) / getDataPoint(ILCFHydroModelobj, 'CDv', x(1), x(2), x(3), x(4))-1;
         x(10) / (RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ x(6))*...
            RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y))-1;
         x(11) / (RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ x(7))*...
            RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x))-1;   
            ];
end