function F = myILCFVrfun(x, RigCylModel2Dobj, ILCFHydroModelobj, Vrn_y, Vrn_x, targetFreqRatio)

    F = [   

         x(1) - (Vrn_y* sqrt((RigCylModel2Dobj.MassRatio_y + x(3))/(RigCylModel2Dobj.MassRatio_y + 1)));
         x(2) - (Vrn_x* sqrt((RigCylModel2Dobj.MassRatio_x + x(4))/(RigCylModel2Dobj.MassRatio_x + 1)));
         x(1) - (targetFreqRatio*x(2));
         x(3) - getDataPoint(ILCFHydroModelobj, 'Cmy', 1, 0.3, x(5), x(1));
         x(4) - getDataPoint(ILCFHydroModelobj, 'Cmx', 1, 0.3, x(5), x(1));
            ];

end