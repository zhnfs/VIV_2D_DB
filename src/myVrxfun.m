function F = myVrxfun(x, Vrn, MassRatio, ILCFHydroModelobj, Ampy, Amyx, theta)

    F = [   
         x(1) - Vrn* sqrt((MassRatio + x(2))/(MassRatio + 1));
         x(2) - getDataPoint(ILCFHydroModelobj, 'Cmx', Ampy, Amyx, theta, x(1));

         ];
end