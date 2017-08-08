function VrResult = CalVrX(RigCylModel2Dobj)


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
    
    %% solution 2 using intersections
    
    Vrn = RigCylModel2Dobj.FluidSpeed / ...
        (RigCylModel2Dobj.NominalNaturalFreq_x*RigCylModel2Dobj.Diameter);

    Increment = 0.01;
    Vr= .1:Increment:10;
    VrTemp =  zeros(size(Vr));
    
    for i = 1: length(Vr)
         CaIL_A0point = getCaPoint(linearCaIlobj, Vr(i));
         VrTemp(i) = Vrn * sqrt((RigCylModel2Dobj.MassRatio_x + CaIL_A0point )/...
        (RigCylModel2Dobj.MassRatio_x+1));
    
    end
    
%     plot(Vr,VrTemp)
%     hold on
%     plot(Vr,Vr)
    
    [VrResult, VrTempResult_1] = intersections(Vr, Vr, Vr, VrTemp,1 );
    
    if(isempty(VrResult))
        VrResult = Vrn;
    end
    VrResult = VrResult(1);
end