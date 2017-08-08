% function ILCFRigidMasterScript()

% clear all
close all
% createCeList();

col = distinguishable_colors(20);
% col = col([13,15,3,8,1,10,14,7,4],:);
shp = ['o';'s';'p';'v';'*';'d';'^';'>';'<';'.'];

utility = Utility.getInstance();
filesep = utility.getFilesep();
[vivaRoot, optRoot, experimentRoot, fatigueRoot, utilityRoot, InlineVIVRoot, CrossflowVIVRoot, mooringRoot, ILCFVIVRoot] = getRoots();

outputDir =[ILCFVIVRoot filesep 'output' filesep 'ICFLMotionHZ1e-5Test'];
utility.conditionalMkdir(outputDir);

%% Build Hydrodynamics Model

% % Inline
CeILobj = CeIL_sfit;
CaILobj = CaIL_sfit;
% CeILobj.plotCeIL_sfit();
% CaILobj.plotCaIL_sfit();
% 
gridside = 0.005;
nonDimFreq =0.2:gridside:0.9;
nonDimAmp = 0:gridside:0.3;

linearCeILobj = LinearCeILModel(CeILobj, nonDimFreq, nonDimAmp);
linearCaILobj = LinearCaILModel(CaILobj, nonDimFreq);
% plotLinearCeIL3DModel(linearCeILobj, CeILobj);
% plotLinearCaIL3DModel(linearCaILobj);
% plotLinearCeCaIL2D(linearCeILobj, linearCaILobj);
%  
% % Crossflow
basic_bare =  load([InlineVIVRoot filesep 'DATA' filesep 'basic_bare_opt']);

gridside = 0.05;
nonDimFreq =basic_bare(:,1);
nonDimAmp = 0:gridside:1.2;

linearCaCFobj = LinearCaCFModel(basic_bare, nonDimFreq);
linearCeCFobj = LinearCeCFModel(basic_bare, nonDimFreq, nonDimAmp);
% plotLinearCeCF3DModel(linearCeCFobj);
% plotLinearCaCF3DModel(linearCaCFobj);
% plotLinearCeCaCF2D(linearCeCFobj, linearCaCFobj);


% % IL-CF 

ILCFHydroModelobj = ILCFHydroModel();


%% Free Virbation Experiment setup 

% %% force vibration
% cylD = 0.0381; % m
% span = 0.6858; % m
% mass = 0.727; % kg
% lengthtosensor = 0.9271; %m
% momarm = lengthtosensor - span/2; % m
% cofm = 0.2318; % m , center of mass below the force sensor
% %U = 0.246; %m/s
% mass = 0.727; % kg
% massRatio = mass/(fluidDensity*Riserlength*pi/4*Diameter^2);   %0.104 lbf/ft = 1.51776591 N/m
% massDensity= fluidDensity * massRatio;  %kg/m^3

Diameter = 0.0762; % m  
Riserlength = 2; % m
aspRatio = Riserlength/Diameter;
fluidDensity = 1000; %kg/m^3 

% my is not the mass in air
% my=[34.5 35.8 35.8 36.8 50.3 51.8]; 
% mx=[30.2 34.7 33.7 32.6 48.02 45.9];
% ky = [697 902 1118 1386 967 1013];
% kx = [609  1301 2130 2850 2580 3235];
% fn_y_air = sqrt(ky./my)/2/pi
% fn_x_air = sqrt(kx./mx)/2/pi

% frequency measured in still water
freqRatio =[1	1.22	1.37	1.52	1.67	1.9];
freqRatioNames ={ '1p0', '1p22', '1p37', '1p52', '1p67', '1p9'};
fn_y = [0.715	0.799	0.894	0.977	0.698	0.704];
fn_x = [0.715	0.974   1.226   1.488   1.167   1.336];

springCons_y = [697  902 1118 1386 967 1013];
springCons_x = [609  1301 2130 2850 2580 3235];

displacedMass = fluidDensity * pi*(Diameter/2)^2*Riserlength;

massRatio_y = (springCons_y./(fn_y*2*pi).^2)/displacedMass-1;
massRatio_x = (springCons_x./(fn_x*2*pi).^2)/displacedMass-1;

% 
% massRatio_y = (springCons_y./(fn_y*2*pi).^2)/displacedMass;
% massRatio_x = (springCons_x./(fn_x*2*pi).^2)/displacedMass;

strucDampRatio_y =0.00001*[2.2	1.3	1.1	1.6	2.6	6.2];
strucDampRatio_x =0.00001*[2.2	1.7	2.5	3.2	2.9	2.5];

% massRatio_x =[3.3 3.8 3.7 3.6 5.3 5]; % this massRatio in paper contains
% normal added mass 1
% massRatio_y =[3.8 3.9 3.9 4   5.5	5.7];

massDensity_x = fluidDensity.*massRatio_x;
massDensity_y = fluidDensity.*massRatio_y;


damperCons_x = strucDampRatio_x*2.*(1+massRatio_x)* displacedMass*2*pi.*fn_x;
damperCons_y = strucDampRatio_y*2.*(1+massRatio_y)* displacedMass*2*pi.*fn_y;

%% load Free Vibration Data
% freqRatioInd is for different frequency ratio 
freeVibDataDir = ['..' filesep 'freeVibJMD' filesep];

% different freq ratio test has different fluid velocity

% FluidVel = [0.161,0.174,0.188,0.201,0.215,0.228,0.241,0.255,0.268,0.282,...
%     0.295,0.308,0.322,0.335,0.349,0.362,0.376,0.389,0.402,0.416,0.429,...
%     0.443,0.456,0.469,0.483,0.496,0.510,0.523,0.536,0.550,0.563,0.577,0.590,0.604,0.617,0.630,0.644];
% FluidVel = FluidVel(1);

% ReyNo = Diameter.*FluidVel/1e-6; 

load([freeVibDataDir 'vel.mat']);

newmatfiles = {'new1p0' 'new1p22' 'new1p37' 'new1p52' 'new1p67' 'new1p9'};
phasematfiles = {'phase1p0' 'phase1p22' 'phase1p37' 'phase1p52' 'phase1p67' 'phase1p9'};
corrmatfiles = {'corr1p0' 'corr1p22' 'corr1p37' 'corr1p52' 'corr1p67' 'corr1p9'};
ampfixmatfiles = {'ampfix1p0' 'ampfix1p22' 'ampfix1p37' 'ampfix1p52' 'ampfix1p67' 'ampfix1p9'};

% compIndx = 3:37;
% compIndx = [9:25];
% compIndx = [9 17 25];
compIndx = [21];
numRun = 100;

% % for freqRatioInd= 4:6 
for freqRatioInd= 6 
    
    eval(['FluidVel = vels' freqRatioNames{freqRatioInd} ';'])
    eval(['load ' freeVibDataDir char(newmatfiles(freqRatioInd)) '.mat;']);
    eval(['load ' freeVibDataDir char(phasematfiles(freqRatioInd)) '.mat;']);
    phase_ind = find(phasextoy > 180);
    phasextoy(phase_ind) = phasextoy(phase_ind)-360;
    eval(['load ' freeVibDataDir char(corrmatfiles(freqRatioInd)) '.mat;']);
    eval(['load ' freeVibDataDir char(ampfixmatfiles(freqRatioInd)) '.mat;']);

    vr_free = Vrn.*fy./(ypeakfreq./2./pi);

    vr_fr(freqRatioInd,:) = vr_free;
    Yad_fr(freqRatioInd,:) = Yad_fix;
    Xad_fr(freqRatioInd,:) = Xad_fix;
    Theta_fr(freqRatioInd,:) = phasextoy;
    vrn_fr(freqRatioInd,:) = Vrn;
% 
    %% Predict free vibrations

    % for i=1:length(FluidVel) % i is for different velocity
    for i=compIndx
       
        fluidSpeed = FluidVel(i);  

     % step 1 CF Initial Prediction    
        RigCylModelCFobj = RigCylModel(Diameter, Riserlength, massDensity_y(freqRatioInd), ...
            fluidSpeed, fluidDensity, springCons_y(freqRatioInd), damperCons_y(freqRatioInd));

        Vrn_y(freqRatioInd,i) = RigCylModelCFobj.FluidSpeed/(RigCylModelCFobj.NominalNaturalFreq * RigCylModelCFobj.Diameter);

        Vr0_y(freqRatioInd,i) = CalVr(RigCylModelCFobj, linearCaCFobj);
        DampingRatio0_y(freqRatioInd,i) = RigCylModelCFobj.DamperCons/(2*(RigCylModelCFobj.MassRatio+feval(linearCaCFobj.CaCF_A0_cfit,1/Vr0_y(freqRatioInd,i)))*...
            RigCylModelCFobj.FluidDensity*pi*RigCylModelCFobj.Diameter^2*RigCylModelCFobj.Length/4*2*pi*RigCylModelCFobj.NominalNaturalFreq);    
        NonDimAmp_linear0_y(freqRatioInd,i) = CalAmp_linear(RigCylModelCFobj, Vr0_y(freqRatioInd,i), Vrn_y(freqRatioInd,i), linearCeCFobj, linearCaCFobj);

     % step 2 IL Initial Prediction 

        RigCylModelILobj = RigCylModel(Diameter, Riserlength, massDensity_x(freqRatioInd), ...
            fluidSpeed, fluidDensity, springCons_x(freqRatioInd), damperCons_x(freqRatioInd));

        Vrn_x(freqRatioInd,i) = RigCylModelILobj.FluidSpeed/(RigCylModelILobj.NominalNaturalFreq * RigCylModelILobj.Diameter);    

        Vr0_x(freqRatioInd,i) = CalVr(RigCylModelILobj, linearCaILobj);
        DampingRatio0_x(freqRatioInd,i) = RigCylModelILobj.DamperCons/(2*(RigCylModelILobj.MassRatio+feval(CaILobj.surface_sfit,1/Vr0_x(freqRatioInd,i), 0 ))*...
            RigCylModelILobj.FluidDensity*pi*RigCylModelILobj.Diameter^2*RigCylModelILobj.Length/4*2*pi*RigCylModelILobj.NominalNaturalFreq);    

        NonDimAmp0_x(freqRatioInd,i) = CalAmp(RigCylModelILobj, Vr0_x(freqRatioInd,i), Vrn_x(freqRatioInd,i), CeILobj, linearCaILobj);
        NonDimAmp_linear0_x(freqRatioInd,i) = CalAmp_linear(RigCylModelILobj, Vr0_x(freqRatioInd,i), Vrn_x(freqRatioInd,i), linearCeILobj, linearCaILobj);
    end


    for i=compIndx
        
        fluidSpeed = FluidVel(i);  
        RigCylModel2Dobj = RigCylModel2D(Diameter, Riserlength, massDensity_x(freqRatioInd), massDensity_y(freqRatioInd), ...
            fluidSpeed, fluidDensity, springCons_x(freqRatioInd), springCons_y(freqRatioInd), damperCons_x(freqRatioInd), damperCons_y(freqRatioInd));
    %     RigCylModel2Dobj.PrintRigCylModel2D()    
        vrn_pred(freqRatioInd,i) = RigCylModel2Dobj.FluidSpeed/(RigCylModel2Dobj.NominalNaturalFreq_y * RigCylModel2Dobj.Diameter);

    % % use free vibration as stating point to solve the whole system, 2 step method 
    %     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse3(RigCylModel2Dobj, ...
    %         ILCFHydroModelobj,  Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i), vr_fr(freqRatioInd,i));

    % % use IL, CF independent solutions as stating point to solve the whole system, 2 step method     
    %     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse3(RigCylModel2Dobj, ...
    %         ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i), NonDimAmp_linear0_x(freqRatioInd,i),Vr0_y(freqRatioInd,i));
    % 
    % % use IL, CF independent solutions as stating point to solve the whole system, looping method   
    %     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse4(RigCylModel2Dobj, ...
    %         ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i), NonDimAmp_linear0_x(freqRatioInd,i),Vr0_y(freqRatioInd,i));


    % use free vibration as initial condition to calcuate theta and vr
    %     [vr_pred(freqRatioInd,i), vr_x_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i), vrerrmin(freqRatioInd,i)] =...
    %         CalVrTheta(RigCylModel2Dobj, ILCFHydroModelobj, Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i),vr_fr(freqRatioInd,i));

    % use IL, CF as initial condition to calcuate theta and vr
    %     [vr_pred(freqRatioInd,i), vr_x_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i), vrerrmin(freqRatioInd,i)] =...
    %         CalVrTheta(RigCylModel2Dobj, ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i),  NonDimAmp_linear0_x(freqRatioInd,i),Vr0_y(freqRatioInd,i));
    %     

    % use fsolve matlab function to solve the whole system

        i
        for j =1:numRun
            j
            x0=zeros(1,11);
            x0(1) = 1.5.*rand(1);
            x0(2) = 1.*rand(1);

    %         x0(1) = NonDimAmp_linear0_y(freqRatioInd,i);
    %         x0(2) = NonDimAmp_linear0_x(freqRatioInd,i);
            x0(3) = 360.*(rand(1)-0.5);
            x0(4) = Vr0_y(freqRatioInd,i);
            x0(5) = Vr0_x(freqRatioInd,i);


    %         x0(1) = 1;
    %         x0(2) = 0.3;
    %         x0(3) = 60;
    %         x0(4) = Vr0_y(freqRatioInd,i);
    %         x0(5) = Vr0_x(freqRatioInd,i);
    %         
            % free vib
    %         x0(1) = 0.9345;
    %         x0(2) = 0.1326;
    %         x0(3) = 62.176;
    %         x0(4) = 5.1480;
    %         x0(5) = x0(4) /2;

    %         x0(1) = 0.4107;
    %         x0(2) = 0.1389;
    %         x0(3) = 59.7928;
    %         x0(4) = 5.3014;
    %         x0(5) = x0(4) /2;
    %         

        %     x0(4) = Vrn_y(freqRatioInd,i);
        %     x0(5) = Vrn_x(freqRatioInd,i);
            x0(6) = getDataPoint(ILCFHydroModelobj, 'Cmy', x0(1), x0(2), x0(3),x0(4));
            x0(7) = getDataPoint(ILCFHydroModelobj, 'Cmx', x0(1), x0(2), x0(3),x0(4));
            x0(8) = getDataPoint(ILCFHydroModelobj, 'CLv', x0(1), x0(2), x0(3),x0(4));
            x0(9) = getDataPoint(ILCFHydroModelobj, 'CDv', x0(1), x0(2), x0(3),x0(4));
            x0(10) = RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ x0(6) )*...
                        RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y);         
            x0(11) = RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ x0(7) )*...
                        RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x);  

            options = optimset( 'Display', 'iter');

            f = @(x) myILCFfun(x, RigCylModel2Dobj, ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i));

%             [x(freqRatioInd,i,:), fval] = fsolve(f, x0, options);
            %         f0 = feval(f, x0);
             lb =[0 0 -180 0 0 -RigCylModel2Dobj.MassRatio_y -RigCylModel2Dobj.MassRatio_x 0 0 0 0];
             ub =[1.5 0.75 180 20 20 5 5 5 5 0.1 0.1];
            [x(freqRatioInd,i,j,:) res(i,j)] = lsqnonlin(f,x0,lb,ub);
        end
    end  


%     % % use fsolve matlab function to solve the whole system, given Vr
%         x0=zeros(1,9);
%     %     x0(1) = NonDimAmp_linear0_y(freqRatioInd,i);
%     %     x0(2) = NonDimAmp_linear0_x(freqRatioInd,i);
%         x0(1) = 1;
%         x0(2) = 0.3;
%         x0(3) = 60;
%         x0(4) = getDataPoint(ILCFHydroModelobj, 'Cmy', x0(1), x0(2), x0(3),vr_fr(freqRatioInd,i));
%         x0(5) = getDataPoint(ILCFHydroModelobj, 'Cmx', x0(1), x0(2), x0(3),vr_fr(freqRatioInd,i));
%         x0(6) = getDataPoint(ILCFHydroModelobj, 'CLv', x0(1), x0(2), x0(3),vr_fr(freqRatioInd,i));
%         x0(7) = getDataPoint(ILCFHydroModelobj, 'CDv', x0(1), x0(2), x0(3),vr_fr(freqRatioInd,i));
%         x0(8) = RigCylModel2Dobj.DamperCons_y/2/((RigCylModel2Dobj.MassRatio_y+ x0(4) )*...
%                     RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_y);         
%         x0(9) = RigCylModel2Dobj.DamperCons_x/2/((RigCylModel2Dobj.MassRatio_x+ x0(5) )*...
%                     RigCylModel2Dobj.FluidDensity*pi*RigCylModel2Dobj.Diameter^2*RigCylModel2Dobj.Length/4*2*pi*RigCylModel2Dobj.NominalNaturalFreq_x);  
%         options = optimset( 'Display', 'iter');
%     
%         f = @(x) myILCFwVrfun(x, RigCylModel2Dobj, ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i), vr_fr(freqRatioInd,i));
%     
%         [x(freqRatioInd,i,:), fval] = fsolve(f, x0, options);
%            
%         Yad_pred(freqRatioInd,i)    = x(freqRatioInd,i,1);
%         Xad_pred(freqRatioInd,i)    = x(freqRatioInd,i,2);
%         Theta_pred(freqRatioInd,i)  = x(freqRatioInd,i,3);
%         vr_pred(freqRatioInd,i)     = vr_fr(freqRatioInd,i);
%         vr_x_pred(freqRatioInd,i)   = 1/2*vr_fr(freqRatioInd,i);

    % Yad_pred(freqRatioInd,compIndx)
    % Xad_pred(freqRatioInd,compIndx)
    % Theta_pred(freqRatioInd,compIndx)
    % vr_pred(freqRatioInd,compIndx)
    % vr_x_pred(freqRatioInd,compIndx)
    
   %% plot comparison

   %% plot uncoupled solution
%     compIndx = 3:37;
    figure; 
    subplot(2,2,1)
    plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
        Vrn_y(freqRatioInd,compIndx),Vr0_y(freqRatioInd,compIndx)./Vrn_y(freqRatioInd,compIndx),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Crossflow Predicted from Forced Vib')
    axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,2)
    plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
        Vrn_y(freqRatioInd,compIndx),2.*Vr0_x(freqRatioInd,compIndx)./Vrn_y(freqRatioInd,compIndx),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', '2*Inline Predicted from Forced Vib')
    axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,3)
    plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',Vrn_y(freqRatioInd,compIndx),NonDimAmp_linear0_y(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,4)
    plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',Vrn_y(freqRatioInd,compIndx),NonDimAmp_linear0_x(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 0.6]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('Uncoupled Prediction vs Free Virbation Experiment Data')
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'unCoupled.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'unCoupled.fig'] );
   
   
   %% minimum residuals
%     for i=compIndx  
%         casenums = [1:numRun];
%         casenums = find(res(i,casenums)<1e-8);
%         figure
%         plot(reshape(x(freqRatioInd,i,casenums,1),1,length(casenums)),'-o')
%         hold on
%         plot(reshape(x(freqRatioInd,i,casenums,2),1,length(casenums)),'-*r')
%         plot(reshape(x(freqRatioInd,i,casenums,3),1,length(casenums))./180,'-*g')
%         plot(reshape(x(freqRatioInd,i,casenums,4),1,length(casenums))./Vrn(i),'-*m')
%         legend('Ay','Ax','theta/180','Vr/Vrn')
%         saveas(gcf, [outputDir filesep 'ResultRandomICcase' num2str(i) '.jpg'] )
%         saveas(gcf, [outputDir filesep 'ResultRandomICcase' num2str(i) '.fig'] )
% 
%         figure
%         plot(res(i,casenums),'-o')
%         title(['Residual Case' num2str(i)])
%         saveas(gcf, [outputDir filesep 'ResidualRandomICcase' num2str(i) '.jpg'] )
%         saveas(gcf, [outputDir filesep 'ResidualRandomICcase' num2str(i) '.fig'] )    
%     end

%     [res_min resmin_ind] = min(res,[],2);
%     for i = compIndx           
%         Yad_pred_res_min(freqRatioInd,i)    = x(freqRatioInd,i,resmin_ind(i),1);
%         Xad_pred_res_min(freqRatioInd,i)    = x(freqRatioInd,i,resmin_ind(i),2);
%         Theta_pred_res_min(freqRatioInd,i)  = x(freqRatioInd,i,resmin_ind(i),3);
%         vr_pred_res_min(freqRatioInd,i)     = x(freqRatioInd,i,resmin_ind(i),4);
%         vr_x_pred_res_min(freqRatioInd,i)   = x(freqRatioInd,i,resmin_ind(i),5);
%     end
%         
%     figure; 
%     subplot(2,2,1)
%     plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
%         vrn_pred(freqRatioInd,compIndx),vr_pred_res_min(freqRatioInd,compIndx)./vrn_pred(freqRatioInd,compIndx),'rs','MarkerSize',9);
%     xlabel('V_{rn}','FontSize',16);
%     ylabel('V_{r}/V_{rn}','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib')
%     axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
%     grid
%     subplot(2,2,2)
%     plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Theta_pred_res_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
%     xlabel('V_{rn}','FontSize',16);
%     ylabel('\theta','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib','Predicted from Forced Vib','Location','Southwest')
%     axis([4 10 -180 180]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
%     grid
%     subplot(2,2,3)
%     plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Yad_pred_res_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
%     xlabel('V_{rn}','FontSize',16);
%     ylabel('A_{y}/D','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
%     axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
%     grid
%     subplot(2,2,4)
%     plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Xad_pred_res_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
%     xlabel('V_{rn}','FontSize',16);
%     ylabel('A_{x}/D','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
%     axis([4 10 0 0.6]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
%     grid
%     suptitle('Min Residual Prediction vs Free Virbation Experiment Data')
%     saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'minRes.jpg'] );
%     saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'minRes.fig'] );

    
    for i=compIndx  
        casenums = [1:numRun];
%         casenums = find(x(freqRatioInd,i,casenums,4)<=8.1);
        casenums = find(res(i,casenums)<1e-5);

        result_conv = reshape(x(freqRatioInd,i,casenums,:),length(casenums),11);
        % scale theta to same scaling as others
        result_conv(:,3) = result_conv(:,3)./360;
        result_uniq = uniquetol(result_conv,0.1,'rows');
        result_uniq(:,3) = result_uniq(:,3).*360;     
        
        Ampy_pred{i}= result_uniq(:,1);
        Ampx_pred{i} = result_uniq(:,2);
        theta2D_pred{i} = result_uniq(:,3);
        vr2D_pred{i} = result_uniq(:,4);
        
        vr2D_pred_x{i} = result_uniq(:,5);
        Cmy_pred{i} = result_uniq(:,6);
        Cmx_pred{i} = result_uniq(:,7);
        CLv_pred{i} = result_uniq(:,8);
        CDv_pred{i} = result_uniq(:,9);
        Damp_pred_y{i} = result_uniq(:,10);
        Damp_pred_x{i} = result_uniq(:,11);
%         f0 = feval(f, result_uniq(1,:));

        for j = 1:length(Cmy_pred{i})
           AvePow_pred{i}(j) = getDataPoint(ILCFHydroModelobj, 'pow', Ampy_pred{i}(j), Ampx_pred{i}(j), theta2D_pred{i}(j), vr2D_pred{i}(j));                  
        end
       
        [pow_min(i), pow_min_ind(i)] = min(abs(AvePow_pred{i}));
        
        Yad_pred_pow_min(freqRatioInd,i)    = result_uniq(pow_min_ind(i),1);
        Xad_pred_pow_min(freqRatioInd,i)    = result_uniq(pow_min_ind(i),2);
        Theta_pred_pow_min(freqRatioInd,i)  = result_uniq(pow_min_ind(i),3);
        vr_pred_pow_min(freqRatioInd,i)     = result_uniq(pow_min_ind(i),4);
        vr_x_pred_pow_min(freqRatioInd,i)   = result_uniq(pow_min_ind(i),5);
    end
    
%% average power
   
   for i = compIndx
       figure
       for j = 1:length(Cmy_pred{i})
            plot(j, AvePow_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerFaceColor',col(j,:),'MarkerSize',9);
            hold on
       end
       xlim([0 length(Cmx_pred{i})+1])
       title(['Average Power Coefficient Case ' num2str(i)])
       saveas(gcf,[outputDir filesep 'PredAvePow' num2str(freqRatio(freqRatioInd)) 'Case' num2str(i) '.jpg'] );
       saveas(gcf,[outputDir filesep 'PredAvePow' num2str(freqRatio(freqRatioInd)) 'Case' num2str(i) '.fig'] );

   end
   
    figure; 
    subplot(2,2,1)
    plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
        vrn_pred(freqRatioInd,compIndx),vr_pred_pow_min(freqRatioInd,compIndx)./vrn_pred(freqRatioInd,compIndx),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,2)
    plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Theta_pred_pow_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib','Predicted from Forced Vib','Location','Southwest')
    axis([4 10 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    subplot(2,2,3)
    plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Yad_pred_pow_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    subplot(2,2,4)
    plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Xad_pred_pow_min(freqRatioInd,compIndx),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
    axis([4 10 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('Min Pow Prediction vs Free Virbation Experiment Data')
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'minPow.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'minPow.fig'] );

            
%% All solutions No color coding
    figure; 
    subplot(2,2,1)
    for i=compIndx
        figh{i} = plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}./vrn_pred(freqRatioInd,i),'rs','MarkerSize',9);
        hold on
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
    legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,2)
    for i=compIndx
        figh{i} = plot(vrn_pred(freqRatioInd,i),theta2D_pred{i},'rs','MarkerSize',9);
        hold on
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
    legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    
    subplot(2,2,3)
    for i=compIndx
        figh{i} = plot(vrn_pred(freqRatioInd,i),Ampy_pred{i},'rs','MarkerSize',9);
        hold on
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,4)
    for i=compIndx
        figh{i} = plot(vrn_pred(freqRatioInd,i),Ampx_pred{i},'rs','MarkerSize',9);
        hold on
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('All Prediction vs Free Virbation Experiment Data')    
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'all.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'all.fig'] );
    
            
%% All solutions color coding
    figure; 
    subplot(2,2,1)
    for i=compIndx
        for j = 1:length(vr2D_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j)./vrn_pred(freqRatioInd,i),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{rn}','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,2)
    for i=compIndx
        for j = 1:length(theta2D_pred{i})
            plot(vrn_pred(freqRatioInd,i),theta2D_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    
    subplot(2,2,3)
    for i=compIndx
        for j = 1:length(Ampy_pred{i})
            plot(vrn_pred(freqRatioInd,i),Ampy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid
    
    subplot(2,2,4)
    for i=compIndx
        for j = 1:length(Ampx_pred{i})
            plot(vrn_pred(freqRatioInd,i),Ampx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    suptitle('All Prediction vs Free Virbation Experiment Data')    
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'allColor.jpg'] );
    saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'allColor.fig'] );
    
        
%% Added mass verification
    
    figure; 
    subplot(2,3,1)
    for i=compIndx
        for j = 1:length(CLv_pred{i})
%             plot(vrn_pred(freqRatioInd,i),CLv_pred{i},'Color', col(j,:),'Marker',shp(j),'MarkerSize',9);
            plot(vrn_pred(freqRatioInd,i),CLv_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('CLv Y','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Southwest')
%     axis([4 10 0 1.5]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    
    subplot(2,3,4)
    for i=compIndx
        for j = 1:length(CDv_pred{i})
            plot(vrn_pred(freqRatioInd,i),CDv_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('CLv X','FontSize',16);
    set(gca,'FontSize',14);
    xlim([4 10])
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
%     axis([4 10 -180 180]);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    
    subplot(2,3,3)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j)./Vrn_y(freqRatioInd,i)./sqrt((massRatio_y(freqRatioInd)+Cmy_pred{i}(j))...
                ./(massRatio_y(freqRatioInd)+1)),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{r}(addedmass) Y','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 0.90 1.1]);
    set(gca,'Xtick',[4 5 6 7 8 9 10],'Ytick',[0.95 1 1.05]);
    grid
    
    subplot(2,3,6)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred_x{i}(j)./Vrn_x(freqRatioInd,i)./sqrt((massRatio_x(freqRatioInd)+Cmx_pred{i}(j))...
                ./(massRatio_x(freqRatioInd)+1)),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr2D_pred_x{i}./vrn_pred(freqRatioInd,i).*sqrt((massRatio_x(freqRatioInd)+Cmy_pred{i})./(massRatio_x(freqRatioInd)+1)),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{r}(addedmass) X','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
    axis([4 10 0.90 1.1]);
    set(gca,'Xtick',[4 5 6 7 8 9 10],'Ytick',[0.95 1 1.05]);
    grid
    
     subplot(2,3,2)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-');
    xlabel('V_{rn}','FontSize',16);
    ylabel('Cmy','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
%     axis([4 10 0.90 1.1]);
    xlim([4 10])
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    
    subplot(2,3,5)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
%     figh{i+1}=plot(vrn_fr(freqRatioInd,compIndx),vr2D_pred_x{i}./vrn_pred(freqRatioInd,i).*sqrt((massRatio_x(freqRatioInd)+Cmy_pred{i})./(massRatio_x(freqRatioInd)+1)),'rs','MarkerSize',9);
    xlabel('V_{rn}','FontSize',16);
    ylabel('Cmx','FontSize',16);
    set(gca,'FontSize',14);
%     legend([figh{i+1}, figh{compIndx(1)}(1)],{'Observed Free Vib', 'Predicted from Forced Vib'},'Location','Northwest')
%     axis([4 10 0.90 1.1]);
    xlim([4 10])
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid   
    
    suptitle('Prediction Verification')    
    saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'all.jpg'] );
    saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'all.fig'] );
        
    
%% Individual Amplitude veification    
   for i= compIndx
        figure; 
        subplot(4,3,1) 
        for j = 1:length(CLv_pred{i})
            plot(j, CLv_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
        ylabel('CLv Y','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid
        subplot(4,3,7)
        for j = 1:length(CDv_pred{i})
            plot(j, CDv_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('CLv X','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid

        subplot(4,3,2)
        for j = 1:length(Cmy_pred{i})
            plot(j,Damp_pred_y{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end  
        ylabel('Damp y','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid

        subplot(4,3,8)
        for j = 1:length(Cmx_pred{i})
            plot(j,Damp_pred_x{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Damp x','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid 

        subplot(4,3,3)
        for j = 1:length(Cmy_pred{i})
            plot(j,Cmy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Cmy','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid

        subplot(4,3,9)
        for j = 1:length(Cmx_pred{i})
            plot(j,Cmx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Cmx','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(CLv_pred{i})+1])
        grid   
        
        subplot(4,3,4) 
        for j = 1:length(vr2D_pred{i})
            plot(j, vr2D_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
        ylabel('Vr Y','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(vr2D_pred{i})+1])
        grid
        subplot(4,3,10)
        for j = 1:length(vr2D_pred_x{i})
            plot(j, vr2D_pred_x{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Vr X','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(vr2D_pred_x{i})+1])
        grid
        subplot(4,3,5) 
        for j = 1:length(Ampy_pred{i})
            plot(j, Ampy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
        ylabel('Amp Y','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(Ampy_pred{i})+1])
        grid
        subplot(4,3,11)
        for j = 1:length(Ampx_pred{i})
            plot(j, Ampx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end 
        ylabel('Amp X','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(Ampx_pred{i})+1])
        grid        
        
        subplot(4,3,6)
        for j = 1:length(Cmy_pred{i})
            plot(j, Ampy_pred{i}(j)./(CLv_pred{i}(j).*vr2D_pred{i}(j).*Vrn_y(freqRatioInd,i)/...
                (4*pi^3*(massRatio_y(freqRatioInd)+Cmy_pred{i}(j)).*Damp_pred_y{i}(j))),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
        ylabel('Amp Y / Amp Y(lift, added mess)','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(Cmy_pred{i})+1])
        ylim([0.95 1.05])
        grid

        subplot(4,3,12)
        for j = 1:length(Cmx_pred{i})
            plot(j, Ampx_pred{i}(j)./(CDv_pred{i}(j).*vr2D_pred_x{i}(j).*Vrn_x(freqRatioInd,i)/...
                (4*pi^3*(massRatio_x(freqRatioInd)+Cmx_pred{i}(j)).*Damp_pred_x{i}(j))),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end   
        ylabel('Amp X / Amp X(lift, added mess)','FontSize',16);
        set(gca,'FontSize',14);
        xlim([0 length(Cmx_pred{i})+1])
        ylim([0.95 1.05])
        grid
        
        
        suptitle(['Prediction Verification Case ' num2str(i)])    
        saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'Case' num2str(i) '.jpg'] );
        saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'Case' num2str(i) '.fig'] );
   end  
   


   %% Added mass verification
    figure
    subplot(2,3,1)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmy_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('Cmy','FontSize',16);
    set(gca,'FontSize',14);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid    
    subplot(2,3,4)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),Cmx_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('Cmx','FontSize',16);
    set(gca,'FontSize',14);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid   
    subplot(2,3,2)
    for i=compIndx
        for j = 1:length(vr2D_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('Vry','FontSize',16);
    set(gca,'FontSize',14);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid    
    subplot(2,3,5)
    for i=compIndx
        for j = 1:length(vr2D_pred_x{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred_x{i}(j),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('Vrx','FontSize',16);
    set(gca,'FontSize',14);
    set(gca,'Xtick',[4 5 6 7 8 9 10]);
    grid
    subplot(2,3,3)
    for i=compIndx
        for j = 1:length(Cmy_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred{i}(j)./Vrn_y(freqRatioInd,i)./sqrt((massRatio_y(freqRatioInd)+Cmy_pred{i}(j))...
                ./(massRatio_y(freqRatioInd)+1)),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{r}(addedmass) Y','FontSize',16);
    set(gca,'FontSize',14);
    axis([4 10 0.90 1.1]);
    set(gca,'Xtick',[4 5 6 7 8 9 10],'Ytick',[0.95 1 1.05]);
    grid
    
    subplot(2,3,6)
    for i=compIndx
        for j = 1:length(Cmx_pred{i})
            plot(vrn_pred(freqRatioInd,i),vr2D_pred_x{i}(j)./Vrn_x(freqRatioInd,i)./sqrt((massRatio_x(freqRatioInd)+Cmx_pred{i}(j))...
                ./(massRatio_x(freqRatioInd)+1)),'Color', col(j,:),'Marker',shp(2),'MarkerSize',9);
            hold on
        end
    end    
    xlabel('V_{rn}','FontSize',16);
    ylabel('V_{r}/V_{r}(addedmass) X','FontSize',16);
    set(gca,'FontSize',14);
    axis([4 10 0.90 1.1]);
    set(gca,'Xtick',[4 5 6 7 8 9 10],'Ytick',[0.95 1 1.05]);
    grid
    suptitle('Prediction Verification Cm Vr')    
    saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'CmVr.jpg'] );
    saveas(gcf,[outputDir filesep 'PredVeri' num2str(freqRatio(freqRatioInd)) 'CmVr.fig'] );
    
%% minimum lift coefficients


    

%     for i=compIndx  
%         Yad_pred_ave(freqRatioInd,i)    = mean(Ampy_pred{i});
%         Xad_pred_ave(freqRatioInd,i)    = mean(Ampx_pred{i});
%         Theta_pred_ave(freqRatioInd,i)  = mean(theta2D_pred{i});
%         vr_pred_ave(freqRatioInd,i)     = mean(vr2D_pred{i});
%     end
%     figure; 
%     subplot(2,2,1)
%     plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',...
%         vrn_pred(freqRatioInd,compIndx),vr_pred_ave(freqRatioInd,compIndx)./vrn_pred(freqRatioInd,compIndx),'rs','MarkerSize',9);
%     xlabel('V_{rn}','FontSize',16);
%     ylabel('V_{r}/V_{rn}','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib')
%     axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
%     grid
%     subplot(2,2,2)
%     plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Theta_pred_ave(freqRatioInd,compIndx),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
%     xlabel('V_{rn}','FontSize',16);
%     ylabel('\theta','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib','Predicted from Forced Vib','Location','Southwest')
%     axis([4 10 -180 180]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
%     grid
%     subplot(2,2,3)
%     plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Yad_pred_ave(freqRatioInd,compIndx),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
%     xlabel('V_{rn}','FontSize',16);
%     ylabel('A_{y}/D','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
%     axis([4 10 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
%     grid
%     subplot(2,2,4)
%     plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Xad_pred_ave(freqRatioInd,compIndx),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
%     xlabel('V_{rn}','FontSize',16);
%     ylabel('A_{x}/D','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
%     axis([4 10 0 0.6]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
%     grid
%     suptitle('Mean Prediction vs Free Virbation Experiment Data')
%     saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'mean.jpg'] );
%     saveas(gcf,[outputDir filesep 'PredvsFree' num2str(freqRatio(freqRatioInd)) 'mean.fig'] );
% 
%         
end

% figure; clf;
% plot(vrn_fr(freqRatioInd,compIndx),vr_fr(freqRatioInd,compIndx)./vrn_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),vr_pred(freqRatioInd,compIndx)./vrn_pred(freqRatioInd,compIndx),'rs','MarkerSize',9);
% %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
% xlabel('V_{rn}','FontSize',16);
% ylabel('V_{r}/V_{rn}','FontSize',16);
% set(gca,'FontSize',14);
% legend('Observed Free Vib', 'Predicted from Forced Vib')
% axis([4 10 0 1.5]);
% set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
% grid
% saveas(gcf,[outputDir filesep 'Vr3' '.jpg'] );
% saveas(gcf,[outputDir filesep 'Vr3' '.fig'] );
% % 
% % figure; clf;
% % plot(vrn_pred(freqRatioInd,compIndx),vrerrmin(freqRatioInd,compIndx),'bo-', vrn_pred(freqRatioInd,compIndx),vr_pred(freqRatioInd,compIndx)./(vr_x_pred(freqRatioInd,compIndx)),'rs','MarkerSize',9);
% % %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
% % xlabel('V_{rn}','FontSize',16);
% % ylabel('V_{ry}/V_{rx}','FontSize',16);
% % set(gca,'FontSize',14);
% % legend( 'Err', 'Vry/Vrx Ratio')
% % axis([4 10 0 2.5]);
% % set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.5 1 1.5 2 2.5]);
% % grid
% % saveas(gcf,[outputDir filesep 'VrXY3' '.jpg'] );
% % saveas(gcf,[outputDir filesep 'VrXY3' '.fig'] );
% 
% figure; clf;
% plot(vrn_fr(freqRatioInd,compIndx),Theta_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Theta_pred(freqRatioInd,compIndx),'rs','MarkerSize',9);
% %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
% xlabel('V_{rn}','FontSize',16);
% ylabel('\theta','FontSize',16);
% set(gca,'FontSize',14);
% legend('Observed Free Vib','Predicted from Forced Vib','Location','Northwest')
% axis([4 10 -180 180]);
% set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
% grid
% saveas(gcf,[outputDir filesep 'Theta3' '.jpg'] );
% saveas(gcf,[outputDir filesep 'Theta3' '.fig'] );
% 
% 
% figure; clf;
% plot(vrn_fr(freqRatioInd,compIndx),Yad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Yad_pred(freqRatioInd,compIndx),'rs','MarkerSize',9);
% %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
% xlabel('V_{rn}','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% set(gca,'FontSize',14);
% legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
% axis([4 10 0 1.5]);
% set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
% grid
% saveas(gcf,[outputDir filesep 'Yad3' '.jpg'] );
% saveas(gcf,[outputDir filesep 'Yad3' '.fig'] );
% 
% figure; clf;
% plot(vrn_fr(freqRatioInd,compIndx),Xad_fr(freqRatioInd,compIndx),'bo-',vrn_pred(freqRatioInd,compIndx),Xad_pred(freqRatioInd,compIndx),'rs','MarkerSize',9);
% %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(freqRatioInd)),'FontSize',16);
% xlabel('V_{rn}','FontSize',16);
% ylabel('A_{x}/D','FontSize',16);
% set(gca,'FontSize',14);
% legend('Observed Free Vib', 'Predicted from Forced Vib','Location','Northwest')
% axis([4 10 0 0.6]);
% set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9 9.5 10],'Ytick',[0 0.15 0.3 0.45 0.6]);
% grid
% saveas(gcf,[outputDir filesep 'Xad3' '.jpg'] );
% saveas(gcf,[outputDir filesep 'Xad3' '.fig'] );




%  % step 3 IL-CF Initial Prediction 
%     
%     RigCylModel2Dobj = RigCylModel2D(Diameter, Riserlength, massDensity_x(freqRatioInd), massDensity_y(freqRatioInd), ...
%         fluidSpeed, fluidDensity, springCons_x(freqRatioInd), springCons_y(freqRatioInd), damperCons_x(freqRatioInd), damperCons_y(freqRatioInd));
%     RigCylModel2Dobj.PrintRigCylModel2D()
%     
%     vrn_pred(freqRatioInd,i) = RigCylModel2Dobj.FluidSpeed/(RigCylModel2Dobj.NominalNaturalFreq_y * RigCylModel2Dobj.Diameter);
% %     Yad_pred(freqRatioInd,i) =0.5;
% %     Xad_pred(freqRatioInd,i) =0.1;
% %     vr_pred(freqRatioInd,i)  =vrn_pred(freqRatioInd,i);
% %     Theta_pred(freqRatioInd,i) = 0;
% 
%     
%     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse(RigCylModel2Dobj, ...
%         ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i), NonDimAmp_linear0_x(freqRatioInd,i), Vr0_y(freqRatioInd,i));
% 
%     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse3(RigCylModel2Dobj, ...
%         ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i), NonDimAmp_linear0_x(freqRatioInd,i), Vr0_y(freqRatioInd,i));
% 
%     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse2(RigCylModel2Dobj, ...
%         ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i), NonDimAmp_linear0_x(freqRatioInd,i), Vr0_y(freqRatioInd,i));
% 
% %     [y_pred, x_pred, Vr_pred, theta_pred] = CalResponse(RigCylModel2Dobj, ILCFHydroModelobj, 0.5, 0.1, 6)    
% 
% 
% 


