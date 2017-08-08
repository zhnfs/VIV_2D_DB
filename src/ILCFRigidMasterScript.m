% function ILCFRigidMasterScript()

clear all
close all
% createCeList();

col = distinguishable_colors(100);
% col = col([13,15,3,8,1,10,14,7,4],:);
shp = ['o';'s';'p';'v';'*';'d';'^';'>';'<';'.';'o' ];
shp = [shp ;shp ;shp ;shp ;shp ;shp ;shp ;shp ;shp ;shp];

utility = Utility.getInstance();
filesep = utility.getFilesep();
[vivaRoot, optRoot, experimentRoot, fatigueRoot, utilityRoot, InlineVIVRoot, CrossflowVIVRoot, mooringRoot, ILCFVIVRoot] = getRoots();

outputDir =[ILCFVIVRoot filesep 'output' filesep 'ICFLMotion_JMD_RobustSmooth_New_run10'];
utility.conditionalMkdir(outputDir);

numRun = 10;

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

ILCFHydroModelobj = ILCFHydroModelJMD();
% ILCFHydroModelobj = ILCFHydroModelHZ();
% ILCFHydroModelobj = ILCFHydroModelZDMiT();


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
% freqRatio = [1	1.22	1.37	1.52	1.67	2];
freqRatioNames ={ '1p0', '1p22', '1p37', '1p52', '1p67', '1p9'};
fn_y = [0.715	0.799	0.894	0.977	0.698	0.704];
fn_x = fn_y .*freqRatio;

% fn_x = [0.715	0.974   1.226   1.488   1.167   1.336];

springCons_y = [697  902 1118 1386 967 1013];
% springCons_x = [609  1301 2130 2850 2580 3235];

springCons_x = springCons_y .* freqRatio.^2;

displacedMass = fluidDensity * pi*(Diameter/2)^2*Riserlength;

massRatio_y_fr = (springCons_y./(fn_y*2*pi).^2)/displacedMass-1;
% massRatio_x = (springCons_x./(fn_x*2*pi).^2)/displacedMass-1;
massRatio_x_fr = massRatio_y_fr;

% % begin adjust mass ratio
mRatioMulti = [1 1 1 1 1 1];
% mRatioMulti = 1./massRatio_y_fr;
% mRatio = massRatio_y(1);
massRatio_y = mRatioMulti .* massRatio_y_fr;
massRatio_x = mRatioMulti .* massRatio_x_fr;
% % end adjust mass ratio
    
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

% targetFreqRatio = [1.5:0.1:2];
targetFreqRatio = 2*ones(1,6);


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

%% Predict free vibrations

% compIndx = 1:37;
% compIndx = 15;
compIndx = [9:25];

for freqRatioInd= 6:-1:2
    % for i=1:length(FluidVel) % i is for different velocity
    
    % load Free Vibration Data
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

    eval(['FluidVel(freqRatioInd,:) = vels' freqRatioNames{freqRatioInd} ';'])

    % begin adjust Fluid velocity to adjust Vrn for mass ratio change
    FluidVel(freqRatioInd,:) = FluidVel(freqRatioInd,:) .* sqrt((massRatio_y_fr(freqRatioInd)+1)/(massRatio_y(freqRatioInd)+1));
    % end added for mass ratio 1

%     % begin extend Vrn to 20    
%     FluidVel = [FluidVel FluidVel(end).*[13:20]./12];
%     % end extend Vrn to 20  

    for i=compIndx
    
     % step 1 CF Initial Prediction    
        RigCylModelCFobj(freqRatioInd,i) = RigCylModel(Diameter, Riserlength, massDensity_y(freqRatioInd), ...
            FluidVel(freqRatioInd,i), fluidDensity, springCons_y(freqRatioInd), damperCons_y(freqRatioInd));

        Vrn_y(freqRatioInd,i) = RigCylModelCFobj(freqRatioInd,i).FluidSpeed/(RigCylModelCFobj(freqRatioInd,i).NominalNaturalFreq * RigCylModelCFobj(freqRatioInd,i).Diameter);

        Vr0_y(freqRatioInd,i) = CalVr(RigCylModelCFobj(freqRatioInd,i), linearCaCFobj);
        DampingRatio0_y(freqRatioInd,i) = RigCylModelCFobj(freqRatioInd,i).DamperCons/(2*(RigCylModelCFobj(freqRatioInd,i).MassRatio+feval(linearCaCFobj.CaCF_A0_cfit,1/Vr0_y(freqRatioInd,i)))*...
            RigCylModelCFobj(freqRatioInd,i).FluidDensity*pi*RigCylModelCFobj(freqRatioInd,i).Diameter^2*RigCylModelCFobj(freqRatioInd,i).Length/4*2*pi*RigCylModelCFobj(freqRatioInd,i).NominalNaturalFreq);    
        NonDimAmp_linear0_y(freqRatioInd,i) = CalAmp_linear(RigCylModelCFobj(freqRatioInd,i), Vr0_y(freqRatioInd,i), Vrn_y(freqRatioInd,i), linearCeCFobj, linearCaCFobj);

     % step 2 IL Initial Prediction 

        RigCylModelILobj(freqRatioInd,i) = RigCylModel(Diameter, Riserlength, massDensity_x(freqRatioInd), ...
            FluidVel(freqRatioInd,i), fluidDensity, springCons_x(freqRatioInd), damperCons_x(freqRatioInd));

        Vrn_x(freqRatioInd,i) = RigCylModelILobj(freqRatioInd,i).FluidSpeed/(RigCylModelILobj(freqRatioInd,i).NominalNaturalFreq * RigCylModelILobj(freqRatioInd,i).Diameter);    

        Vr0_x(freqRatioInd,i) = CalVr(RigCylModelILobj(freqRatioInd,i), linearCaILobj);
        DampingRatio0_x(freqRatioInd,i) = RigCylModelILobj(freqRatioInd,i).DamperCons/(2*(RigCylModelILobj(freqRatioInd,i).MassRatio+feval(CaILobj.surface_sfit,1/Vr0_x(freqRatioInd,i), 0 ))*...
            RigCylModelILobj(freqRatioInd,i).FluidDensity*pi*RigCylModelILobj(freqRatioInd,i).Diameter^2*RigCylModelILobj(freqRatioInd,i).Length/4*2*pi*RigCylModelILobj(freqRatioInd,i).NominalNaturalFreq);    

        NonDimAmp0_x(freqRatioInd,i) = CalAmp(RigCylModelILobj(freqRatioInd,i), Vr0_x(freqRatioInd,i), Vrn_x(freqRatioInd,i), CeILobj, linearCaILobj);
        NonDimAmp_linear0_x(freqRatioInd,i) = CalAmp_linear(RigCylModelILobj(freqRatioInd,i), Vr0_x(freqRatioInd,i), Vrn_x(freqRatioInd,i), linearCeILobj, linearCaILobj);

         
        RigCylModel2Dobj(freqRatioInd,i) = RigCylModel2D(Diameter, Riserlength, massDensity_x(freqRatioInd), massDensity_y(freqRatioInd), ...
            FluidVel(freqRatioInd,i), fluidDensity, springCons_x(freqRatioInd), springCons_y(freqRatioInd), damperCons_x(freqRatioInd), damperCons_y(freqRatioInd));
    %     RigCylModel2Dobj(freqRatioInd,i).PrintRigCylModel2D()    
        vrn_pred(freqRatioInd,i) = RigCylModel2Dobj(freqRatioInd,i).FluidSpeed/(RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_y * RigCylModel2Dobj(freqRatioInd,i).Diameter);
    
    % % use free vibration as stating point to solve the whole system, 2 step method 
    %     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse3(RigCylModel2Dobj(freqRatioInd,i), ...
    %         ILCFHydroModelobj,  Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i), vr_fr(freqRatioInd,i));

    % % use IL, CF independent solutions as stating point to solve the whole system, 2 step method     
    %     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse3(RigCylModel2Dobj(freqRatioInd,i), ...
    %         ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i), NonDimAmp_linear0_x(freqRatioInd,i),Vr0_y(freqRatioInd,i));
    % 
    % % use IL, CF independent solutions as stating point to solve the whole system, looping method   
    %     [Yad_pred(freqRatioInd,i), Xad_pred(freqRatioInd,i), vr_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i)] = CalResponse4(RigCylModel2Dobj(freqRatioInd,i), ...
    %         ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i), NonDimAmp_linear0_x(freqRatioInd,i),Vr0_y(freqRatioInd,i));


    % use free vibration as initial condition to calcuate theta and vr
    %     [vr_pred(freqRatioInd,i), vr_x_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i), vrerrmin(freqRatioInd,i)] =...
    %         CalVrTheta(RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i),vr_fr(freqRatioInd,i));

    % use IL, CF as initial condition to calcuate theta and vr
    %     [vr_pred(freqRatioInd,i), vr_x_pred(freqRatioInd,i), Theta_pred(freqRatioInd,i), vrerrmin(freqRatioInd,i)] =...
    %         CalVrTheta(RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, NonDimAmp_linear0_y(freqRatioInd,i),  NonDimAmp_linear0_x(freqRatioInd,i),Vr0_y(freqRatioInd,i));
    %   
    
 
    % use fsolve matlab function to solve the whole system

        i
        for j =1:numRun
%             j
            x0=zeros(1,11);
            x0(1) = 1.5.*rand(1);
            x0(2) =1.*rand(1);
    %         x0(1) = NonDimAmp_linear0_y(freqRatioInd,i);
    %         x0(2) = NonDimAmp_linear0_x(freqRatioInd,i);
            x0(3) = 360./180.*(rand(1)-0.5);
%             x0(3) = 360.*(rand(1)-0.5);
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
            x0(10) = RigCylModel2Dobj(freqRatioInd,i).DamperCons_y/2/((RigCylModel2Dobj(freqRatioInd,i).MassRatio_y+ x0(6) )*...
                        RigCylModel2Dobj(freqRatioInd,i).FluidDensity*pi*RigCylModel2Dobj(freqRatioInd,i).Diameter^2*RigCylModel2Dobj(freqRatioInd,i).Length/4*2*pi*RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_y);         
            x0(11) = RigCylModel2Dobj(freqRatioInd,i).DamperCons_x/2/((RigCylModel2Dobj(freqRatioInd,i).MassRatio_x+ x0(7) )*...
                        RigCylModel2Dobj(freqRatioInd,i).FluidDensity*pi*RigCylModel2Dobj(freqRatioInd,i).Diameter^2*RigCylModel2Dobj(freqRatioInd,i).Length/4*2*pi*RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_x);  

            options = optimset( 'Display', 'iter');

            f = @(x) myILCFfun(x, RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i), targetFreqRatio(freqRatioInd));

%             [x(freqRatioInd,i,:), fval] = fsolve(f, x0, options);
            %         f0 = feval(f, x0);
            
%              lb =[0.25 0 -1 4.5 2 -RigCylModel2Dobj(freqRatioInd,i).MassRatio_y -RigCylModel2Dobj(freqRatioInd,i).MassRatio_x 0 0 0 0];
%              ub =[1.5 0.75 1 20 20 5 5 5 5 0.1 0.1];
             lb =[0 0 -1 0 0 -RigCylModel2Dobj(freqRatioInd,i).MassRatio_y -RigCylModel2Dobj(freqRatioInd,i).MassRatio_x 0 0 0 0];
             ub =[1.5 0.75 1 10 10 5 5 5 5 0.1 0.1];
             
            [x(freqRatioInd,i,j,:) res(freqRatioInd, i,j)] = lsqnonlin(f,x0,lb,ub);
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
%         x0(8) = RigCylModel2Dobj(freqRatioInd,i).DamperCons_y/2/((RigCylModel2Dobj(freqRatioInd,i).MassRatio_y+ x0(4) )*...
%                     RigCylModel2Dobj(freqRatioInd,i).FluidDensity*pi*RigCylModel2Dobj(freqRatioInd,i).Diameter^2*RigCylModel2Dobj(freqRatioInd,i).Length/4*2*pi*RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_y);         
%         x0(9) = RigCylModel2Dobj(freqRatioInd,i).DamperCons_x/2/((RigCylModel2Dobj(freqRatioInd,i).MassRatio_x+ x0(5) )*...
%                     RigCylModel2Dobj(freqRatioInd,i).FluidDensity*pi*RigCylModel2Dobj(freqRatioInd,i).Diameter^2*RigCylModel2Dobj(freqRatioInd,i).Length/4*2*pi*RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_x);  
%         options = optimset( 'Display', 'iter');
%     
%         f = @(x) myILCFwVrfun(x, RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i), vr_fr(freqRatioInd,i));
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
        
         for i=compIndx
            casenums = [1:numRun];
    %         casenums = find(x(freqRatioInd,i,casenums,4)<=8.1);
    %         casenums = find(res(freqRatioInd, i,casenums)<1e-3);

            result_conv = reshape(x(freqRatioInd,i,casenums,:),length(casenums),11);
            % scale theta to same scaling as others
    %         result_conv(:,3) = result_conv(:,3)./360;
            result_uniq = uniquetol(result_conv,0.1,'rows');
    %         result_uniq(:,3) = result_uniq(:,3).*360;     

           % remove solutions due to numerical erros
           j=1;
           f2 = @(x) myILCFfun2(x, RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i), targetFreqRatio(freqRatioInd));
           while j <=size(result_uniq,1)
                f0 = feval(f2, result_uniq(j,:));

                if max(abs(f0(1:5)))>0.1
                    result_uniq(j,:) = [];
                else
                    j = j+1;
                end
           end

            f = @(x) myILCFfun(x, RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i), targetFreqRatio(freqRatioInd));
           j=1;
           while j <=size(result_uniq,1)
                f0 = feval(f, result_uniq(j,:));

                if max(abs(f0(1:5)))>0.1
                    result_uniq(j,:) = [];
                else
                    j = j+1;
                end
           end
    % %        
           % remove solutions has Amp Y <0.15  for HZ  0.25 for JMD,
    %        result_uniq(result_uniq(:,1)<0.15,:) =[];
    %        result_uniq(result_uniq(:,2)<0.1,:) =[];

    %        % incase no solution satisified all requirement, use the one with min residuals     
    %        if isempty(result_uniq)
    %            casenums = [1:numRun];
    %            [temp, casenum_minpow] = min(res(freqRatioInd, i,casenums));
    %            result_uniq =  reshape(x(freqRatioInd,i,casenum_minpow,:),1,11);
    %        end
    %     
           % sort solutions by the Amp Y

           result_uniq = sortrows(result_uniq,-1);

    %        %% remove solutions has Vr out of scope
    %        result_uniq(result_uniq(:,4)>8,:) =[];

           % save solutions
            if ~isempty(result_uniq)
                Ampy_pred{i}= result_uniq(:,1);
                Ampx_pred{i} = result_uniq(:,2);
                theta2D_pred{i} = result_uniq(:,3).*180;
                vr2D_pred{i} = result_uniq(:,4);

                vr2D_pred_x{i} = result_uniq(:,5);
                Cmy_pred{i} = result_uniq(:,6);
                Cmx_pred{i} = result_uniq(:,7);
                CLv_pred{i} = result_uniq(:,8);
                CDv_pred{i} = result_uniq(:,9);
                Damp_pred_y{i} = result_uniq(:,10);
                Damp_pred_x{i} = result_uniq(:,11);
            else
                Ampy_pred{i}= NaN;
                Ampx_pred{i} = NaN;
                theta2D_pred{i} = NaN;
                vr2D_pred{i} = NaN;

                vr2D_pred_x{i} = NaN;
                Cmy_pred{i} = NaN;
                Cmx_pred{i} = NaN;
                CLv_pred{i} = NaN;
                CDv_pred{i} = NaN;
                Damp_pred_y{i} = NaN;
                Damp_pred_x{i} = NaN;
            end



            % Min Ave Power
            AvePow_pred{i} = [];
            for j = 1:length(Ampy_pred{i})
               AvePow_pred{i}(j) = getDataPoint(ILCFHydroModelobj, 'pow', Ampy_pred{i}(j), Ampx_pred{i}(j), theta2D_pred{i}(j)./180, vr2D_pred{i}(j));                  
            end

            [pow_min(i), pow_min_ind(i)] = min(abs(AvePow_pred{i}));

            Yad_pred_pow_min(freqRatioInd,i)    = Ampy_pred{i}(pow_min_ind(i));
            Xad_pred_pow_min(freqRatioInd,i)    = Ampx_pred{i}(pow_min_ind(i));
            Theta_pred_pow_min(freqRatioInd,i)  = theta2D_pred{i}(pow_min_ind(i));
            vr_pred_pow_min(freqRatioInd,i)     = vr2D_pred{i}(pow_min_ind(i));
            vr_x_pred_pow_min(freqRatioInd,i)   = vr2D_pred_x{i}(pow_min_ind(i));

            % Max Response Amplitude in CF

            [MaxAmpy(i), MaxAmpy_ind(i)] = max(Ampy_pred{i});

            Yad_pred_MaxAmpy(freqRatioInd,i)    = Ampy_pred{i}(MaxAmpy_ind(i));
            Xad_pred_MaxAmpy(freqRatioInd,i)    = Ampx_pred{i}(MaxAmpy_ind(i));
            Theta_pred_MaxAmpy(freqRatioInd,i)  = theta2D_pred{i}(MaxAmpy_ind(i));
            vr_pred_MaxAmpy(freqRatioInd,i)    = vr2D_pred{i}(MaxAmpy_ind(i));
            vr_x_pred_MaxAmpy(freqRatioInd,i)   = vr2D_pred_x{i}(MaxAmpy_ind(i));

         end
    
        % Min Residual
        [res_min resmin_ind] = min(res,[],3);
        for i = compIndx           
            Yad_pred_res_min(freqRatioInd,i)    = x(freqRatioInd,i,resmin_ind(i),1);
            Xad_pred_res_min(freqRatioInd,i)    = x(freqRatioInd,i,resmin_ind(i),2);
            Theta_pred_res_min(freqRatioInd,i)  = x(freqRatioInd,i,resmin_ind(i),3);
            vr_pred_res_min(freqRatioInd,i)     = x(freqRatioInd,i,resmin_ind(i),4);
            vr_x_pred_res_min(freqRatioInd,i)   = x(freqRatioInd,i,resmin_ind(i),5);
        end

       %% plot comparison
    %     compIndx = [9:25];
    
       % plot uncoupled solution
        plot_uncouple(vrn_fr, vr_fr, Yad_fr, Xad_fr, Vrn_y, Vr0_y,Vr0_x, NonDimAmp_linear0_y,NonDimAmp_linear0_x, freqRatio, freqRatioInd, compIndx, outputDir)
       close
       
       % minimum residuals
        plot_min_res(vrn_fr, vr_fr, Yad_fr, Xad_fr, Theta_fr,vrn_pred, vr_pred_res_min, Yad_pred_res_min, Xad_pred_res_min, Theta_pred_res_min, x, res, numRun, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

       % plot vr
        plot_vr(vrn_fr, vr_fr, vrn_pred, vr2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

       % min average power
        plot_min_pow(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred, AvePow_pred, vr_pred_pow_min, Yad_pred_pow_min, Xad_pred_pow_min, Theta_pred_pow_min, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)
        close

       % plot max amp y solution   
        plot_maxY(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred, vr_pred_MaxAmpy, Yad_pred_MaxAmpy, Xad_pred_MaxAmpy, Theta_pred_MaxAmpy, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

       % All solutions vs vrn
        plot_All_vrn(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred,vr2D_pred, Ampy_pred, Ampx_pred, theta2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

       % All solutions vs vr
        plot_All_vr(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred,vr2D_pred, Ampy_pred, Ampx_pred, theta2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)
        close all

       % Added mass verification
        plot_addedMass(CLv_pred, CDv_pred, Cmy_pred, Cmx_pred, vr2D_pred_x, Vrn_y, Vrn_x, ...
        massRatio_y, massRatio_x, vrn_pred,vr2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

       % Individual Amplitude veification               
        plot_Amp_ver(CLv_pred, CDv_pred, Cmy_pred, Cmx_pred, vr2D_pred_x, Vrn_y, Vrn_x, massRatio_y, massRatio_x,...
        Damp_pred_y, Ampy_pred, Ampx_pred, Damp_pred_x, vr2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)
        close all

    %    % plot combined max amp y solution
    %     plot_combine(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, Vrn_y, Vr0_y,Vr0_x, NonDimAmp_linear0_y,NonDimAmp_linear0_x, ...
    %         vrn_pred, vr_pred_MaxAmpy, Yad_pred_MaxAmpy, Xad_pred_MaxAmpy, Theta_pred_MaxAmpy, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)
    %  
end
     
save ICFLMotionJMD_RobustSmooth_run10.mat


