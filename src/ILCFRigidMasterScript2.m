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

ILCFVIVRoot = '/Users/haining/Dropbox/viv/src/ILCFprediction';

outputDir =[ILCFVIVRoot filesep 'output' filesep 'ICFLMotion_ZDMiT_20130908_run10'];
utility.conditionalMkdir(outputDir);

compIndx = [9:25];

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

% ILCFHydroModelobj = ILCFHydroModelJMD();
% ILCFHydroModelobj = ILCFHydroModelHZ();
ILCFHydroModelobj = ILCFHydroModelZDMiT();


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
springCons_x = [609  1301 2130 2850 2580 3235];

% springCons_x = springCons_y .* freqRatio.^2;

displacedMass = fluidDensity * pi*(Diameter/2)^2*Riserlength;

massRatio_y_fr = (springCons_y./(fn_y*2*pi).^2)/displacedMass-1;
massRatio_x_fr = (springCons_x./(fn_x*2*pi).^2)/displacedMass-1;
% massRatio_x_fr = massRatio_y_fr;

% % begin adjust mass ratio
mRatioMulti = [1 1 1 1 1 1]; % original settings
% mRatioMulti = 1./massRatio_y_fr;
% mRatio = massRatio_y(1);
massRatio_y = mRatioMulti .* massRatio_y_fr;
massRatio_x = mRatioMulti .* massRatio_x_fr;
% % end adjust mass ratio
    
% massRatio_y = (springCons_y./(fn_y*2*pi).^2)/displacedMass;
% massRatio_x = (springCons_x./(fn_x*2*pi).^2)/displacedMass;
% 
% strucDampRatio_y =0.00001*[2.2	1.3	1.1	1.6	2.6	6.2];
% strucDampRatio_x =0.00001*[2.2	1.7	2.5	3.2	2.9	2.5];

% Original Damping Ratio
strucDampRatio_y =0.01*[2.2	1.3	1.1	1.6	2.6	6.2];
strucDampRatio_x =0.01*[2.2	1.7	2.5	3.2	2.9	2.5];

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

for freqRatioInd= 6
    % for i=1:length(FluidVel) % i is for different velocity
    
    % load Free Vibration Data
    eval(['load ' freeVibDataDir char(newmatfiles(freqRatioInd)) '.mat;']);
    eval(['load ' freeVibDataDir char(phasematfiles(freqRatioInd)) '.mat;']);
    phase_ind = find(phasextoy > 180);
    phasextoy(phase_ind) = phasextoy(phase_ind)-360;
    eval(['load ' freeVibDataDir char(corrmatfiles(freqRatioInd)) '.mat;']);
    eval(['load ' freeVibDataDir char(ampfixmatfiles(freqRatioInd)) '.mat;']); 

    vr_free = Vrn.*fy./(ypeakfreq./2./pi);

    vr_fr(freqRatioInd,:)  = vr_free;
    Yad_fr(freqRatioInd,:) = Yad_fix;
    Xad_fr(freqRatioInd,:) = Xad_fix;
    Theta_fr(freqRatioInd,:) = phasextoy;
    vrn_fr(freqRatioInd,:) = Vrn;
    
    CL3_fr(freqRatioInd,:) = CLmag3;
    CL1_fr(freqRatioInd,:) = CLmag1;
    
    CL3_pred_free(freqRatioInd,i) = getDataPoint(ILCFHydroModelobj, 'CL3', ...
        Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i), Theta_fr(freqRatioInd,i)./180, vr_fr(freqRatioInd,i));
    CL1_pred_free(freqRatioInd,i) = getDataPoint(ILCFHydroModelobj, 'CL1', ...
        Yad_fr(freqRatioInd,i), Xad_fr(freqRatioInd,i), Theta_fr(freqRatioInd,i)./180, vr_fr(freqRatioInd,i));

    eval(['FluidVel(freqRatioInd,:) = vels' freqRatioNames{freqRatioInd} ';'])
% 
%     % begin adjust Fluid velocity to adjust Vrn for mass ratio change
%     FluidVel(freqRatioInd,:) = FluidVel(freqRatioInd,:) .* sqrt((massRatio_y_fr(freqRatioInd)+1)/(massRatio_y(freqRatioInd)+1));
%     % end added for mass ratio 1

%     % begin extend Vrn to 20    
%     FluidVel = [FluidVel FluidVel(end).*[13:20]./12];
%     % end extend Vrn to 20  
    
    %% step 1 Solve Vr,  use fsolve matlab function to solve Vry, Vrx, Theta
    for i= compIndx 
        i
        RigCylModel2Dobj(freqRatioInd,i) = RigCylModel2D(Diameter, Riserlength, massDensity_x(freqRatioInd), massDensity_y(freqRatioInd), ...
            FluidVel(freqRatioInd,i), fluidDensity, springCons_x(freqRatioInd), springCons_y(freqRatioInd), damperCons_x(freqRatioInd), damperCons_y(freqRatioInd));    
        vrn_pred(freqRatioInd,i) = RigCylModel2Dobj(freqRatioInd,i).FluidSpeed/(RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_y * RigCylModel2Dobj(freqRatioInd,i).Diameter);
        Vrn_y(freqRatioInd,i) = vrn_pred(freqRatioInd,i);
        Vrn_x(freqRatioInd,i) = RigCylModel2Dobj(freqRatioInd,i).FluidSpeed/(RigCylModel2Dobj(freqRatioInd,i).NominalNaturalFreq_x * RigCylModel2Dobj(freqRatioInd,i).Diameter);
        for j =1:numRun
            x0=zeros(1,5);
            x0(1) = Vrn_y(freqRatioInd,i);
            x0(2) = Vrn_x(freqRatioInd,i);
            x0(5) = 360./180.*(rand(1)-0.5);
            x0(3) = getDataPoint(ILCFHydroModelobj, 'Cmy', 1, 0.3, 0, x0(1));
            x0(4) = getDataPoint(ILCFHydroModelobj, 'Cmx', 1, 0.3, 0, x0(1));
                     
            options = optimset( 'Display', 'iter');

            f = @(x) myILCFVrfun(x, RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i), targetFreqRatio(freqRatioInd));

            lb =[3 3 -RigCylModel2Dobj(freqRatioInd,i).MassRatio_y  -RigCylModel2Dobj(freqRatioInd,i).MassRatio_x  -1];
            ub =[10 10 5 5 1];
             
            [x(freqRatioInd,i,j,:) res(freqRatioInd, i,j)] = lsqnonlin(f,x0,lb,ub);
        end 
    end
    for i= compIndx   
        casenums = [1:numRun];
        result_conv = reshape(x(freqRatioInd,i,casenums,:),length(casenums),5);
        res_conv = reshape(res(freqRatioInd,i,casenums,:),length(casenums),1);

        [result_uniq, indx] = uniquetol(result_conv,0.1,'rows');
        res_uniq = res_conv(indx);     
       % save solutions
        if ~isempty(result_uniq)
            theta2D_pred{i} = result_uniq(:,5).*180;
            vr2D_pred{i} = result_uniq(:,1);
            vr2D_pred_x{i} = result_uniq(:,2);
            Cmy_pred{i} = result_uniq(:,3);
            Cmx_pred{i} = result_uniq(:,4);
        else
            theta2D_pred{i} = NaN;
            vr2D_pred{i} =  NaN;
            vr2D_pred_x{i} =  NaN;
            Cmy_pred{i} =  NaN;
            Cmx_pred{i} =  NaN;
        end
        % Min Residual
        [res_min(i), res_min_ind(i)] = min(res_uniq);
        Theta_pred_res_min(freqRatioInd,i)  = theta2D_pred{i}(res_min_ind(i));
        vr_pred_res_min(freqRatioInd,i)     = vr2D_pred{i}(res_min_ind(i));
        vr_x_pred_res_min(freqRatioInd,i)   = vr2D_pred_x{i}(res_min_ind(i));
        
        vr2D_pred_res_min{i}    = vr2D_pred{i}(res_min_ind(i));
        vr2D_pred_x_res_min{i} = vr2D_pred_x{i}(res_min_ind(i));
        Cmy_pred_res_min{i} = Cmy_pred{i}(res_min_ind(i));
        Cmx_pred_res_min{i} = Cmx_pred{i}(res_min_ind(i));
        
        Yad_pred_res_min(freqRatioInd,i)    = NaN;
        Xad_pred_res_min(freqRatioInd,i)    = NaN;
   
        CLv_pred_res_min{i} = NaN;
        CDv_pred_res_min{i} =NaN;
        Damp_pred_y_res_min{i} = NaN;
        Damp_pred_x_res_min{i} = NaN;

        CL3_pred_res_min(freqRatioInd,i) = NaN;
        CL1_pred_res_min(freqRatioInd,i) = NaN;
    end
   % plot minimum residual solution   
    plot_oneSol(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, CL1_fr, CL3_fr, CL1_pred_free, CL3_pred_free, vrn_pred, vr_pred_res_min, Yad_pred_res_min, Xad_pred_res_min, Theta_pred_res_min, CL1_pred_res_min, CL3_pred_res_min, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   % Coefficeints verification for min res    
    plot_Coef(CLv_pred_res_min, CDv_pred_res_min, Cmy_pred_res_min, Cmx_pred_res_min, vrn_pred, ...
        CLv2, CDv2, Cmy, Cmx, vr_fr, Yad_fr, Xad_fr, Theta_fr, RigCylModel2Dobj(freqRatioInd, i), ILCFHydroModelobj, freqRatio, freqRatioInd, 9:24,  Vrn_y, Vrn_x, outputDir, 'minres')

        
%% step 1 Solve Vr,  use fsolve matlab function to solve Ay, Ax

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
        
%%    %  remove solutions with constraints
     for i= compIndx
        casenums = [1:numRun];

        result_conv = reshape(x(freqRatioInd,i,casenums,:),length(casenums),11);
        res_conv = reshape(res(freqRatioInd,i,casenums,:),length(casenums),1);
        
        % scale theta to same scaling as others
%         result_conv(:,3) = result_conv(:,3)./360;
        [result_uniq, indx] = uniquetol(result_conv,0.1,'rows');
        res_uniq = res_conv(indx);
%         result_uniq(:,3) = result_uniq(:,3).*360;     
% 
%        % remove solutions due to numerical erros
%        j=1;
%        f2 = @(x) myILCFfun2(x, RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i), targetFreqRatio(freqRatioInd));
%        while j <=size(result_uniq,1)
%             f0 = feval(f2, result_uniq(j,:))
%             if max(abs(f0(1:5)))>0.1
%                 result_uniq(j,:) = [];
%                 res_uniq(j) = [];
%             else
%                 j = j+1;
%             end
%        end
% 
%        f = @(x) myILCFfun(x, RigCylModel2Dobj(freqRatioInd,i), ILCFHydroModelobj, Vrn_y(freqRatioInd,i), Vrn_x(freqRatioInd,i), targetFreqRatio(freqRatioInd));
%        j=1;
%        while j <=size(result_uniq,1)
%             f0 = feval(f, result_uniq(j,:))
%             if max(abs(f0(1:5)))>0.1
%                 result_uniq(j,:) = [];
%                 res_uniq(j) = [];
%             else
%                 j = j+1;
%             end
%        end
%        
%       % remove solutions has Amp Y <0.15  for HZ  0.25 for JMD
%         
%        indxx = result_uniq(:,1)<0.15;  
%        result_uniq(indxx,:) =[];
%        res_uniq(indxx) =[];
%        
%        indxxx = result_uniq(:,2)<0.01;
%        result_uniq(indxxx,:) =[];
%        res_uniq(indxxx) =[];
%        
%      % remove solutions has Vr >8     
%        indx1 = find(result_uniq(:,4) > 8);  
%        temp = result_uniq;
%        temp(indx1,:) =[];
%        if ~isempty(temp)
%            result_uniq(indx1,:) =[];
%            res_uniq(indx1) =[];
%        end
%                      
%        % remove solutions has res > 1e-3
%        indx1 = find(res_uniq > 1e-3);
%        temp = result_uniq;
%        temp(indx1,:) =[];
%        if ~isempty(temp)       
%            result_uniq(indx1,:) =[];
%            res_uniq(indx1) =[];
%        end
%         
%        
%        % sort solutions by the Amp Y
%        [result_uniq indx4]= sortrows(result_uniq,-1);
%        res_uniq = res_uniq(indx4);
        
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

        % Min Residual
        [res_min(i), res_min_ind(i)] = min(res_uniq);
                   
        Yad_pred_res_min(freqRatioInd,i)    = Ampy_pred{i}(res_min_ind(i));
        Xad_pred_res_min(freqRatioInd,i)    = Ampx_pred{i}(res_min_ind(i));
        Theta_pred_res_min(freqRatioInd,i)  = theta2D_pred{i}(res_min_ind(i));
        vr_pred_res_min(freqRatioInd,i)     = vr2D_pred{i}(res_min_ind(i));
        vr_x_pred_res_min(freqRatioInd,i)   = vr2D_pred_x{i}(res_min_ind(i));

        vr2D_pred_res_min{i}    = vr2D_pred{i}(res_min_ind(i));
        vr2D_pred_x_res_min{i} = vr2D_pred_x{i}(res_min_ind(i));
        Cmy_pred_res_min{i} = Cmy_pred{i}(res_min_ind(i));
        Cmx_pred_res_min{i} = Cmx_pred{i}(res_min_ind(i));
        CLv_pred_res_min{i} = CLv_pred{i}(res_min_ind(i));
        CDv_pred_res_min{i} = CDv_pred{i}(res_min_ind(i));
        Damp_pred_y_res_min{i} = Damp_pred_y{i}(res_min_ind(i));
        Damp_pred_x_res_min{i} = Damp_pred_x{i}(res_min_ind(i));

        CL3_pred_res_min(freqRatioInd,i) = getDataPoint(ILCFHydroModelobj, 'CL3', ...
            Ampy_pred{i}(res_min_ind(i)), Ampx_pred{i}(res_min_ind(i)), theta2D_pred{i}(res_min_ind(i))./180, vr2D_pred{i}(res_min_ind(i)));
        CL1_pred_res_min(freqRatioInd,i) = getDataPoint(ILCFHydroModelobj, 'CL1', ...
            Ampy_pred{i}(res_min_ind(i)), Ampx_pred{i}(res_min_ind(i)), theta2D_pred{i}(res_min_ind(i))./180, vr2D_pred{i}(res_min_ind(i)));

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
        
        vr2D_pred_MaxAmpy{i}    = vr2D_pred{i}(MaxAmpy_ind(i));
        vr2D_pred_x_MaxAmpy{i} = vr2D_pred_x{i}(MaxAmpy_ind(i));
        Cmy_pred_MaxAmpy{i} = Cmy_pred{i}(MaxAmpy_ind(i));
        Cmx_pred_MaxAmpy{i} = Cmx_pred{i}(MaxAmpy_ind(i));
        CLv_pred_MaxAmpy{i} = CLv_pred{i}(MaxAmpy_ind(i));
        CDv_pred_MaxAmpy{i} = CDv_pred{i}(MaxAmpy_ind(i));
        Damp_pred_y_MaxAmpy{i} = Damp_pred_y{i}(MaxAmpy_ind(i));
        Damp_pred_x_MaxAmpy{i} = Damp_pred_x{i}(MaxAmpy_ind(i));

        CL3_pred_MaxAmpy(freqRatioInd,i) = getDataPoint(ILCFHydroModelobj, 'CL3', ...
            Ampy_pred{i}(MaxAmpy_ind(i)), Ampx_pred{i}(MaxAmpy_ind(i)), theta2D_pred{i}(MaxAmpy_ind(i))./180, vr2D_pred{i}(MaxAmpy_ind(i)));
        CL1_pred_MaxAmpy(freqRatioInd,i) = getDataPoint(ILCFHydroModelobj, 'CL1', ...
            Ampy_pred{i}(MaxAmpy_ind(i)), Ampx_pred{i}(MaxAmpy_ind(i)), theta2D_pred{i}(MaxAmpy_ind(i))./180, vr2D_pred{i}(MaxAmpy_ind(i)));

%         % HH
%         CL1_pred{i} = [];
%         CL3_pred{i} = [];
%         CL5_pred{i} = [];
%         CD2_pred{i} = [];
%         CD4_pred{i} = [];
%         for j = 1:length(Ampy_pred{i})
%            CL1_pred{i}(j) = getDataPoint(ILCFHydroModelobj, 'CL1', Ampy_pred{i}(j), Ampx_pred{i}(j), theta2D_pred{i}(j)./180, vr2D_pred{i}(j));
%            CL3_pred{i}(j) = getDataPoint(ILCFHydroModelobj, 'CL3', Ampy_pred{i}(j), Ampx_pred{i}(j), theta2D_pred{i}(j)./180, vr2D_pred{i}(j));
%            CL5_pred{i}(j) = getDataPoint(ILCFHydroModelobj, 'CL5', Ampy_pred{i}(j), Ampx_pred{i}(j), theta2D_pred{i}(j)./180, vr2D_pred{i}(j));
%            CD2_pred{i}(j) = getDataPoint(ILCFHydroModelobj, 'CD2', Ampy_pred{i}(j), Ampx_pred{i}(j), theta2D_pred{i}(j)./180, vr2D_pred{i}(j));
%            CD4_pred{i}(j) = getDataPoint(ILCFHydroModelobj, 'CD4', Ampy_pred{i}(j), Ampx_pred{i}(j), theta2D_pred{i}(j)./180, vr2D_pred{i}(j));
%         end
     end
      
   %% plot comparison

   % plot uncoupled solution
    plot_uncouple(vrn_fr, vr_fr, Yad_fr, Xad_fr, Vrn_y, Vr0_y,Vr0_x, NonDimAmp_linear0_y,NonDimAmp_linear0_x, freqRatio, freqRatioInd, compIndx, outputDir)

   % minimum residuals
    plot_min_res(vrn_fr, vr_fr, Yad_fr, Xad_fr, Theta_fr,vrn_pred, vr_pred_res_min, Yad_pred_res_min, Xad_pred_res_min, Theta_pred_res_min, x, res, numRun, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   % plot vr
    plot_vr(vrn_fr, vr_fr, vrn_pred, vr2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   % min average power
    plot_min_pow(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred, AvePow_pred, vr_pred_pow_min, Yad_pred_pow_min, Xad_pred_pow_min, Theta_pred_pow_min, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)
    close

   % plot max amp y solution   
    plot_maxY(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, CL1_fr, CL3_fr, CL1_pred_free, CL3_pred_free, vrn_pred, vr_pred_MaxAmpy, Yad_pred_MaxAmpy, Xad_pred_MaxAmpy, Theta_pred_MaxAmpy, CL1_pred_MaxAmpy, CL3_pred_MaxAmpy, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   % plot minimum residual solution   
    plot_oneSol(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, CL1_fr, CL3_fr, CL1_pred_free, CL3_pred_free, vrn_pred, vr_pred_res_min, Yad_pred_res_min, Xad_pred_res_min, Theta_pred_res_min, CL1_pred_res_min, CL3_pred_res_min, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   % Coefficeints verification for min res    
    plot_Coef(CLv_pred_res_min, CDv_pred_res_min, Cmy_pred_res_min, Cmx_pred_res_min, vrn_pred, ...
        CLv2, CDv2, Cmy, Cmx, vr_fr, Yad_fr, Xad_fr, Theta_fr, RigCylModel2Dobj(freqRatioInd, i), ILCFHydroModelobj, freqRatio, freqRatioInd, 9:24,  Vrn_y, Vrn_x, outputDir, 'minres')

   % Coefficeints verification for max amp y solution   
    plot_Coef(CLv_pred_MaxAmpy, CDv_pred_MaxAmpy, Cmy_pred_MaxAmpy, Cmx_pred_MaxAmpy, vrn_pred, ...
        CLv2, CDv2, Cmy, Cmx, vr_fr, Yad_fr, Xad_fr, Theta_fr, RigCylModel2Dobj(freqRatioInd, i), ILCFHydroModelobj, freqRatio, freqRatioInd, 9:24,  Vrn_y, Vrn_x, outputDir, 'MaxAmpy')
  
    
    
    
   % All solutions vs vrn
    plot_All_vrn(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred,vr2D_pred, Ampy_pred, Ampx_pred, theta2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   % All solutions vs vr
    plot_All_vr(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, vrn_pred,vr2D_pred, Ampy_pred, Ampx_pred, theta2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   % Added mass verification
    plot_addedMass(CLv_pred, CDv_pred, Cmy_pred, Cmx_pred, vr2D_pred_x, Vrn_y, Vrn_x, ...
    massRatio_y, massRatio_x, vrn_pred,vr2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

   % Individual Amplitude veification               
    plot_Amp_ver(CLv_pred, CDv_pred, Cmy_pred, Cmx_pred, vr2D_pred_x, Vrn_y, Vrn_x, massRatio_y, massRatio_x,...
    Damp_pred_y, Ampy_pred, Ampx_pred, Damp_pred_x, vr2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)
    close all

        % plot 3rd harmonic CL3
%         plot_hh_vrn(vrn_fr, CL3_fr, CL5_fr, vrn_pred,vr2D_pred, Ampy_pred, Ampx_pred, theta2D_pred, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)

        
    %    % plot combined max amp y solution
    %     plot_combine(vrn_fr, vr_fr,Yad_fr, Xad_fr, Theta_fr, Vrn_y, Vr0_y,Vr0_x, NonDimAmp_linear0_y,NonDimAmp_linear0_x, ...
    %         vrn_pred, vr_pred_res_min, Yad_pred_MaxAmpy, Xad_pred_MaxAmpy, Theta_pred_MaxAmpy, freqRatio, freqRatioInd, compIndx, col, shp, outputDir)
    %  
end
     
save ICFLMotion_ZDMiT_RobustSmooth_20130908_run10.mat