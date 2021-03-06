classdef ILCFHydroModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

   
    properties
        vr_set    
        theta_set
        X_set
        Y_set
        phase_vec
        X_vec
        Y_vec
        vr_vec
        bigX
        bigY
        bigtheta
        bigvr
        bigCmx
        bigCmy
        bigCLv
        bigCDv
        bigpow
    end
    
    methods
         function obj =  ILCFHydroModel()
            
            obj.vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
            obj.theta_set = {'-pi','-3/4pi','-1/2pi','-1/4pi','0','1/4pi','1/2pi','3/4pi'};
            obj.X_set = {'0','0.15','0.3','0.45','0.6','0.75'};
            obj.Y_set = {'0.25','0.5','0.75','1','1.25','1.5'};
            obj.phase_vec = [-180 -135 -90 -45 0 45 90 135 180];
            obj.X_vec = [0:0.15:0.75];
            obj.Y_vec = [0.25:0.25:1.5];
            obj.vr_vec = [4.5:0.5:8];
            dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseJMD/';
            expName = 'jason';
            for p = 1:length(obj.vr_set)
                eval(['load ' dataPath 'vr' char(obj.vr_set(p)) '.mat']);

                [ind1,ind2] = find(isnan(Cmx));
                Cmx(ind1,ind2) = 1;%% just for ploting reason

                for i = 1:(length(obj.theta_set))
                  for j = 1:length(obj.Y_vec)
                     ind = find(allYad(i,:) == obj.Y_vec(j));
                     for k = ind
                         newX(j,1:length(ind),i) = allXad(i,ind);% 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
                         newY(j,1:length(ind),i) = allYad(i,ind);% 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
                         newtheta(j,1:length(ind),i) = alltheta(i,ind);% theta rotating, 8 Phase, from -180 to 180 degrees in increments of 45 degrees.
                         newvr(j,1:length(ind),i) = Vr(i,ind);
                         newCmy(j,1:length(ind),i) = Cmy(i,ind);
                         newCmx(j,1:length(ind),i) = Cmx(i,ind);
                         newCLv(j,1:length(ind),i) = CLv(i,ind);
                         newCDv(j,1:length(ind),i) = CDv(i,ind);
                         newCDv(j,1:length(ind),i) = CDv(i,ind);
                         if strcmp(expName, 'jason')
                %%             special for jason's experiment
                            if p == 1 | p == 2 | p == 3 % for vr less than 6
                                % The Reynolds number is held constant at a value of 8760 for reduced velocities greater than 
                                % or equal to 6 and the Reynolds number is 6860 for reduced velocities less than 6. The Reynolds
                                % number cannot be held constant for all motions due to frequency constraints on the experimental apparatus.

                                % The diameter of the cylinder was 0.0381 m (1.5 inches) in order to maximize Reynolds number 
                                % and minimize the forcing frequency for forced motions. The span was 0.6858 m (27 inches).
                                % current vel :0.18m/s =>Re=6858
                                % current vel :0.23m/s =>Re=8763
                                avgPow(j,1:length(ind),i) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.18.^3*0.6858*0.0381);
                            else
                                avgPow(j,1:length(ind),i) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.23.^3*0.6858*0.0381);
                            end
                        elseif strcmp(expName, 'haining')
                            avgPow(j,1:length(ind),i) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.2.^3*0.64135*0.0381);
                        end
                     end
                  end    
                end
                
                newX(:,:,9) = newX(:,:,1);  % 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
                newY(:,:,9) = newY(:,:,1);  % 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
                newtheta(:,:,9) = 180;      % theta rotating? 8 Phase, from -180 to 180 degrees in increments of 45 degrees.
                newvr(:,:,9) = newvr(:,:,1);      % 8 reduced velocity, from 4.5 to 8 in increments of 0.5 
                newCmy(:,:,9) = newCmy(:,:,1);
                newCmx(:,:,9) = newCmx(:,:,1);
                newCDv(:,:,9) = newCDv(:,:,1);
                newCLv(:,:,9) = newCLv(:,:,1);
                avgPow(:,:,9) = avgPow(:,:,1);
                 
                obj.bigX(:,:,:,p) = newX;
                obj.bigY(:,:,:,p) = newY;
                obj.bigtheta(:,:,:,p) = newtheta;
                obj.bigvr(:,:,:,p) = newvr;
                obj.bigCmy(:,:,:,p) = newCmy;
                obj.bigCmx(:,:,:,p) = newCmx;
                obj.bigCLv(:,:,:,p) = newCLv;
                obj.bigCDv(:,:,:,p) = newCDv;
                obj.bigpow(:,:,:,p) = avgPow;
            end
            obj.bigCDv(isnan(obj.bigCDv)) = 0 

        end
          
         function [value] = getDataPoint(obj, targetpar, y, x, theta, vr)
             int = 'spline';
             if strcmp(targetpar, 'CLv')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCLv,y, x, theta, vr,int);
             elseif strcmp(targetpar, 'CDv')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCDv,y, x, theta, vr,int);
             elseif strcmp(targetpar, 'Cmx')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmx,y, x, theta, vr,int);
             elseif strcmp(targetpar, 'Cmy')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmy,y, x, theta, vr,int);
             elseif strcmp(targetpar, 'pow')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigpow,y, x, theta, vr,int);
             end
         end
         
         function plot1D_vr(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                targetValue(i) = getDataPoint(obj, targetpar, y, x, theta, vr(i));
             end
             figure
             plot(vr, targetValue,'-o')
             xlabel('vr')
             title([targetpar ' vs ' 'vr' ' @ x=' num2str(x) ', theta=' num2str(theta) ', y=' num2str(y)])
         end
         function plot1D_theta(obj, targetpar, y, x, theta, vr)
             for i =1:length(theta)
                targetValue(i) = getDataPoint(obj, targetpar, y, x, theta(i), vr);
             end
             figure
             plot(theta, targetValue,'-o')
             xlabel('theta')
             title([targetpar ' vs ' 'theta' ' @ x=' num2str(x) ', y=' num2str(y) ', vr=' num2str(vr)])
         end
         function plot1D_y(obj, targetpar, y, x, theta, vr)
             for i =1:length(y)
                targetValue(i) = getDataPoint(obj, targetpar, y(i), x, theta, vr);
             end
             figure
             plot(y, targetValue,'-o')
             xlabel('y')
             title([targetpar ' vs ' 'y' ' @ x=' num2str(x) ', theta=' num2str(theta) ', vr=' num2str(vr)])
         end
         function plot1D_x(obj, targetpar, y, x, theta, vr)
             for i =1:length(x)
                targetValue(i) = getDataPoint(obj, targetpar, y, x(i), theta, vr);
             end
             figure
             plot(x, targetValue,'-o')
             xlabel('x')
             title([targetpar ' vs ' 'x' ' @ y=' num2str(y) ', theta=' num2str(theta) ', vr=' num2str(vr)])
         end                  
    end
end