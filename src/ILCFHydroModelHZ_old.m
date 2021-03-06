classdef ILCFHydroModelHZ_old
    % Use Zeros to fill the empty exp data instead of NaNs
    % load data directly
   
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
        bigCmx_smooth
        bigCmy_smooth
        bigCLv_smooth
        bigCDv_smooth
        bigpow_smooth
    end
    
    methods
         function obj =  ILCFHydroModelHZ_old()
       
            obj.vr_set = {'4' '4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
            obj.theta_set = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
                 '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
            obj.X_set = {'0.05','0.1','0.15','0.2','0.25','0.3'};
            obj.Y_set = {'0.15','0.25','0.5','0.75', '1'};
            obj.phase_vec = [-180:22.5:-45, -30:15:30, 45:22.5:180]./180;
	        phase_vec_temp = [-360:15:-330, -315:22.5:-45, -30:15:30, 45:22.5:315, 330:15:345];
            obj.X_vec = [0.05:0.05:0.3];
            obj.Y_vec = [0.15,0.25,0.5,0.75,1];
            obj.vr_vec = [4:0.5:8];
            compType = computer;
            if strcmp(compType, 'PCWIN') == 1
                dataPath =  'D:\Dropbox\VIV\src\ILCFprediction\databaseHZ\';
            else
                dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseHZ/';
            end
            expName = 'haining';
            for p = 1:length(obj.vr_set)
                eval(['load ' dataPath 'vr' char(obj.vr_set(p)) '.mat']);

                [ind1,ind2] = find(isnan(Cmx));
                Cmx(ind1,ind2) = 1;%% just for ploting reason

                for i = 1:(length(obj.theta_set))
                  for j = 1:length(obj.Y_vec)
                     ind = find(allYad(i,:) == obj.Y_vec(j));
                     for k = ind
                         newX(j,1:length(ind),i+9) = allXad(i,ind);% 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
                         newY(j,1:length(ind),i+9) = allYad(i,ind);% 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
                         newtheta(j,1:length(ind),i+9) = alltheta(i,ind)./180;% theta rotating, 8 Phase, from -180 to 180 degrees in increments of 45 degrees.
%                          newtheta(j,1:length(ind),i) = alltheta(i,ind);
                         newvr(j,1:length(ind),i+9) = Vr(i,ind);
                         newCmy(j,1:length(ind),i+9) = Cmy(i,ind);
                         newCmx(j,1:length(ind),i+9) = Cmx(i,ind);
                         newCLv(j,1:length(ind),i+9) = CLv(i,ind);
                         newCDv(j,1:length(ind),i+9) = CDv(i,ind);
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
                                avgPow(j,1:length(ind),i+9) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.18.^3*0.6858*0.0381);
                            else
                                avgPow(j,1:length(ind),i+9) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.23.^3*0.6858*0.0381);
                            end
                        elseif strcmp(expName, 'haining')
                            avgPow(j,1:length(ind),i+9) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.2.^3*0.64135*0.0381);
                        end
                     end
                  end    
                end
                
                for i = 1:9
                    newX(:,:,i) = newX(:,:,i+18);  % 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
                    newY(:,:,i) = newY(:,:,i+18);  % 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
                    newtheta(:,:,i) = phase_vec_temp(i)./180;      % theta rotating? 8 Phase, from -180 to 180 degrees in increments of 45 degrees.
                    newvr(:,:,i) = newvr(:,:,i+18);      % 8 reduced velocity, from 4 to 8 in increments of 0.5 
                    newCmy(:,:,i) = newCmy(:,:,i+18);
                    newCmx(:,:,i) = newCmx(:,:,i+18);
                    newCDv(:,:,i) = newCDv(:,:,i+18);
                    newCLv(:,:,i) = newCLv(:,:,i+18);
                    avgPow(:,:,i) = avgPow(:,:,i+18);
                end
                for i = 28:36
                    newX(:,:,i) = newX(:,:,i-18);  % 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
                    newY(:,:,i) = newY(:,:,i-18);  % 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
                    newtheta(:,:,i) = phase_vec_temp(i)./180;      % theta rotating? 8 Phase, from -180 to 180 degrees in increments of 45 degrees.
                    newvr(:,:,i) = newvr(:,:,i-18);      % 8 reduced velocity, from 4 to 8 in increments of 0.5 
                    newCmy(:,:,i) = newCmy(:,:,i-18);
                    newCmx(:,:,i) = newCmx(:,:,i-18);
                    newCDv(:,:,i) = newCDv(:,:,i-18);
                    newCLv(:,:,i) = newCLv(:,:,i-18);
                    avgPow(:,:,i) = avgPow(:,:,i-18);
                end
                 
                obj.bigX(:,:,:,p) = newX;
                obj.bigY(:,:,:,p) = newY;
                obj.bigtheta(:,:,:,p) = newtheta;
                obj.bigvr(:,:,:,p) = newvr;
                obj.bigCmy(:,:,:,p) = -newCmy;
                obj.bigCmx(:,:,:,p) = newCmx;
                obj.bigCLv(:,:,:,p) = -newCLv;
                obj.bigCDv(:,:,:,p) = newCDv;
                obj.bigpow(:,:,:,p) = avgPow;
            end
            obj.bigCDv(isnan(obj.bigCDv)) = 0 ;
            
            [bigCmx_smooth_temp, s, exitflag] = smoothn(obj.bigCmx,'robust');
            [bigCmy_smooth_temp, s, exitflag] = smoothn(obj.bigCmy,'robust');
            [bigCDv_smooth_temp, s, exitflag] = smoothn(obj.bigCDv,'robust');
            [bigCLv_smooth_temp, s, exitflag] = smoothn(obj.bigCLv,'robust');
            [bigpow_smooth_temp, s, exitflag] = smoothn(obj.bigpow,'robust');
            
            obj.bigX = obj.bigX(:,:,10:28,:);
            obj.bigY = obj.bigY(:,:,10:28,:);
            obj.bigtheta = obj.bigtheta(:,:,10:28,:);
            obj.bigvr = obj.bigvr(:,:,10:28,:);
            obj.bigCmx = obj.bigCmx(:,:,10:28,:);
            obj.bigCmy = obj.bigCmy(:,:,10:28,:);
            obj.bigCDv = obj.bigCDv(:,:,10:28,:);
            obj.bigCLv = obj.bigCLv(:,:,10:28,:);
            obj.bigpow = obj.bigpow(:,:,10:28,:);
            obj.bigCmx_smooth = bigCmx_smooth_temp(:,:,10:28,:);
            obj.bigCmy_smooth = bigCmy_smooth_temp(:,:,10:28,:);
            obj.bigCDv_smooth = bigCDv_smooth_temp(:,:,10:28,:);
            obj.bigCLv_smooth = bigCLv_smooth_temp(:,:,10:28,:);
            obj.bigpow_smooth = bigpow_smooth_temp(:,:,10:28,:);
        end
          
         function [value] = getDataPoint(obj, targetpar, y, x, theta, vr)
%              int = 'spline';
             int = 'linear';
            if abs(theta)>10
                display('Theta Input Format Error')
            end
%             
%             if y >1.5
%                 y=1.5;
%             end
%             if x >0.75
%                 x= 0.75;
%             end 
%             
            if vr <4
                vr =4;
            elseif vr>8
                vr =8;
            end
%             
%             if theta<-1
%                 theta = theta+2;
%             elseif theta >1
%                 theta = theta-2;
%             end
               
             if strcmp(targetpar, 'CLv')
%                  if vr <4
%                      value= 0;
%                  elseif vr>8
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCLv_smooth,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'CDv')
%                  if vr <4
%                      value= 0;
%                  elseif vr>8
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCDv_smooth,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'Cmx')
%                  if vr <4
%                     value= 1;
%                  elseif vr>8
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmx_smooth,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'Cmy')
%                  if vr <4
%                     value= 1;
%                  elseif vr>8
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmy_smooth,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'pow')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigpow_smooth,y, x, theta, vr,int,0);
             end
         end
          
         function [value] = getDataPoint_ori(obj, targetpar, y, x, theta, vr)
%              int = 'spline';
%              int = 'linear';
             int = 'nearest';
            if abs(theta)>10
                display('Theta Input Format Error')
            end
%             
%             if y >1.5
%                 y=1.5;
%             end
%             if x >0.75
%                 x= 0.75;
%             end 
%             
            if vr <4
                vr =4;
            elseif vr>8
                vr =8;
            end
%             
%             if theta<-1
%                 theta = theta+2;
%             elseif theta >1
%                 theta = theta-2;
%             end
               
             if strcmp(targetpar, 'CLv')
%                  if vr <4
%                      value= 0;
%                  elseif vr>8
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCLv,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'CDv')
%                  if vr <4
%                      value= 0;
%                  elseif vr>8
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCDv,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'Cmx')
%                  if vr <4
%                     value= 1;
%                  elseif vr>8
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmx,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'Cmy')
%                  if vr <4
%                     value= 1;
%                  elseif vr>8
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmy,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'pow')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigpow,y, x, theta, vr,int,0);
             end
         end
         
         function plot3D_theta(obj, targetpar, y, x, theta, vr)
                    
             color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
             for i =1:length(vr)
                 for j = 1:length(x)
                     for k =  1:length(y)
                        targetValue(k,j,i) = getDataPoint(obj, targetpar, y(k), x(j), theta/180, vr(i));
                     end
                 end
             end
            
            pX(:,:,:) = obj.bigX(:,:,:,1);
            pY(:,:,:) = obj.bigY(:,:,:,1);
            pvr(:,:,:)= obj.bigvr(:,:,:,1);

            minP = min(min(min(targetValue)));
            maxP = max(max(max(targetValue)));;                    
            isoRange = linspace(0, maxP-minP,5);
            isoRange = minP + isoRange(2:4);
            isoRange_set = cell(1,3);
            for m = 1:3
                 isoRange_set{m} = num2str(isoRange(m),2);
            end                      
            for j = 1:length(isoRange)
                figh(i,j) = patch(isosurface(y,x,vr, targetValue, isoRange(j) ));           
                set(figh(i,j),'FaceColor',color_set{j});
                set(figh(i,j),'FaceAlpha',0.3)
                ylabel('Inline Motion')
                xlabel('Crossflow Motion')
                zlabel('Reduced Velocity')
                hold on
            end
            legend(isoRange_set)
            title([targetpar ' at theta ' num2str(theta) ])
            grid on
            view(-20,30)
            hold off
         end
         
         function plot2D_vr(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end
        
             figure
             surface(vr, theta, targetValue)
             xlabel('Reduced Velocity')
             ylabel('Theta')
             view([65 30])
             title([targetpar ' vs ' 'Vr' ' @ x=' num2str(x) ', y=' num2str(y)])
             
             figure
             surface(1./vr, theta, targetValue)
             xlabel('Non-Dimensional Frequency')
             ylabel('Theta')
             view([65 30])
             title([targetpar ' vs ' 'f*' ' @ x=' num2str(x) ', y=' num2str(y)])
         end
         
         function h = plot1D_f(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                targetValue(i) = getDataPoint(obj, targetpar, y, x, theta/180, vr(i));
             end
%              figure
             h = plot(1./vr, targetValue,'-o');
             xlabel('f*')
             xlim([0.1 0.25])
%              title([targetpar ' vs ' 'f*' ' @ x=' num2str(x) ', theta=' num2str(theta) ', y=' num2str(y)])
%              title([ 'theta=' num2str(theta) ])
         end
         
         function h = plot1D_vr(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                targetValue(i) = getDataPoint(obj, targetpar, y, x, theta/180, vr(i));
             end
%              figure
             h = plot(vr, targetValue,'-o');
             xlabel('vr')
%              title([targetpar ' vs ' 'vr' ' @ x=' num2str(x) ', theta=' num2str(theta) ', y=' num2str(y)])
             title([ 'theta=' num2str(theta) ])
         end
         
         function plot1D_theta(obj, targetpar, y, x, theta, vr)
             for i =1:length(theta)
                targetValue(i) = getDataPoint(obj, targetpar, y, x, theta(i)/180, vr);
             end
             figure
             plot(theta*180, targetValue,'-o')
             xlabel('theta')
             title([targetpar ' vs ' 'theta' ' @ x=' num2str(x) ', y=' num2str(y) ', vr=' num2str(vr)])
         end
         function plot1D_y(obj, targetpar, y, x, theta, vr)
             for i =1:length(y)
                targetValue(i) = getDataPoint(obj, targetpar, y(i), x, theta/180, vr);
             end
             figure
             plot(y, targetValue,'-o')
             xlabel('y')
             title([targetpar ' vs ' 'y' ' @ x=' num2str(x) ', theta=' num2str(theta*180,3) ', vr=' num2str(vr)])
         end
         function plot1D_x(obj, targetpar, y, x, theta, vr)
             for i =1:length(x)
                targetValue(i) = getDataPoint(obj, targetpar, y, x(i), theta/180, vr);
             end
             figure
             plot(x, targetValue,'-o')
             xlabel('x')
             title([targetpar ' vs ' 'x' ' @ y=' num2str(y) ', theta=' num2str(theta*180,3) ', vr=' num2str(vr)])
         end
         
         function plot2D_vr_ori(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint_ori(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end
        
             figure
             surface(vr, theta, targetValue)
             xlabel('Reduced Velocity')
             ylabel('Theta')
             view([65 30])
             title([targetpar ' vs ' 'Vr' ' @ x=' num2str(x) ', y=' num2str(y)])
             
             figure
             surface(1./vr, theta, targetValue)
             xlabel('Non-Dimensional Frequency')
             ylabel('Theta')
             view([65 30])
             title([targetpar ' vs ' 'f*' ' @ x=' num2str(x) ', y=' num2str(y)])
         end
         
         function h = plot1D_f_ori(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                targetValue(i) = getDataPoint_ori(obj, targetpar, y, x, theta/180, vr(i));
             end
%              figure
             h = plot(1./vr, targetValue,'-+r');
             xlabel('f*')
             xlim([0.1 0.25])
%              title([targetpar ' vs ' 'f*' ' @ x=' num2str(x) ', theta=' num2str(theta) ', y=' num2str(y)])
%              title([ 'theta=' num2str(theta) ])
         end         
         function h = plot1D_vr_ori(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                targetValue(i) = getDataPoint_ori(obj, targetpar, y, x, theta/180, vr(i));
             end
%              figure
             h = plot(vr, targetValue,'-+r');
             xlabel('vr')
%              title([targetpar ' vs ' 'vr' ' @ x=' num2str(x) ', theta=' num2str(theta) ', y=' num2str(y)])
             title([ 'theta=' num2str(theta) ])
         end
         
         function h =  plot1D_theta_ori(obj, targetpar, y, x, theta, vr)
             for i =1:length(theta)
                targetValue(i) = getDataPoint_ori(obj, targetpar, y, x, theta(i)/180, vr);
             end
%              figure
              h = plot(theta, targetValue,'-+r');
%              xlabel('theta')
%              title([targetpar ' vs ' 'theta' ' @ x=' num2str(x) ', y=' num2str(y) ', vr=' num2str(vr)])
         end
         function h =  plot1D_y_ori(obj, targetpar, y, x, theta, vr)
             for i =1:length(y)
                targetValue(i) = getDataPoint_ori(obj, targetpar, y(i), x, theta/180, vr);
             end
%              figure
              h = plot(y, targetValue,'-+r');
%              xlabel('y')
%              title([targetpar ' vs ' 'y' ' @ x=' num2str(x) ', theta=' num2str(theta*180,3) ', vr=' num2str(vr)])
         end
         function  h = plot1D_x_ori(obj, targetpar, y, x, theta, vr)
             for i =1:length(x)
                targetValue(i) = getDataPoint_ori(obj, targetpar, y, x(i), theta/180, vr);
             end
%              figure
              h = plot(x, targetValue,'-+r');
             xlabel('x')
             title([targetpar ' vs ' 'x' ' @ y=' num2str(y) ', theta=' num2str(theta*180,3) ', vr=' num2str(vr)])
         end 

    end
end