classdef ILCFHydroModelZDMiT_NormalSmooth
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
        bigCmy_int
        bigCmx_int
        bigCLv_int
        bigCDv_int
        bigpow_int
        bigCmx_smooth
        bigCmy_smooth
        bigCLv_smooth
        bigCDv_smooth
        bigpow_smooth
    end
    
    methods
         function obj =  ILCFHydroModelZDMiT_NormalSmooth()

            obj.vr_set = {'3' '3.5' '4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8' '9' '10' '15' '20' '40' '100'}; 
            obj.X_set = {'0','0.05','0.1','0.15','0.2','0.25','0.3','0.45','0.6','0.75'};
            obj.Y_set = {'0.15','0.25','0.5','0.75', '1', '1.25', '1.5'};

            obj.theta_set = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
                 '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
            obj.Y_vec = [0.15, 0.25, 0.5, 0.75, 1, 1.25, 1.5 ];
            obj.X_vec = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.45,0.6,0.75 ];
            obj.phase_vec = [-180:22.5:-45, -30:15:30, 45:22.5:180]./180;
            obj.vr_vec = [3:0.5:8,9,10,15,20,40,100 ];
            phase_vec_temp = [-360:15:-330, -315:22.5:-45, -30:15:30, 45:22.5:315, 330:15:345]./180;

            [rawDatabaseStrucHZ] = loadDatabaseHZ();

            [rawDatabaseStrucJMD] = loadDatabaseJMD();

            rawDatabaseStruc = combineDatabase(rawDatabaseStrucHZ, rawDatabaseStrucJMD);           

            obj.bigX = rawDatabaseStruc.X;
            obj.bigY = rawDatabaseStruc.Y;
            obj.bigtheta = rawDatabaseStruc.Theta;
            obj.bigvr = rawDatabaseStruc.Vr;
            obj.bigCmy = rawDatabaseStruc.Cmy;
            obj.bigCmx = rawDatabaseStruc.Cmx;
            obj.bigCLv = rawDatabaseStruc.CLv;
            obj.bigCDv = rawDatabaseStruc.CDv;
            obj.bigpow = rawDatabaseStruc.Pow;
            
            obj.bigCmy_int = inpaintn(obj.bigCmy);
            obj.bigCmx_int = inpaintn(obj.bigCmx);
            obj.bigCLv_int = inpaintn(obj.bigCLv);
            obj.bigCDv_int = inpaintn(obj.bigCDv);
            obj.bigpow_int = inpaintn(obj.bigpow);
            
%             [bigCmx_smooth_temp, s, exitflag] = smoothn(obj.bigCmx_int,'robust');
%             [bigCmy_smooth_temp, s, exitflag] = smoothn(obj.bigCmy_int,'robust');
%             [bigCDv_smooth_temp, s, exitflag] = smoothn(obj.bigCDv_int,'robust');
%             [bigCLv_smooth_temp, s, exitflag] = smoothn(obj.bigCLv_int,'robust');
%             [bigpow_smooth_temp, s, exitflag] = smoothn(obj.bigpow_int,'robust');
            
            smoPar = 0.1;
            [bigCmx_smooth_temp, s, exitflag] = smoothn(obj.bigCmx_int,smoPar);
            [bigCmy_smooth_temp, s, exitflag] = smoothn(obj.bigCmy_int,smoPar);
            [bigCDv_smooth_temp, s, exitflag] = smoothn(obj.bigCDv_int,smoPar);
            [bigCLv_smooth_temp, s, exitflag] = smoothn(obj.bigCLv_int,smoPar);
            [bigpow_smooth_temp, s, exitflag] = smoothn(obj.bigpow_int,smoPar);
            
            
%             [bigCmx_smooth_temp, s, exitflag] = smoothn(obj.bigCmx,'robust');
%             [bigCmy_smooth_temp, s, exitflag] = smoothn(obj.bigCmy,'robust');
%             [bigCDv_smooth_temp, s, exitflag] = smoothn(obj.bigCDv,'robust');
%             [bigCLv_smooth_temp, s, exitflag] = smoothn(obj.bigCLv,'robust');
%             [bigpow_smooth_temp, s, exitflag] = smoothn(obj.bigpow,'robust');
%             
            obj.bigX = obj.bigX(:,:,10:28,:);
            obj.bigY = obj.bigY(:,:,10:28,:);
            obj.bigtheta = obj.bigtheta(:,:,10:28,:);
            obj.bigvr = obj.bigvr(:,:,10:28,:);
            
            obj.bigCmx = obj.bigCmx(:,:,10:28,:);
            obj.bigCmy = obj.bigCmy(:,:,10:28,:);
            obj.bigCDv = obj.bigCDv(:,:,10:28,:);
            obj.bigCLv = obj.bigCLv(:,:,10:28,:);
            obj.bigpow = obj.bigpow(:,:,10:28,:);
            
            obj.bigCmx_int = obj.bigCmx_int(:,:,10:28,:);
            obj.bigCmy_int = obj.bigCmy_int(:,:,10:28,:);
            obj.bigCDv_int = obj.bigCDv_int(:,:,10:28,:);
            obj.bigCLv_int = obj.bigCLv_int(:,:,10:28,:);
            obj.bigpow_int = obj.bigpow_int(:,:,10:28,:);
            
            obj.bigCmx_smooth = bigCmx_smooth_temp(:,:,10:28,:);
            obj.bigCmy_smooth = bigCmy_smooth_temp(:,:,10:28,:);
            obj.bigCDv_smooth = bigCDv_smooth_temp(:,:,10:28,:);
            obj.bigCLv_smooth = bigCLv_smooth_temp(:,:,10:28,:);
            obj.bigpow_smooth = bigpow_smooth_temp(:,:,10:28,:);
            
                  
            % make sure phase has no influence when x = 0
            for i = 1:length(obj.Y_vec)
                for j = 1:1:length(obj.vr_vec)
                    obj.bigCmy(i,1,:,j) = mean(obj.bigCmy(i,1,1:18,j));
                    obj.bigCmx(i,1,:,j) = mean(obj.bigCmx(i,1,1:18,j));
                    obj.bigCLv(i,1,:,j) = mean(obj.bigCLv(i,1,1:18,j));
                    obj.bigCDv(i,1,:,j) = mean(obj.bigCDv(i,1,1:18,j));
                    obj.bigpow(i,1,:,j) = mean(obj.bigpow(i,1,1:18,j));
                    obj.bigCmy_int(i,1,:,j) = mean(obj.bigCmy_int(i,1,1:18,j));
                    obj.bigCmx_int(i,1,:,j) = mean(obj.bigCmx_int(i,1,1:18,j));
                    obj.bigCLv_int(i,1,:,j) = mean(obj.bigCLv_int(i,1,1:18,j));
                    obj.bigCDv_int(i,1,:,j) = mean(obj.bigCDv_int(i,1,1:18,j));
                    obj.bigpow_int(i,1,:,j) = mean(obj.bigpow_int(i,1,1:18,j));
                    obj.bigCmy_smooth(i,1,:,j) = mean(obj.bigCmy_int(i,1,1:18,j));
                    obj.bigCmx_smooth(i,1,:,j) = mean(obj.bigCmx_int(i,1,1:18,j));
                    obj.bigCLv_smooth(i,1,:,j) = mean(obj.bigCLv_int(i,1,1:18,j));
                    obj.bigCDv_smooth(i,1,:,j) = mean(obj.bigCDv_int(i,1,1:18,j));
                    obj.bigpow_smooth(i,1,:,j) = mean(obj.bigpow_int(i,1,1:18,j));
                end
            end
            
             %%  force values for low and high Vr 
             for i = 1:length(obj.Y_set)
                for j= 1:length(obj.X_set)
                    for m = 1:length(obj.theta_set)
                       obj.bigCmy_smooth(i,j,m,1:3) = linspace(1,obj.bigCmy_smooth(i,j,m,3),3);
                       obj.bigCmx_smooth(i,j,m,1:3) = linspace(1,obj.bigCmx_smooth(i,j,m,3),3);
                       obj.bigCLv_smooth(i,j,m,1:3) = linspace(0,obj.bigCLv_smooth(i,j,m,3),3);
                       obj.bigCDv_smooth(i,j,m,1:3) = linspace(0,obj.bigCDv_smooth(i,j,m,3),3);
                       obj.bigpow_smooth(i,j,m,1:3) = linspace(0,obj.bigpow_smooth(i,j,m,3),3);   
                       obj.bigCmy_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigCmy_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigCmy_smooth(i,j,:,length(obj.vr_set)-6)),7);
                       obj.bigCmx_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigCmx_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigCmx_smooth(i,j,:,length(obj.vr_set)-6)),7);
                       obj.bigCLv_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigCLv_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigCLv_smooth(i,j,:,length(obj.vr_set)-6)),7);
                       obj.bigCDv_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigCDv_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigCDv_smooth(i,j,:,length(obj.vr_set)-6)),7);
                       obj.bigpow_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigpow_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigpow_smooth(i,j,:,length(obj.vr_set)-6)),7);
                    end
                end
             end

%              %% force values for low and high Vr 
%              for i = 1:length(obj.Y_set)
%                 for j= 1:length(obj.X_set)
%                     for m = 1:length(obj.theta_set)
%                         for n = 1:4
%                            obj.bigCmy_smooth(i,j,m,n) = 1;
%                            obj.bigCmx_smooth(i,j,m,n) = 1;
%                            obj.bigCLv_smooth(i,j,m,n) = 0;
%                            obj.bigCDv_smooth(i,j,m,n) = 0;
%                            obj.bigpow_smooth(i,j,m,n) = 0;
%                         end
%                         for n = (length(obj.vr_set)-3): length(obj.vr_set)
%                            obj.bigCmy_smooth(i,j,m,n) = -0.5;
%                            obj.bigCmx_smooth(i,j,m,n) = -0.5;
%                            obj.bigCLv_smooth(i,j,m,n) = 0;
%                            obj.bigCDv_smooth(i,j,m,n) = 0;
%                            obj.bigpow_smooth(i,j,m,n) = 0;     
%                         end
%                     end
%                 end
%             end
   
            
        end
          
         function [value] = getDataPoint(obj, targetpar, y, x, theta, vr)
%              int = 'spline';
             int = 'linear';
            if abs(theta)>10
                display('Theta Input Format Error')
            end            
             
%             if y >1.5
%                 y=1.5;
%             end
%             if x >0.75
%                 x= 0.75;
%             end 
%             
            if vr <3
                vr = 3;
            elseif vr>100
                vr =100;
            end
%             
%             if theta<-1
%                 theta = theta+2;
%             elseif theta >1
%                 theta = theta-2;
%             end
               
             if strcmp(targetpar, 'CLv')
%                  if vr <3
%                      value= 0;
%                  elseif vr>100
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCLv_smooth,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'CDv')
%                  if vr <3
%                      value= 0;
%                  elseif vr>100
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCDv_smooth,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'Cmx')
%                     % begin, fix theta to 45 degree
%                     theta = 45/180;
%                     % end, fix theta to 45 degree
% 
%                     % begin, fix y,x to 0.75, 0.15 
%                     y = 0.75;
%                     x = 0.15;
%                     % end, fix y,x to 0.75, 0.15 
                    
%                  if vr <3
%                     value= 1;
%                  elseif vr>100
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmx_smooth,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'Cmy')
                 
%                     % begin, fix theta to 45 degree
%                     theta = 45/180;
%                     % end, fix theta to 45 degree
% 
%                     % begin, fix y,x to 0.75, 0.15 
%                     y = 0.75;
%                     x = 0.15;
%                     % end, fix y,x to 0.75, 0.15 
                    
%                  if vr <3
%                     value= 1;
%                  elseif vr>100
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
            if vr <3
                vr = 3;
            elseif vr>100
                vr =100;
            end
%             
%             if theta<-1
%                 theta = theta+2;
%             elseif theta >1
%                 theta = theta-2;
%             end
               
             if strcmp(targetpar, 'CLv')
%                  if vr <3
%                      value= 0;
%                  elseif vr>100
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCLv,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'CDv')
%                  if vr <3
%                      value= 0;
%                  elseif vr>100
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCDv,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'Cmx')
%                  if vr <3
%                     value= 1;
%                  elseif vr>100
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmx,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'Cmy')
%                  if vr <3
%                     value= 1;
%                  elseif vr>100
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmy,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'pow')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigpow,y, x, theta, vr,int,0);
             end
         end
         
                   
         function [value] = getDataPoint_int(obj, targetpar, y, x, theta, vr)
%              int = 'spline';
             int = 'linear';
%              int = 'nearest';
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
            if vr <3
                vr = 3;
            elseif vr>100
                vr =100;
            end
%             
%             if theta<-1
%                 theta = theta+2;
%             elseif theta >1
%                 theta = theta-2;
%             end
               
             if strcmp(targetpar, 'CLv')
%                  if vr <3
%                      value= 0;
%                  elseif vr>100
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCLv_int,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'CDv')
%                  if vr <3
%                      value= 0;
%                  elseif vr>100
%                      value= 0;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCDv_int,y, x, theta, vr,int,0);
%                  end
%                  if value <0
%                      value =0;
%                  end
             elseif strcmp(targetpar, 'Cmx')
%                  if vr <3
%                     value= 1;
%                  elseif vr>100
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmx_int,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'Cmy')
%                  if vr <3
%                     value= 1;
%                  elseif vr>100
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigCmy_int,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'pow')
                value= interpn(obj.Y_vec,obj.X_vec, obj.phase_vec, obj.vr_vec, obj.bigpow_int,y, x, theta, vr,int,0);
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
            maxP = max(max(max(targetValue)));                  
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
        
%              figure
             surface(vr, theta, targetValue)
             xlabel('Reduced Velocity')
             ylabel('Theta')
             view([65 30])
             title([targetpar ' vs ' 'Vr' ' @ x=' num2str(x) ', y=' num2str(y)])  
         end
         
         function plot2D_vr_int(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint_int(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end
        
%              figure
             surface(vr, theta, targetValue)
             xlabel('Reduced Velocity')
             ylabel('Theta')
             view([65 30])
             title([targetpar ' vs ' 'Vr' ' @ x=' num2str(x) ', y=' num2str(y)])  
         end
         
         function plot2D_f(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end
%              figure
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
         function h = plot1D_y(obj, targetpar, y, x, theta, vr)
             for i =1:length(y)
                targetValue(i) = getDataPoint(obj, targetpar, y(i), x, theta/180, vr);
             end
             figure
             h = plot(y, targetValue,'-o');
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
         
         
         function h = plot1D_f_int(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                targetValue(i) = getDataPoint_int(obj, targetpar, y, x, theta/180, vr(i));
             end
%              figure
             h = plot(1./vr, targetValue,'-sg');
             xlabel('f*')
             
%              title([targetpar ' vs ' 'f*' ' @ x=' num2str(x) ', theta=' num2str(theta) ', y=' num2str(y)])
%              title([ 'theta=' num2str(theta) ])
         end         
         function h = plot1D_vr_int(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                targetValue(i) = getDataPoint_int(obj, targetpar, y, x, theta/180, vr(i));
             end
%              figure
             h = plot(vr, targetValue,'-+r');
             xlabel('vr')
%              title([targetpar ' vs ' 'vr' ' @ x=' num2str(x) ', theta=' num2str(theta) ', y=' num2str(y)])
             title([ 'theta=' num2str(theta) ])
         end
         
         function h =  plot1D_theta_int(obj, targetpar, y, x, theta, vr)
             for i =1:length(theta)
                targetValue(i) = getDataPoint_int(obj, targetpar, y, x, theta(i)/180, vr);
             end
%              figure
              h = plot(theta, targetValue,'-+r');
%              xlabel('theta')
%              title([targetpar ' vs ' 'theta' ' @ x=' num2str(x) ', y=' num2str(y) ', vr=' num2str(vr)])
         end
         function h =  plot1D_y_int(obj, targetpar, y, x, theta, vr)
             for i =1:length(y)
                targetValue(i) = getDataPoint_int(obj, targetpar, y(i), x, theta/180, vr);
             end
%              figure
              h = plot(y, targetValue,'-+r');
%              xlabel('y')
%              title([targetpar ' vs ' 'y' ' @ x=' num2str(x) ', theta=' num2str(theta*180,3) ', vr=' num2str(vr)])
         end
         function  h = plot1D_x_int(obj, targetpar, y, x, theta, vr)
             for i =1:length(x)
                targetValue(i) = getDataPoint_int(obj, targetpar, y, x(i), theta/180, vr);
             end
%              figure
              h = plot(x, targetValue,'-+r');
             xlabel('x')
             title([targetpar ' vs ' 'x' ' @ y=' num2str(y) ', theta=' num2str(theta*180,3) ', vr=' num2str(vr)])
         end 
    end
end