classdef ILCFHydroModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

   
    properties
        vr_set    
        theta_set
        X_set
        Y_set
        theta_vec
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
         function obj =  ILCFHydroModel()            
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCLv_smooth,y, x, theta, vr,int,0);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCDv_smooth,y, x, theta, vr,int,0);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCmx_smooth,y, x, theta, vr,int,1);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCmy_smooth,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'pow')
                value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigpow_smooth,y, x, theta, vr,int,0);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCLv,y, x, theta, vr,int,0);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCDv,y, x, theta, vr,int,0);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCmx,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'Cmy')
%                  if vr <3
%                     value= 1;
%                  elseif vr>100
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCmy,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'pow')
                value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigpow,y, x, theta, vr,int,0);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCLv_int,y, x, theta, vr,int,0);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCDv_int,y, x, theta, vr,int,0);
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
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCmx_int,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'Cmy')
%                  if vr <3
%                     value= 1;
%                  elseif vr>100
%                     value= -0.5;
%                  else
                    value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigCmy_int,y, x, theta, vr,int,1);
%                  end
             elseif strcmp(targetpar, 'pow')
                value= interpn(obj.Y_vec,obj.X_vec, obj.theta_vec, obj.vr_vec, obj.bigpow_int,y, x, theta, vr,int,0);
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
         
         function plot2D_theta_vr(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end
%              figure
             surface(vr, theta, targetValue)
             xlabel('Reduced Velocity')
             ylabel('Theta')
%              %              view([65 30])
%              title([targetpar ' vs ' 'Vr' ' @ x=' num2str(x) ', y=' num2str(y)])  
         end
         function plot2D_theta_vr_int(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint_int(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end
%              figure
             surface(vr, theta, targetValue)
             xlabel('Reduced Velocity')
             ylabel('Theta')
             %              view([65 30])
%              title([targetpar ' vs ' 'Vr' ' @ x=' num2str(x) ', y=' num2str(y)])  
         end
         function plot2D_theta_vr_ori(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint_ori(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end        
%              figure
             surface(vr, theta, targetValue)
             xlabel('Reduced Velocity')
             ylabel('Theta')
             %              view([65 30])
%              title([targetpar ' vs ' 'Vr' ' @ x=' num2str(x) ', y=' num2str(y)])  
         end
         
         function plot2D_theta_f(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end
%              figure
             surface(1./vr, theta, targetValue)
             xlabel('Non-Dimensional Frequency')
             ylabel('Theta')
             %              view([65 30])
%              title([targetpar ' vs ' 'f*' ' @ x=' num2str(x) ', y=' num2str(y)])
         end
         function plot2D_theta_f_ori(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                 for j = 1:length(theta)
                    targetValue(j,i) = getDataPoint_ori(obj, targetpar, y, x, theta(j)/180, vr(i));
                 end
             end
%              figure
             surface(1./vr, theta, targetValue)
             xlabel('Non-Dimensional Frequency')
             ylabel('Theta')
             %              view([65 30])
%              title([targetpar ' vs ' 'f*' ' @ x=' num2str(x) ', y=' num2str(y)])
         end
         
         function h = plot2D_theta_f_int(obj, targetpar, y, x, theta, vr)
             for i =1:length(vr)
                targetValue(i) = getDataPoint_int(obj, targetpar, y, x, theta/180, vr(i));
             end
%              figure
             h = plot(1./vr, targetValue,'-o');
             xlabel('f*')
%              title([targetpar ' vs ' 'f*' ' @ x=' num2str(x) ', theta=' num2str(theta) ', y=' num2str(y)])
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
         
         function h = plot1D_theta(obj, targetpar, y, x, theta, vr)
             for i =1:length(theta)
                targetValue(i) = getDataPoint(obj, targetpar, y, x, theta(i)/180, vr);
             end
             figure
             h = plot(theta*180, targetValue,'-o');
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
         function h = plot1D_x(obj, targetpar, y, x, theta, vr)
             for i =1:length(x)
                targetValue(i) = getDataPoint(obj, targetpar, y, x(i), theta/180, vr);
             end
             figure
             h = plot(x, targetValue,'-o');
             xlabel('x')
             title([targetpar ' vs ' 'x' ' @ y=' num2str(y) ', theta=' num2str(theta*180,3) ', vr=' num2str(vr)])
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