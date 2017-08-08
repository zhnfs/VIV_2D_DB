classdef ILCFHydroModelZDMiT < ILCFHydroModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

   
    properties
    end
    
    methods
         function obj =  ILCFHydroModelZDMiT()

%             obj.vr_set = {'3' '3.5' '4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8' '9' '10' '15' '20' '40' '100'}; 
            obj.vr_set = {'3' '3.5' '4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8' '8.5' '9' '9.5' '10' '20' '100'};  
            obj.X_set = {'0','0.05','0.1','0.15','0.2','0.25','0.3','0.45','0.6','0.75'};
            obj.Y_set = {'0.15','0.25','0.5','0.75', '1', '1.25', '1.5'};

%             obj.theta_set = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
%                  '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
%              
            obj.theta_set = {'180n' '157.5n' '135n' '120n' '105n' '90n' '75n' '60n'  '45n' '22.5n' '0'...
                '22.5' '45' '67.5' '90' '112.5' '135' '157.5' '180'};
 
            obj.Y_vec = [0.15, 0.25, 0.5, 0.75, 1, 1.25, 1.5 ];
            obj.X_vec = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.45,0.6,0.75 ];
%             obj.theta_vec = [-180:22.5:-45, -30:15:30, 45:22.5:180]./180;
            obj.theta_vec =  [-180:22.5:-135, -120:15:-45, -22.5:22.5:180]./180;
            
            obj.vr_vec = [3:0.5:10,20,100];
            theta_vec_ext = [-360:22.5:-135, -120:15:-45, -22.5:22.5:225, 240:15:315, 337.5,360]./180;

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
            
%             smoPar = 0.1;
%             [bigCmx_smooth_temp, s, exitflag] = smoothn(obj.bigCmx_int,smoPar);
%             [bigCmy_smooth_temp, s, exitflag] = smoothn(obj.bigCmy_int,smoPar);
%             [bigCDv_smooth_temp, s, exitflag] = smoothn(obj.bigCDv_int,smoPar);
%             [bigCLv_smooth_temp, s, exitflag] = smoothn(obj.bigCLv_int,smoPar);
%             [bigpow_smooth_temp, s, exitflag] = smoothn(obj.bigpow_int,smoPar);
%             
            
            [bigCmx_smooth_temp, s, exitflag] = smoothn(obj.bigCmx,'robust');
            [bigCmy_smooth_temp, s, exitflag] = smoothn(obj.bigCmy,'robust');
            [bigCDv_smooth_temp, s, exitflag] = smoothn(obj.bigCDv,'robust');
            [bigCLv_smooth_temp, s, exitflag] = smoothn(obj.bigCLv,'robust');
            [bigpow_smooth_temp, s, exitflag] = smoothn(obj.bigpow,'robust');
            
            obj.bigX = obj.bigX(:,:,9:27,:);
            obj.bigY = obj.bigY(:,:,9:27,:);
            obj.bigtheta = obj.bigtheta(:,:,9:27,:);
            obj.bigvr = obj.bigvr(:,:,9:27,:);
            
            obj.bigCmx = obj.bigCmx(:,:,9:27,:);
            obj.bigCmy = obj.bigCmy(:,:,9:27,:);
            obj.bigCDv = obj.bigCDv(:,:,9:27,:);
            obj.bigCLv = obj.bigCLv(:,:,9:27,:);
            obj.bigpow = obj.bigpow(:,:,9:27,:);
            
            obj.bigCmx_int = obj.bigCmx_int(:,:,9:27,:);
            obj.bigCmy_int = obj.bigCmy_int(:,:,9:27,:);
            obj.bigCDv_int = obj.bigCDv_int(:,:,9:27,:);
            obj.bigCLv_int = obj.bigCLv_int(:,:,9:27,:);
            obj.bigpow_int = obj.bigpow_int(:,:,9:27,:);
            
            obj.bigCmx_smooth = bigCmx_smooth_temp(:,:,9:27,:);
            obj.bigCmy_smooth = bigCmy_smooth_temp(:,:,9:27,:);
            obj.bigCDv_smooth = bigCDv_smooth_temp(:,:,9:27,:);
            obj.bigCLv_smooth = bigCLv_smooth_temp(:,:,9:27,:);
            obj.bigpow_smooth = bigpow_smooth_temp(:,:,9:27,:);
            
                  
            % make sure theta has no influence when x = 0
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
                       obj.bigCmy_smooth(i,j,m,(end-1):end) =  -0.5;
                       obj.bigCmx_smooth(i,j,m,(end-1):end) =  -0.5;
                       obj.bigCLv_smooth(i,j,m,(end-1):end) = 0;
                       obj.bigCDv_smooth(i,j,m,(end-1):end) = 0;
                       obj.bigpow_smooth(i,j,m,(end-1):end) = 0;
%                        obj.bigCmy_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigCmy_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigCmy_smooth(i,j,:,length(obj.vr_set)-6)),7);
%                        obj.bigCmx_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigCmx_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigCmx_smooth(i,j,:,length(obj.vr_set)-6)),7);
%                        obj.bigCLv_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigCLv_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigCLv_smooth(i,j,:,length(obj.vr_set)-6)),7);
%                        obj.bigCDv_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigCDv_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigCDv_smooth(i,j,:,length(obj.vr_set)-6)),7);
%                        obj.bigpow_smooth(i,j,m,(length(obj.vr_set)-6): length(obj.vr_set)) = linspace(obj.bigpow_smooth(i,j,m,length(obj.vr_set)-6), mean(obj.bigpow_smooth(i,j,:,length(obj.vr_set)-6)),7);
                       obj.bigCmy_smooth(i,j,m,(length(obj.vr_set)-6): (length(obj.vr_set)-2)) = linspace(obj.bigCmy_smooth(i,j,m,length(obj.vr_set)-6),-0.5,5);
                       obj.bigCmx_smooth(i,j,m,(length(obj.vr_set)-6): (length(obj.vr_set)-2)) = linspace(obj.bigCmx_smooth(i,j,m,length(obj.vr_set)-6),-0.5,5);
                       obj.bigCLv_smooth(i,j,m,(length(obj.vr_set)-6): (length(obj.vr_set)-2)) = linspace(obj.bigCLv_smooth(i,j,m,length(obj.vr_set)-6), 0,5);
                       obj.bigCDv_smooth(i,j,m,(length(obj.vr_set)-6): (length(obj.vr_set)-2)) = linspace(obj.bigCDv_smooth(i,j,m,length(obj.vr_set)-6), 0,5);
                       obj.bigpow_smooth(i,j,m,(length(obj.vr_set)-6): (length(obj.vr_set)-2)) = linspace(obj.bigpow_smooth(i,j,m,length(obj.vr_set)-6), 0,5);
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
    end
end