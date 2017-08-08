classdef ILCFHydroModelHZ < ILCFHydroModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
   
    properties
    end
    
    methods
         function obj =  ILCFHydroModelHZ()

            obj.vr_set = {'4','4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
%             obj.theta_set = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
%                  '15' '30' '45' '67.5' '90' '112.5' '135' '157.5' '180'};
            obj.theta_set = {'180n' '157.5n' '135n' '120n' '105n' '90n' '75n' '60n'  '45n' '22.5n' '0'...
                '22.5' '45' '67.5' '90' '112.5' '135' '157.5' '180'};
            obj.X_set = {'0.05','0.1','0.15','0.2','0.25','0.3'};
            obj.Y_set = {'0.15','0.25','0.5','0.75', '1'};
%             obj.theta_vec = [-180:22.5:-45, -30:15:30, 45:22.5:180]./180;
% 	        theta_vec_temp = [-360:15:-330, -315:22.5:-45, -30:15:30, 45:22.5:315, 330:15:345];
            obj.theta_vec = [-180:22.5:-135, -120:15:-45, -22.5:22.5:180]./180;
            theta_vec_ext = [-360:22.5:-135, -120:15:-45, -22.5:22.5:240, 255:15:315, 337.5,360]./180;

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
            
            [rawDatabaseStruc] = loadDatabaseHZ();        

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
            
            [bigCmx_smooth_temp, s, exitflag] = smoothn(obj.bigCmx_int,'robust');
            [bigCmy_smooth_temp, s, exitflag] = smoothn(obj.bigCmy_int,'robust');
            [bigCDv_smooth_temp, s, exitflag] = smoothn(obj.bigCDv_int,'robust');
            [bigCLv_smooth_temp, s, exitflag] = smoothn(obj.bigCLv_int,'robust');
            [bigpow_smooth_temp, s, exitflag] = smoothn(obj.bigpow_int,'robust');
            
            obj.bigX = obj.bigX(:,:,9:27,:);
            obj.bigY = obj.bigY(:,:,9:27,:);
            obj.bigtheta = obj.bigtheta(:,:,9:27,:);
            obj.bigvr = obj.bigvr(:,:,9:27,:);
            obj.bigCmx = obj.bigCmx(:,:,9:27,:);
            obj.bigCmy = obj.bigCmy(:,:,9:27,:);
            obj.bigCDv = obj.bigCDv(:,:,9:27,:);
            obj.bigCLv = obj.bigCLv(:,:,9:27,:);
            obj.bigpow = obj.bigpow(:,:,9:27,:);           
            
            obj.bigCmx_int =  obj.bigCmx_int(:,:,9:27,:);
            obj.bigCmy_int =  obj.bigCmy_int(:,:,9:27,:);
            obj.bigCDv_int =  obj.bigCDv_int(:,:,9:27,:);
            obj.bigCLv_int =  obj.bigCLv_int(:,:,9:27,:);
            obj.bigpow_int =  obj.bigpow_int(:,:,9:27,:);
            
            obj.bigCmx_smooth = bigCmx_smooth_temp(:,:,9:27,:);
            obj.bigCmy_smooth = bigCmy_smooth_temp(:,:,9:27,:);
            obj.bigCDv_smooth = bigCDv_smooth_temp(:,:,9:27,:);
            obj.bigCLv_smooth = bigCLv_smooth_temp(:,:,9:27,:);
            obj.bigpow_smooth = bigpow_smooth_temp(:,:,9:27,:);
         end
    end
end