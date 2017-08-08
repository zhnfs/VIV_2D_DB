classdef ILCFHydroModelJMD < ILCFHydroModel
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

   
    properties
    end
    
    methods
         function obj =  ILCFHydroModelJMD()

            obj.vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
            obj.theta_set = {'-pi','-3/4pi','-1/2pi','-1/4pi','0','1/4pi','1/2pi','3/4pi'};
            obj.X_set = {'0','0.15','0.3','0.45','0.6','0.75'};
            obj.Y_set = {'0.25','0.5','0.75','1','1.25','1.5'};
            obj.theta_vec = [-180 -135 -90 -45 0 45 90 135 180]./180;
	        theta_vec_temp = (-360 + [0:15]*45)./180;
            obj.X_vec = [0:0.15:0.75];
            obj.Y_vec = [0.25:0.25:1.5];
            obj.vr_vec = [4.5:0.5:8];
            compType = computer;
            if strcmp(compType, 'PCWIN') == 1
                dataPath =  'D:\Dropbox\VIV\src\ILCFprediction\databaseJMD\';
            else
                dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseJMD/';
%                 dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseHZ/';
            end
            expName = 'jason';            
          
            [rawDatabaseStruc] = loadDatabaseJMD();        

            obj.bigX = rawDatabaseStruc.X;
            obj.bigY = rawDatabaseStruc.Y;
            obj.bigtheta = rawDatabaseStruc.Theta;
            obj.bigvr = rawDatabaseStruc.Vr;
            obj.bigCmy = rawDatabaseStruc.Cmy;
            obj.bigCmx = rawDatabaseStruc.Cmx;
            obj.bigCLv = rawDatabaseStruc.CLv;
            obj.bigCDv = rawDatabaseStruc.CDv;
            obj.bigpow = rawDatabaseStruc.Pow;
            
            [bigCmx_smooth_temp, s, exitflag] = smoothn(obj.bigCmx,'robust');
            [bigCmy_smooth_temp, s, exitflag] = smoothn(obj.bigCmy,'robust');
            [bigCDv_smooth_temp, s, exitflag] = smoothn(obj.bigCDv,'robust');
            [bigCLv_smooth_temp, s, exitflag] = smoothn(obj.bigCLv,'robust');
            [bigpow_smooth_temp, s, exitflag] = smoothn(obj.bigpow,'robust');
            
            obj.bigX = obj.bigX(:,:,5:13,:);
            obj.bigY = obj.bigY(:,:,5:13,:);
            obj.bigtheta = obj.bigtheta(:,:,5:13,:);
            obj.bigvr = obj.bigvr(:,:,5:13,:);
            obj.bigCmx = obj.bigCmx(:,:,5:13,:);
            obj.bigCmy = obj.bigCmy(:,:,5:13,:);
            obj.bigCDv = obj.bigCDv(:,:,5:13,:);
            obj.bigCLv = obj.bigCLv(:,:,5:13,:);
            obj.bigpow = obj.bigpow(:,:,5:13,:);
            
            obj.bigCmx_smooth = bigCmx_smooth_temp(:,:,5:13,:);
            obj.bigCmy_smooth = bigCmy_smooth_temp(:,:,5:13,:);
            obj.bigCDv_smooth = bigCDv_smooth_temp(:,:,5:13,:);
            obj.bigCLv_smooth = bigCLv_smooth_temp(:,:,5:13,:);
            obj.bigpow_smooth = bigpow_smooth_temp(:,:,5:13,:);
        end

    end
end