function [rawDatabaseStrucJMD] = loadDatabaseJMD()
    
    vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
    theta_set_ori = {'-pi','-3/4pi','-1/2pi','-1/4pi','0','1/4pi','1/2pi','3/4pi'};
    theta_set_ext = {'360n' '315n' '270n' '225n' '180n' '135n' '90n' '45n' '0' '45' '90' '135' '180' '225' '270' '315' '360'};
    X_set = {'0','0.15','0.3','0.45','0.6','0.75'};
    Y_set = {'0.25','0.5','0.75','1','1.25','1.5'};
    theta_vec = [-180 -135 -90 -45 0 45 90 135 180]./180;
    theta_vec_ext = [-360 + [0:16]*45]./180;
    X_vec = [0:0.15:0.75];
    Y_vec = [0.25:0.25:1.5];
    vr_vec = [4.5:0.5:8];
    compType = computer;
    if strcmp(compType, 'PCWIN') == 1
        dataPath =  'D:\Dropbox\VIV\src\ILCFprediction\databaseJMD\';
    else
        dataPath =  '/Users/haining/VIV/src/ILCFprediction/databaseJMD/';
    end
    expName = 'jason';
    
    bigX        = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));  % Ay,Ax,theta, vr
    bigY        = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigtheta    = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCL1      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCL3      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCD2      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCmy      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCmx      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigvr       = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigpow      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCDv      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCDa      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCLa      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    bigCLv      = NaN(length(Y_set),length(X_set),length(theta_set_ext),length(vr_set));
    
    newX        = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newY        = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newtheta    = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCL3      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCL1      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCD2      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCmy      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCmx      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCDv      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCDa      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCLv      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newCLa      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    newvr       = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    avgPow      = NaN(length(Y_set),length(X_set),length(theta_set_ext));
    
    for p = 1:length(vr_set)
        eval(['load ' dataPath 'vr' char(vr_set(p)) '.mat']);

        [ind1,ind2] = find(isnan(Cmx));
        Cmx(ind1,ind2) = 1;%% just for ploting reason

        for i = 1:(length(theta_set_ori))
          for j = 1:length(Y_vec)
             ind = find(allYad(i,:) == Y_vec(j));
             for k = ind
                 newX(j,1:length(ind),i+4) = allXad(i,ind);% 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
                 newY(j,1:length(ind),i+4) = allYad(i,ind);% 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
                 newtheta(j,1:length(ind),i+4) = alltheta(i,ind)./180;% theta rotating, 8 theta, from -180 to 180 degrees in increments of 45 degrees.
%                          newtheta(j,1:length(ind),i) = alltheta(i,ind);
                 newvr(j,1:length(ind),i+4) = Vr(i,ind);
                 newCmy(j,1:length(ind),i+4) = Cmy(i,ind);
                 newCmx(j,1:length(ind),i+4) = Cmx(i,ind);
                 newCLv(j,1:length(ind),i+4) = CLv(i,ind);
                 newCDv(j,1:length(ind),i+4) = CDv(i,ind);
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
                        avgPow(j,1:length(ind),i+4) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.18.^3*0.6858*0.0381);
                    else
                        avgPow(j,1:length(ind),i+4) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.23.^3*0.6858*0.0381);
                    end
                elseif strcmp(expName, 'haining')
                    avgPow(j,1:length(ind),i+4) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.2.^3*0.64135*0.0381);
                end
             end
          end    
        end

        for i = 1:4
            newX(:,:,i) = newX(:,:,i+8);  % 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
            newY(:,:,i) = newY(:,:,i+8);  % 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
            newtheta(:,:,i) = theta_vec_ext(i)./180;      % theta rotating? 8 theta, from -180 to 180 degrees in increments of 45 degrees.
            newvr(:,:,i) = newvr(:,:,i+8);      % 8 reduced velocity, from 4.5 to 8 in increments of 0.5 
            newCmy(:,:,i) = newCmy(:,:,i+8);
            newCmx(:,:,i) = newCmx(:,:,i+8);
            newCDv(:,:,i) = newCDv(:,:,i+8);
            newCLv(:,:,i) = newCLv(:,:,i+8);
            avgPow(:,:,i) = avgPow(:,:,i+8);
        end
        for i = 13:17
            newX(:,:,i) = newX(:,:,i-8);  % 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
            newY(:,:,i) = newY(:,:,i-8);  % 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
            newtheta(:,:,i) = theta_vec_ext(i)./180;      % theta rotating? 8 theta, from -180 to 180 degrees in increments of 45 degrees.
            newvr(:,:,i) = newvr(:,:,i-8);      % 8 reduced velocity, from 4.5 to 8 in increments of 0.5 
            newCmy(:,:,i) = newCmy(:,:,i-8);
            newCmx(:,:,i) = newCmx(:,:,i-8);
            newCDv(:,:,i) = newCDv(:,:,i-8);
            newCLv(:,:,i) = newCLv(:,:,i-8);
            avgPow(:,:,i) = avgPow(:,:,i-8);
        end

        bigX(:,:,:,p) = newX;
        bigY(:,:,:,p) = newY;
        bigtheta(:,:,:,p) = newtheta;
        bigvr(:,:,:,p) = newvr;
        bigCmy(:,:,:,p) = newCmy;
        bigCmx(:,:,:,p) = newCmx;
        bigCLv(:,:,:,p) = newCLv;
        bigCDv(:,:,:,p) = newCDv;
        bigpow(:,:,:,p) = avgPow;
    end
    
        bigCDv(isnan(bigCDv)) = 0 ;
        
    for i = 1:length(Y_set)
        for j = 1:1:length(vr_set)
            bigCmy(i,1,:,j) = mean(bigCmy(i,1,1:8,j));
            bigCmx(i,1,:,j) = mean(bigCmx(i,1,1:8,j));
            bigCLv(i,1,:,j) = mean(bigCLv(i,1,1:8,j));
            bigCDv(i,1,:,j) = mean(bigCDv(i,1,1:8,j));
            bigpow(i,1,:,j) = mean(bigpow(i,1,1:8,j));
        end
    end
            
%     [bigY_2 , bigX_2, bigtheta_2, bigvr_2] = ndgrid(Y_vec, X_vec, theta_vec_ext,vr_vec);
    
    
    [bigY , bigX, bigtheta, bigvr] = ndgrid(Y_vec, X_vec, theta_vec_ext,vr_vec);

%     [bigCmx_smooth_temp, s,exitflag] = smoothn(bigCmx,'robust');
%     [bigCmy_smooth_temp, s,exitflag] = smoothn(bigCmy,'robust');
%     [bigCDv_smooth_temp, s,exitflag] = smoothn(bigCDv,'robust');
%     [bigCLv_smooth_temp, s,exitflag] = smoothn(bigCLv,'robust');
%     [bigpow_smooth_temp, s,exitflag] = smoothn(bigpow,'robust');
% 
%     bigX = bigX(:,:,5:13,:);
%     bigY = bigY(:,:,5:13,:);
%     bigtheta = bigtheta(:,:,5:13,:);
%     bigvr = bigvr(:,:,5:13,:);
%     bigCmx = bigCmx(:,:,5:13,:);
%     bigCmy = bigCmy(:,:,5:13,:);
%     bigCDv = bigCDv(:,:,5:13,:);
%     bigCLv = bigCLv(:,:,5:13,:);
%     bigpow = bigpow(:,:,5:13,:);
%     bigCmx_smooth = bigCmx_smooth_temp(:,:,5:13,:);
%     bigCmy_smooth = bigCmy_smooth_temp(:,:,5:13,:);
%     bigCDv_smooth = bigCDv_smooth_temp(:,:,5:13,:);
%     bigCLv_smooth = bigCLv_smooth_temp(:,:,5:13,:);
%     bigpow_smooth = bigpow_smooth_temp(:,:,5:13,:);
%     
%     totLength = size(bigX,1)* size(bigX,2)* size(bigX,3)* size(bigX,4);
% 
%     allX_vec = reshape(bigX, totLength, 1);
%     allY_vec = reshape(bigY, totLength, 1);
%     alltheta_vec = reshape(bigtheta, totLength, 1);
%     allvr_vec = reshape(bigvr, totLength, 1);
%     allCmx_vec = reshape(bigCmx, totLength, 1);
%     allCmy_vec = reshape(bigCmy, totLength, 1);
%     allCLv_vec = reshape(bigCLv, totLength, 1);
%     allCDv_vec = reshape(bigCDv, totLength, 1);
%     allpow_vec = reshape(bigpow, totLength, 1);
% 
%     datasetMatrix = [allY_vec,allX_vec,alltheta_vec,allvr_vec,allCmy_vec,allCmx_vec,allCLv_vec,allCDv_vec, allpow_vec];
%     
%     rawDatabaseDS = mat2dataset(datasetMatrix,'VarNames',{'Y','X','Theta','Vr',  'Cmy','Cmx', 'CLv', 'CDv','Pow'});

    rawDatabaseStrucJMD = struct('Y',bigY,'X',bigX,'Theta',bigtheta,'Vr',bigvr,'Cmy',bigCmy,'Cmx', bigCmx,'CLv',bigCLv, 'CDv',bigCDv,'Pow',bigpow);
    
end