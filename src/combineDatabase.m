function rawDatabaseStruc = combineDatabase(rawDatabaseStrucHZ, rawDatabaseStrucJMD)
    % return -360 to 360 extended database
    % extend vr from 4 to 3, 8 to 100
    
    
% vr_set_comb = {'3' '3.5' '4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8' '9' '10' '15' '20' '40' '100'}; 
vr_set_comb = {'3' '3.5' '4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8' '8.5' '9' '9.5' '10' '20' '100'}; 
X_set_comb = {'0','0.05','0.1','0.15','0.2','0.25','0.3','0.45','0.6','0.75'};
Y_set_comb = {'0.15','0.25','0.5','0.75', '1', '1.25', '1.5'};
% theta_set_comb = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
%      '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
theta_set_comb = {'180n' '157.5n' '135n' '120n' '105n' '90n' '75n' '60n'  '45n' '22.5n' '0'...
     '22.5' '45' '67.5' '90' '112.5' '135' '157.5' '180n'};
theta_set_ext = {'360n' '337.5n' '315n' '292.5n' '270n' '247.5n' '225n' '202.5n' '180n' '157.5n' '135n' '120n' '105n' '90n' '75n' '60n'  '45n' '22.5n' '0'...
        '22.5' '45' '67.5' '90' '112.5' '135' '157.5' '180' '202.5' '225' '240' '255' '270' '285' '300' '315' '337.5' '360'};  

Y_vec_comb = [0.15, 0.25, 0.5, 0.75, 1, 1.25, 1.5 ];
X_vec_comb = [0,0.05,0.1,0.15,0.2,0.25,0.3,0.45,0.6,0.75 ];
theta_vec_comb = [-180:22.5:-135, -120:15:-45, -22.5:22.5:180]./180;
% theta_vec_ext = [-360:15:-330, -315:22.5:-45, -30:15:30, 45:22.5:315, 330:15:345]./180;
theta_vec_ext = [-360:22.5:-135, -120:15:-45, -22.5:22.5:225, 240:15:315, 337.5,360]./180;

vr_vec_comb = [3:0.5:10,20,100 ];

% bigX        = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));  % Ay,Ax,theta, vr
% bigY        = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));
% bigtheta    = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));
% bigvr       = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));
bigCmy      = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));
bigCmx      = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));
bigCLv      = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));
bigCDv      = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));
bigpow      = NaN(length(Y_set_comb),length(X_set_comb),length(theta_set_ext),length(vr_set_comb));

[bigY , bigX, bigtheta, bigvr] = ndgrid(Y_vec_comb, X_vec_comb, theta_vec_ext,vr_vec_comb);

% insert JMD into big matrix

% 6*6*8*8 = 2304  6*6*9*8 = 2592
Y_set_JMD = {'0.25','0.5','0.75','1', '1.25', '1.5'};
X_set_JMD = {'0','0.15','0.3','0.45','0.6','0.75'};
% theta_set_JMD = {'180n' '135n' '90n' '45n' '0' '45' '90' '135'};
theta_set_JMD_ext = {'360n' '315n' '270n' '225n' '180n' '135n' '90n' '45n' '0' '45' '90' '135' '180' '225' '270' '315' '360'};
vr_set_JMD = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'}; 

% 
% AmpYIndex_JMD = [2:7];
% AmpXIndex_JMD = [1,4,7:10];
% % JasonThetaIndex = [1, 3, 5, 7, 10, 13, 15, 17];
% ThetaIndex_JMD = [1 4 6 8 10 12 14 16 19 22 24 26 28 30 32 34];
% VrIndex_JMD = [2:9];

for i = 1:length(Y_set_JMD)
    for j= 1:length(X_set_JMD)
        for m = 1:length(theta_set_JMD_ext)
            for n = 1:length(vr_set_JMD)
               ind_1 = find(abs(Y_vec_comb - rawDatabaseStrucJMD.Y(i,j,m,n))<1e-5);
               ind_2 = find(abs(X_vec_comb - rawDatabaseStrucJMD.X(i,j,m,n))<1e-5);
               ind_3 = find(abs(theta_vec_ext - rawDatabaseStrucJMD.Theta(i,j,m,n))<1e-5);
               [temp, ind_4] = min(abs((vr_vec_comb-rawDatabaseStrucJMD.Vr(i,j,m,n))));              
%                ind_1 = find(Y_vec_comb==rawDatabaseStrucJMD.Y(i,j,m,n));
%                ind_2 = find(X_vec_comb==rawDatabaseStrucJMD.X(i,j,m,n));
%                ind_3 = find(theta_vec_ext==rawDatabaseStrucJMD.Theta(i,j,m,n));
%                [temp, ind_4] = min(abs((vr_vec_comb-rawDatabaseStrucJMD.Vr(i,j,m,n))));
               bigX(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.X(i,j,m,n);
               bigY(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.Y(i,j,m,n);
               bigtheta(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.Theta(i,j,m,n);
               bigvr(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.Vr(i,j,m,n);
               bigCmy(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.Cmy(i,j,m,n);
               bigCmx(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.Cmx(i,j,m,n);
               bigCLv(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.CLv(i,j,m,n);
               bigCDv(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.CDv(i,j,m,n);
               bigpow(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucJMD.Pow(i,j,m,n);               
            end
        end
    end
end

% bigX(AmpYIndex_JMD, AmpXIndex_JMD,ThetaIndex_JMD, VrIndex_JMD)=rawDatabaseStrucJMD.X;




%% insert HZ into big matrix
    
vr_set_HZ = {'4' '4.5' '5' '5.5' '6' '6.5' '7' '7.5' '8'}; 
X_set_HZ = {'0.05','0.1','0.15','0.2','0.25','0.3'};
Y_set_HZ = {'0.15','0.25','0.5','0.75', '1'};
% theta_set_HZ = {'180n' '157.5n' '135n' '112.5n' '90n' '67.5n' '45n' '30n' '15n' '0'...
%      '15' '30' '45' '67.5' '90' '112.5' '135' '157.5'};
% theta_set_HZ =  {'180n' '157.5n' '135n' '120n' '105n' '90n' '75n' '60n'  '45n' '22.5n' '0'...
%      '22.5' '45' '67.5' '90' '112.5' '135' '157.5'};
theta_set_HZ_ext = {'360n' '337.5n' '315n' '292.5n' '270n' '247.5n' '225n' '202.5n' '180n' '157.5n' '135n' '120n' '105n' '90n' '75n' '60n'  '45n' '22.5n' '0'...
        '22.5' '45' '67.5' '90' '112.5' '135' '157.5' '180' '202.5' '225' '240' '255' '270' '285' '300' '315' '337.5' '360'};  

 for i = 1:length(Y_set_HZ)
    for j= 1:length(X_set_HZ)
        for m = 1:length(theta_set_HZ_ext)
            for n = 1:length(vr_set_HZ)
               ind_1 = find(abs(Y_vec_comb - rawDatabaseStrucHZ.Y(i,j,m,n))<1e-5);
               ind_2 = find(abs(X_vec_comb - rawDatabaseStrucHZ.X(i,j,m,n))<1e-5);
               ind_3 = find(abs(theta_vec_ext - rawDatabaseStrucHZ.Theta(i,j,m,n))<1e-5);
               [temp, ind_4] = min(abs((vr_vec_comb-rawDatabaseStrucHZ.Vr(i,j,m,n))));
               if isnan(bigCmy(ind_1,ind_2,ind_3,ind_4))
                   bigX(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.X(i,j,m,n);
                   bigY(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.Y(i,j,m,n);
                   bigtheta(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.Theta(i,j,m,n);
                   bigvr(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.Vr(i,j,m,n);
                   bigCmy(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.Cmy(i,j,m,n);
                   bigCmx(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.Cmx(i,j,m,n);
                   bigCLv(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.CLv(i,j,m,n);
                   bigCDv(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.CDv(i,j,m,n);
                   bigpow(ind_1,ind_2,ind_3,ind_4) = rawDatabaseStrucHZ.Pow(i,j,m,n);
               else
                   bigX(ind_1,ind_2,ind_3,ind_4) = (bigX(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.X(i,j,m,n))/2;
                   bigY(ind_1,ind_2,ind_3,ind_4) = (bigY(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.Y(i,j,m,n))/2;
                   bigtheta(ind_1,ind_2,ind_3,ind_4) = (bigtheta(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.Theta(i,j,m,n))/2;
                   bigvr(ind_1,ind_2,ind_3,ind_4) = (bigvr(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.Vr(i,j,m,n))/2;
                   bigCmy(ind_1,ind_2,ind_3,ind_4) = (bigCmy(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.Cmy(i,j,m,n))/2;
                   bigCmx(ind_1,ind_2,ind_3,ind_4) = (bigCmx(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.Cmx(i,j,m,n))/2;
                   bigCLv(ind_1,ind_2,ind_3,ind_4) = (bigCLv(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.CLv(i,j,m,n))/2;
                   bigCDv(ind_1,ind_2,ind_3,ind_4) = (bigCDv(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.CDv(i,j,m,n))/2;
                   bigpow(ind_1,ind_2,ind_3,ind_4) = (bigpow(ind_1,ind_2,ind_3,ind_4) + rawDatabaseStrucHZ.Pow(i,j,m,n))/2;                   
               end
            end
        end
    end
 end

 %%  force values for low and high Vr 
 for i = 1:length(Y_set_comb)
    for j= 1:length(X_set_comb)
        for m = 1:length(theta_set_ext)
           bigCmy(i,j,m,1:3) = linspace(1,bigCmy(i,j,m,3),3);
           bigCmx(i,j,m,1:3) = linspace(1,bigCmx(i,j,m,3),3);
           bigCLv(i,j,m,1:3) = linspace(0,bigCLv(i,j,m,3),3);
           bigCDv(i,j,m,1:3) = linspace(0,bigCDv(i,j,m,3),3);
           bigpow(i,j,m,1:3) = linspace(0,bigpow(i,j,m,3),3);
           bigCmy(i,j,m,(end-1):end) = -0.5;
           bigCmx(i,j,m,(end-1):end) = -0.5;
           bigCLv(i,j,m,(end-1):end) = 0;
           bigCDv(i,j,m,(end-1):end) = 0;
           bigpow(i,j,m,(end-1):end) = 0;  
           bigCmy(i,j,m,(length(vr_set_comb)-6): (length(vr_set_comb)-2)) = linspace(bigCmy(i,j,m,length(vr_set_comb)-6), -0.5,5);
           bigCmx(i,j,m,(length(vr_set_comb)-6): (length(vr_set_comb)-2)) = linspace(bigCmx(i,j,m,length(vr_set_comb)-6), -0.5,5);
           bigCLv(i,j,m,(length(vr_set_comb)-6): (length(vr_set_comb)-2)) = linspace(bigCLv(i,j,m,length(vr_set_comb)-6), 0,5);
           bigCDv(i,j,m,(length(vr_set_comb)-6): (length(vr_set_comb)-2)) = linspace(bigCDv(i,j,m,length(vr_set_comb)-6), 0,5);
           bigpow(i,j,m,(length(vr_set_comb)-6): (length(vr_set_comb)-2)) = linspace(bigpow(i,j,m,length(vr_set_comb)-6), 0,5);
%            bigCmy(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = linspace(bigCmy(i,j,m,length(vr_set_comb)-6), mean(bigCmy(i,j,:,length(vr_set_comb)-6)),5);
%            bigCmx(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = linspace(bigCmx(i,j,m,length(vr_set_comb)-6), mean(bigCmx(i,j,:,length(vr_set_comb)-6)),5);
%            bigCLv(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = linspace(bigCLv(i,j,m,length(vr_set_comb)-6), mean(bigCLv(i,j,:,length(vr_set_comb)-6)),5);
%            bigCDv(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = linspace(bigCDv(i,j,m,length(vr_set_comb)-6), mean(bigCDv(i,j,:,length(vr_set_comb)-6)),5);
%            bigpow(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = linspace(bigpow(i,j,m,length(vr_set_comb)-6), mean(bigpow(i,j,:,length(vr_set_comb)-6)),5);
        end
    end
 end

%   %%  force values for low and high Vr 
%  for i = 1:length(Y_set_comb)
%     for j= 1:length(X_set_comb)
%         for m = 1:length(theta_set_ext)
%            bigCmy(i,j,m,1:3) = linspace(1,bigCmy(i,j,m,3),3);
%            bigCmx(i,j,m,1:3) = linspace(1,bigCmx(i,j,m,3),3);
%            bigCLv(i,j,m,1:3) = linspace(0,bigCLv(i,j,m,3),3);
%            bigCDv(i,j,m,1:3) = linspace(0,bigCDv(i,j,m,3),3);
%            bigpow(i,j,m,1:3) = linspace(0,bigpow(i,j,m,3),3); 
%             
%            bigCmy(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log10(bigCmy(i,j,m,length(vr_set_comb)-6)-1.1*min(bigCmy(i,j,m,length(vr_set_comb)-6),-0.5)), log10(-0.5+(-1.1*min(bigCmy(i,j,m,length(vr_set_comb)-6),-0.5))),7)+1.1*min(bigCmy(i,j,m,length(vr_set_comb)-6),-0.5);
%            bigCmx(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log10(bigCmx(i,j,m,length(vr_set_comb)-6)-1.1*min(bigCmx(i,j,m,length(vr_set_comb)-6),-0.5)), log10(-0.5+(-1.1*min(bigCmx(i,j,m,length(vr_set_comb)-6),-0.5))),7)+1.1*min(bigCmx(i,j,m,length(vr_set_comb)-6),-0.5);
%            bigCLv(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log10(bigCLv(i,j,m,length(vr_set_comb)-6)-1.1*min(bigCLv(i,j,m,length(vr_set_comb)-6),0)), log10(-0.5+(-1.1*min(bigCLv(i,j,m,length(vr_set_comb)-6),0))),7)+1.1*min(bigCLv(i,j,m,length(vr_set_comb)-6),0);
%            bigCDv(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log10(bigCDv(i,j,m,length(vr_set_comb)-6)-1.1*min(bigCDv(i,j,m,length(vr_set_comb)-6),0)), log10(-0.5+(-1.1*min(bigCDv(i,j,m,length(vr_set_comb)-6),0))),7)+1.1*min(bigCDv(i,j,m,length(vr_set_comb)-6),0);
%            bigpow(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log10(bigpow(i,j,m,length(vr_set_comb)-6)-1.1*min(bigpow(i,j,m,length(vr_set_comb)-6),0)), log10(-0.5+(-1.1*min(bigpow(i,j,m,length(vr_set_comb)-6),0))),7)+1.1*min(bigpow(i,j,m,length(vr_set_comb)-6),0);
%            
% %            bigCmy(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log(bigCmy(i,j,m,length(vr_set_comb)-6)), -0.5,7);
% %            bigCmx(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log(bigCmx(i,j,m,length(vr_set_comb)-6)), -0.5,7);
% %            bigCLv(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log(bigCLv(i,j,m,length(vr_set_comb)-6)), 0,7);
% %            bigCDv(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log(bigCDv(i,j,m,length(vr_set_comb)-6)), 0,7);
% %            bigpow(i,j,m,(length(vr_set_comb)-6): length(vr_set_comb)) = logspace(log(bigpow(i,j,m,length(vr_set_comb)-6)), mean(bigpow(i,j,:,length(vr_set_comb)-6)),7);
%         end
%     end
%  end

%   %% insert Vr 3.5, Vr 8.5
% for i = 1:length(Y_set_comb)
%     for j= 1:length(X_set_comb)
%         for m = 1:length(theta_set_ext)
%             for n = 1:4
%                bigCmy(i,j,m,n) = 1;
%                bigCmx(i,j,m,n) = 1;
%                bigCLv(i,j,m,n) = 0;
%                bigCDv(i,j,m,n) = 0;
%                bigpow(i,j,m,n) = 0;
%             end
%             for n = length(vr_set_comb)-3: length(vr_set_comb)
%                bigCmy(i,j,m,n) = -0.5;
%                bigCmx(i,j,m,n) = -0.5;
%                bigCLv(i,j,m,n) = 0;
%                bigCDv(i,j,m,n) = 0;
%                bigpow(i,j,m,n) = 0;     
%             end
%         end
%     end
% end  
rawDatabaseStruc = struct('Y',bigY,'X',bigX,'Theta',bigtheta,'Vr',bigvr,'Cmy',bigCmy,'Cmx', bigCmx,'CLv',bigCLv, 'CDv',bigCDv,'Pow',bigpow);
