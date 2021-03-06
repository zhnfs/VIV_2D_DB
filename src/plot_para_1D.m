function plot_para_1D(dataPath, outputPath, fixedpar_1, fixedindex_1, fixedpar_2, fixedindex_2,fixedpar_3, fixedindex_3, targetpar,vr_set, theta_set, X_set, Y_set, all_Y  )
% clear
% plot_para_1D('theta',[1:2],' vr',[1:4],'CL1',1)

% vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
% theta_set = {'-pi','-3/4pi','-1/2pi','-1/4pi','0','1/4pi','1/2pi','3/4pi','pi' };
% X_set = {'0','0.15','0.3','0.45','0.6','0.75'};
% Y_set = {'0.25','0.5','0.75','1','1.25','1.5'};
% color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
% 
% targetpar_set = {'CL1','CL3','CD2','Cmy', 'Cmx','pow','CDv'};

bigX        = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));  % Ay,Ax,theta, vr
bigY        = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigtheta    = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCL1      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCL3      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCD2      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCmy      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCmx      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigvr       = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigpow      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCDv      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCDa      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCLa      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));
bigCLv      = NaN(length(Y_set),length(X_set),length(theta_set),length(vr_set));

newX = NaN(length(Y_set),length(X_set),length(theta_set));
newY = NaN(length(Y_set),length(X_set),length(theta_set));
newtheta = NaN(length(Y_set),length(X_set),length(theta_set));
newCL3 = NaN(length(Y_set),length(X_set),length(theta_set));
newCL1 = NaN(length(Y_set),length(X_set),length(theta_set));
newCD2 = NaN(length(Y_set),length(X_set),length(theta_set));
newCmy = NaN(length(Y_set),length(X_set),length(theta_set));
newCmx = NaN(length(Y_set),length(X_set),length(theta_set));
newCDv = NaN(length(Y_set),length(X_set),length(theta_set));
newCDa = NaN(length(Y_set),length(X_set),length(theta_set));
newCLv = NaN(length(Y_set),length(X_set),length(theta_set));
newCLa = NaN(length(Y_set),length(X_set),length(theta_set));
newvr     = NaN(length(Y_set),length(X_set),length(theta_set));
avgPow = NaN(length(Y_set),length(X_set),length(theta_set));


for p = 1:length(vr_set)
    eval(['load ' dataPath 'vr' char(vr_set(p)) '.mat']);
    
    [ind1,ind2] = find(isnan(Cmx));
    Cmx(ind1,ind2) = 1;%% just for ploting reason
%     [ind1,ind2] = find(isnan(CDv));
%     CDv(ind1,ind2) = 0;%% just for ploting reason
    
   for i = 1:(length(theta_set)-1)
       for j = 1:length(all_Y)
            ind = find(allYad(i,:) == all_Y(j));

            for k = ind
                newX(j,1:length(ind),i) = allXad(i,ind);
                newY(j,1:length(ind),i) = allYad(i,ind);
                newtheta(j,1:length(ind),i) = alltheta(i,ind);
                newCL3(j,1:length(ind),i) = CL3(i,ind);
                newCL1(j,1:length(ind),i) = CL1(i,ind);
                newCD2(j,1:length(ind),i) = CD2(i,ind);
                newCmy(j,1:length(ind),i) = Cmy(i,ind);
                newCmx(j,1:length(ind),i) = Cmx(i,ind);
                newCDv(j,1:length(ind),i) = CDv(i,ind);
                newCDa(j,1:length(ind),i) = CDa(i,ind);
                newCLv(j,1:length(ind),i) = CLv(i,ind);
                newCLa(j,1:length(ind),i) = CLa(i,ind);
                newvr(j,1:length(ind),i) = Vr(i,ind);
                avgPow(j,1:length(ind),i) = (Plift(i,ind)+Pdrag(i,ind))./(0.5*1000*0.2.^3*0.64135*0.0381);
            end
       end    
            
%         k = 1;
%         for j = 1:length(all_Y)
%                 
%             incre = 1;
%             while allYad(i,k+incre) == allYad(i,k)
%                 incre = incre +1;
%                 if k+incre > size(allYad,2)
%                     break
%                 end
%             end
%             incre = incre -1;
% 
%             newX(j,1:1+incre,i) = allXad(i,k:k+incre);
%             newY(j,1:1+incre,i) = allYad(i,k:k+incre);
%             newtheta(j,1:1+incre,i) = alltheta(i,k:k+incre);
%             newCL3(j,1:1+incre,i) = CL3(i,k:k+incre);
%             newCL1(j,1:1+incre,i) = CL1(i,k:k+incre);
%             newCD2(j,1:1+incre,i) = CD2(i,k:k+incre);
%             newCmy(j,1:1+incre,i) = Cmy(i,k:k+incre);
%             newCmx(j,1:1+incre,i) = Cmx(i,k:k+incre);
%             newCDv(j,1:1+incre,i) = CDv(i,k:k+incre);
%             newCDa(j,1:1+incre,i) = CDa(i,k:k+incre);
%             newCLv(j,1:1+incre,i) = CLv(i,k:k+incre);
%             newCLa(j,1:1+incre,i) = CLa(i,k:k+incre);
%             newvr(j,1:1+incre,i) = Vr(i,k:k+incre);
% 
%             if strcmp(expName, 'jason')
%     %%             special for jason's experiment
%                 if p == 1 | p == 2 | p == 3 % for vr less than 6
%                     % The Reynolds number is held constant at a value of 8760 for reduced velocities greater than 
%                     % or equal to 6 and the Reynolds number is 6860 for reduced velocities less than 6. The Reynolds
%                     % number cannot be held constant for all motions due to frequency constraints on the experimental apparatus.
% 
%                     % The diameter of the cylinder was 0.0381 m (1.5 inches) in order to maximize Reynolds number 
%                     % and minimize the forcing frequency for forced motions. The span was 0.6858 m (27 inches).
%                     % current vel :0.18m/s =>Re=6858
%                     % current vel :0.23m/s =>Re=8763
%                     avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.18.^3*0.6858*0.0381);
%                 else
%                     avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
%                 end
%             elseif strcmp(expName, 'haining')
%                 avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.2.^3*0.64135*0.0381);
%             end
%             k = k+incre+1;
%         end
   end
    
    newX(:,:,9) = newX(:,:,1);  % 6 in-line amplitudes, from 0 to 0.75 in increments of 0.15 
    newY(:,:,9) = newY(:,:,1);  % 6 transverse amplitudes, from 0.25 to 1.5 in increments of 0.25 
    newtheta(:,:,9) = 180;      % theta rotating? 8 Phase, from -180 to 180 degrees in increments of 45 degrees.
    newCL3(:,:,9) = newCL3(:,:,1);
    newCL1(:,:,9) = newCL1(:,:,1);
    newCD2(:,:,9) = newCD2(:,:,1);
    newCmy(:,:,9) = newCmy(:,:,1);
    newCmx(:,:,9) = newCmx(:,:,1);
    newCDv(:,:,9) = newCDv(:,:,1);
    newCDa(:,:,9) = newCDa(:,:,1);
    newCLv(:,:,9) = newCLv(:,:,1);
    newCLa(:,:,9) = newCLa(:,:,1);
    newvr(:,:,9) = newvr(:,:,1);      % 8 reduced velocity, from 4.5 to 8 in increments of 0.5 
    avgPow(:,:,9) = avgPow(:,:,1);
    
    bigX(:,:,:,p) = newX;
    bigY(:,:,:,p) = newY;
    bigtheta(:,:,:,p) = newtheta;
    bigvr(:,:,:,p) = newvr;
    bigCL1(:,:,:,p) = newCL1;
    bigCL3(:,:,:,p) = newCL3;
    bigCD2(:,:,:,p) = newCD2;
    bigpow(:,:,:,p) = avgPow;
    bigCmy(:,:,:,p) = newCmy;
    bigCmx(:,:,:,p) = newCmx;
    bigCDv(:,:,:,p) = newCDv;
    bigCDa(:,:,:,p) = newCDa;
    bigCLv(:,:,:,p) = newCLv;
    bigCLa(:,:,:,p) = newCLa;
end




% %% interpolation
% ifac = 2;
% int = 'spline';
% bigXi = interpn(bigX,ifac,int);
% bigYi = interpn(bigY,ifac,int);
% bigthetai = interpn(bigtheta,ifac,int);
% bigCL1i = interpn(bigCL1,ifac,int);
% bigCL3i = interpn(bigCL3,ifac,int);
% bigCD2i = interpn(bigCD2,ifac,int);
% bigvri = interpn(bigvr,ifac,int);
% bigpowi = interpn(bigpow,ifac,int);
% bigCmyi = interpn(bigCmy,ifac,int);
% bigCmxi = interpn(bigCmx,ifac,int);

% display(['Reduced frequency: ' vr_set{1:end}])
% display(['Theta: ' theta_set{1:end}])
% display(['Inline displacement: ' X_set{1:end}])
% display(['Crossflow displcament: ' Y_set{1:end}])
% 
% display('choose one of the parameter at a fixed value')

scrsz = get(0,'ScreenSize');

%% each fixedindex_1 corresponding to one figure which may have subplot of different fixedindex_2
if length(fixedindex_2)==1
    subplotSize_1 = 1;
    subplotSize_2 = 1;
else
    subplotSize_2 = 2;
    subplotSize_1 = ceil(length(fixedindex_2)/2);
end

switch fixedpar_1
     case 'vr'
        switch fixedpar_2
            case 'theta'
               switch fixedpar_3
                    case 'X'
                    pX(:,:) = bigX(1,fixedindex_3, fixedindex_2,fixedindex_1);
                    pY(:,:) = reshape(bigY(:,fixedindex_3, fixedindex_2,fixedindex_1),length(Y_set),1);
                    pvr(:,:) = bigvr(1,fixedindex_3, fixedindex_2,fixedindex_1);
                    ptheta(:,:) = bigtheta(1,fixedindex_3, fixedindex_2,fixedindex_1);

                    eval(['p' targetpar '(:,:) = big' targetpar '(:,:,curInd_2,curInd_1);']);
                    eval(['[C h] = contour( pX, pY, p' targetpar ',[-10:0.5:10]);'])

                    set(h,'ShowText','on')
                    if n == subplotSize_1 * subplotSize_2 || n == subplotSize_1 * subplotSize_2-1
                        xlabel('Inline Motion')
                    end
                    ylabel('Crossflow Motion')

                    colorbar
                    title(strcat('\theta: ', theta_set{curInd_2}, '    Vr: ', vr_set{curInd_1}))
                    n = n+1;

                    suptitle(targetpar)
                    saveas(gcf, [outputPath '1D ' ' vr' vr_set{curInd_1} 'theta' theta_set{fixedindex_2(1)} ' to ' theta_set{fixedindex_2(end)} targetpar '.jpg'])
                    saveas(gcf, [outputPath '1D ' ' vr' vr_set{curInd_1} 'theta' theta_set{fixedindex_2(1)} ' to ' theta_set{fixedindex_2(end)} targetpar '.fig'])
            end
        end
  

    case 'theta'
        switch fixedpar_2
            case ' X'
            for i = 1:length(fixedindex_1)
                curInd_1 = fixedindex_1(i);
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                n =1;
                for j = 1:length(fixedindex_2)
                    curInd_2 = fixedindex_2(j);
                    subplot(subplotSize_1,subplotSize_2,n)

                    pX(:,:) = bigX(:,curInd_2,curInd_1,:);
                    pY(:,:) = bigY(:,curInd_2,curInd_1,:);
                    pvr(:,:) = bigvr(:,curInd_2,curInd_1,:);
                    ptheta(:,:) = bigtheta(:,curInd_2,curInd_1,:);

                    eval(['p' targetpar '(:,:) = big' targetpar '(:,curInd_2,curInd_1,:);']);
                    eval(['[C h] = contour( pvr,pY, p' targetpar ');'])

                    set(h,'ShowText','on')
%                             set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
%                             colormap cool
                    if n == subplotSize_1 * subplotSize_2 || n == subplotSize_1 * subplotSize_2-1
                    xlabel('Reduced Velocity')
                    end
                    ylabel('Crossflow Motion')

%                             colorbar
                    title(strcat('\theta: ', theta_set{curInd_1}, '    X: ', X_set{curInd_2}))
                    n = n+1;
                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '1D ' 'theta' theta_set{curInd_1} ' X' X_set{fixedindex_2(1)} ' to ' theta_set{fixedindex_2(end)} targetpar '.jpg'])
%                         saveas(gcf, [outputPath '1D ' 'theta' theta_set{curInd_1} ' X' X_set{fixedindex_2(1)} ' to ' theta_set{fixedindex_2(end)} targetpar '.fig'])
            end

            case ' vr'
            for i = 1:length(fixedindex_1)
                curInd_1 = fixedindex_1(i);
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                n =1;
                for j = 1:length(fixedindex_2)
                    curInd_2 = fixedindex_2(j);
                    subplot(subplotSize_1,subplotSize_2,n)

                    pX(:,:) = bigX(:,:,curInd_1,curInd_2);
                    pY(:,:) = bigY(:,:,curInd_1,curInd_2);
                    pvr(:,:) = bigvr(:,:,curInd_1,curInd_2);
                    ptheta(:,:) = bigtheta(:,:,curInd_1,curInd_2);

                    eval(['p' targetpar '(:,:) = big' targetpar '(:,:,curInd_1,curInd_2);']);
                    eval(['[C h] = contour( pX, pY, p' targetpar ');'])

                    set(h,'ShowText','on')
                    if n == subplotSize_1 * subplotSize_2 || n == subplotSize_1 * subplotSize_2-1
                        xlabel('Inline Motion')
                    end
                    ylabel('Crossflow Motion')

                    colorbar
                    title(strcat('\theta: ', theta_set{curInd_1}, '    Vr: ', vr_set{curInd_2}))
                    n = n+1;
                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '1D ' 'theta' theta_set{curInd_1} ' vr' vr_set{fixedindex_2(1)} ' to ' vr_set{fixedindex_2(end)} targetpar '.jpg'])
%                         saveas(gcf, [outputPath '1D ' 'theta' theta_set{curInd_1} ' vr' vr_set{fixedindex_2(1)} ' to ' vr_set{fixedindex_2(end)} targetpar '.fig'])
            end


        end

    case 'X'
        switch fixedpar_2
            case 'theta'
            for i = 1:length(fixedindex_1)
                curInd_1 = fixedindex_1(i);
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                n =1;
                for j = 1:length(fixedindex_2)
                    curInd_2 = fixedindex_2(j);
                    subplot(subplotSize_1,subplotSize_2,n)

                    pX(:,:) = bigX(:,curInd_1,curInd_2,:);
                    pY(:,:) = bigY(:,curInd_1,curInd_2,:);
                    pvr(:,:) = bigvr(:,curInd_1,curInd_2,:);
                    ptheta(:,:) = bigtheta(:,curInd_1,curInd_2,:);

                    eval(['p' targetpar '(:,:) = big' targetpar '(:,curInd_1,curInd_2,:);']);
                    eval(['[C h] = contour( pvr,pY, p' targetpar ');'])

                    set(h,'ShowText','on')
                    if n == subplotSize_1 * subplotSize_2 || n == subplotSize_1 * subplotSize_2-1
                        xlabel('Reduced Velocity')
                    end
                    ylabel('Crossflow Motion')

                    colorbar
                    title(strcat('\theta: ', theta_set{curInd_2}, '    X: ', X_set{curInd_1}))
                    n = n+1;
                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '1D ' ' X' X_set{curInd_1} 'theta' theta_set{fixedindex_2(1)} ' to ' theta_set{fixedindex_2(end)} targetpar '.jpg'])
%                         saveas(gcf, [outputPath '1D ' ' X' X_set{curInd_1} 'theta' theta_set{fixedindex_2(1)} ' to ' theta_set{fixedindex_2(end)} targetpar '.fig'])

            end
            case 'Y'
            for i = 1:length(fixedindex_1)
                curInd_1 = fixedindex_1(i);
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                n =1;
                for j = 1:length(fixedindex_2)
                    curInd_2 = fixedindex_2(j);
                    subplot(subplotSize_1,subplotSize_2,n)

                    pX(:,:) = bigX(curInd_2,curInd_1,:,:);
                    pY(:,:) = bigY(curInd_2,curInd_1,:,:);
                    pvr(:,:) = bigvr(curInd_2,curInd_1,:,:);
                    ptheta(:,:) = bigtheta(curInd_2,curInd_1,:,:);

                    eval(['p' targetpar '(:,:) = big' targetpar '(curInd_2,curInd_1,:,:);']);
                    eval(['[C h] = contour( pvr,ptheta, p' targetpar ');'])

                    set(h,'ShowText','on')
                    if n == subplotSize_1 * subplotSize_2 || n == subplotSize_1 * subplotSize_2-1
                        xlabel('Reduced Velocity')
                    end
                    ylabel('Theta')

                    colorbar
                    title(strcat('Y: ', Y_set{curInd_2}, '    X: ', X_set{curInd_1}))
                    n = n+1;
                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '1D ' ' X' X_set{curInd_1} 'Y' Y_set{fixedindex_2(1)} ' to '  Y_set{fixedindex_2(end)} targetpar '.jpg'])
                saveas(gcf, [outputPath '1D ' ' X' X_set{curInd_1} 'Y' Y_set{fixedindex_2(1)} ' to '  Y_set{fixedindex_2(end)} targetpar '.fig'])

                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                n =1;
                for j = 1:length(fixedindex_2)
                    curInd_2 = fixedindex_2(j);
                    subplot(subplotSize_1,subplotSize_2,n)

                    pX(:,:) = bigX(curInd_2,curInd_1,:,:);
                    pY(:,:) = bigY(curInd_2,curInd_1,:,:);
                    pvr(:,:) = bigvr(curInd_2,curInd_1,:,:);
                    ptheta(:,:) = bigtheta(curInd_2,curInd_1,:,:);

                    eval(['p' targetpar '(:,:) = big' targetpar '(curInd_2,curInd_1,:,:);']);
                    eval(['[h] = surface( pvr,ptheta, p' targetpar ');'])

                    if n == subplotSize_1 * subplotSize_2 || n == subplotSize_1 * subplotSize_2-1
                        xlabel('Reduced Velocity')
                    end
                    ylabel('Theta')

                    colorbar
                    title(strcat('Y: ', Y_set{curInd_2}, '    X: ', X_set{curInd_1}))
                    n = n+1;
                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '3D ' ' X' X_set{curInd_1} 'Y' Y_set{fixedindex_2(1)} ' to '  Y_set{fixedindex_2(end)} targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' ' X' X_set{curInd_1} 'Y' Y_set{fixedindex_2(1)} ' to '  Y_set{fixedindex_2(end)} targetpar '.fig'])

            end
        end


end
