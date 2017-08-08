function plot_para(dataPath, outputPath, fixedpar, fixedindex, targetpar, plotMode, vr_set, theta_set, X_set, Y_set, color_set,all_Y )
% plot_para('vr',[1:6],'CL3',1)
% fixedpar = 'theta';
% fixedindex = [2:4];
% plotMode = 1;

% mode 1 draw all figures in one picture
% mode 2 draw each figures in one picture
% mode 3 draw all iso surfaces in one picture
% mode 4 draw each iso surface in one picture

% bigX = NaN(6,6,9,8);  % Ay,Ax,theta, vr

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
    
   for i = 1:length(theta_set)
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

% isoRange = [-0.3:0.1:0.3];
% isoRange_set = {'-0.3','-0.2','-0.1','0','0.1','0.2','0.3'};

scrsz = get(0,'ScreenSize');


switch plotMode
    %% case 1
%     case 1
% 
%         subplotSize_1 = 2;
%         subplotSize_2 = ceil(length(fixedindex)/subplotSize_1);
% 
%         switch fixedpar
%             case 'vr'
%                 figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
% 
%                 for i = 1:length(fixedindex)
%                     curInd = fixedindex(i);
%                     pX(:,:,:) = bigX(:,:,:,curInd);
%                     pY(:,:,:) = bigY(:,:,:,curInd);
%                     pvr(:,:,:) = bigvr(:,:,:,curInd);
%                     ptheta(:,:,:) = bigtheta(:,:,:,curInd);
% 
%                     eval(['p' targetpar '(:,:,:) = big' targetpar '(:,:,:,curInd);']);
%                     subplot(subplotSize_1,subplotSize_2,i)
%                     for j = 1:length(isoRange)
% 
%                         eval([ 'figh(i,j) = patch(isosurface(pX,pY,ptheta,p' targetpar ',' num2str(isoRange(j)) '));']);
% 
%                         set(figh(i,j),'FaceColor',color_set{j});
%                         set(figh(i,j),'FaceAlpha',0.3)
%                         xlabel('Inline Motion')
%                         ylabel('Crossflow Motion')
%                         zlabel('\Theta')
%                         hold on
%                     end
%                     legend(isoRange_set)
%                     title(vr_set{curInd})
%                     grid on
%                     view(-20,30)
%                     hold off
%                 end
%                 saveas(gcf, [outputPath '3D ' 'vr' vr_set{curInd} ' ' targetpar '.jpg'])
%                 saveas(gcf, [outputPath '3D ' 'vr' vr_set{curInd} ' ' targetpar '.fig'])
% 
%                 
% 
%             case 'theta'
%                 figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
%                 for i = 1:length(fixedindex)
%                     curInd = fixedindex(i);
%                     pX(:,:,:) = bigX(:,:,curInd,:);
%                     pY(:,:,:) = bigY(:,:,curInd,:);
%                     pvr(:,:,:) = bigvr(:,:,curInd,:);
%                     ptheta(:,:,:) = bigtheta(:,:,curInd,:);
% 
%                     eval(['p' targetpar '(:,:,:) = big' targetpar '(:,:,curInd,:);']);
%                     subplot(subplotSize_1,subplotSize_2,i)
%                     for j = 1:length(isoRange)
%                         eval([ 'figh(i,j) = patch(isosurface(pX,pY,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);           
% 
%                         set(figh(i,j),'FaceColor',color_set{j});
%                         set(figh(i,j),'FaceAlpha',0.3)
%                         xlabel('Inline Motion')
%                         ylabel('Crossflow Motion')
%                         zlabel('Reduced Velocity')
%                         hold on
%                     end
%                     legend(isoRange_set)
%                     title(theta_set{curInd})
%                     grid on
%                     view(-20,30)
%                     hold off
% 
%                 end
%                 saveas(gcf, [outputPath '3D ' 'theta' theta_set{curInd} ' ' targetpar '.jpg'])
%                 saveas(gcf, [outputPath '3D ' 'theta' theta_set{curInd} ' ' targetpar '.fig'])
% 
% 
%             case 'X'
%                 figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])   
%                 for i = 1:length(fixedindex)
%                     curInd = fixedindex(i);
%                     pX(:,:,:) = bigX(:,curInd,:,:);
%                     pY(:,:,:) = bigY(:,curInd,:,:);
%                     pvr(:,:,:) = bigvr(:,curInd,:,:);
%                     ptheta(:,:,:) = bigtheta(:,curInd,:,:);
% 
%                     eval(['p' targetpar '(:,:,:) = big' targetpar '(i,:,:,:);']);
%                     subplot(subplotSize_1,subplotSize_2,i)
%                     for j = 1:length(isoRange)
%                         eval([ 'figh(i) = patch(isosurface(ptheta,pY,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);
% 
%                         set(figh(i),'FaceColor',color_set{j});
%                         set(figh(i),'FaceAlpha',0.3)
%                         xlabel('\Theta')
%                         ylabel('Crossflow Motion')
%                         zlabel('Reduced Velocity')
%                         hold on
%                     end
%                     legend(isoRange_set)
%                     title(X_set{curInd})
%                     grid on
%                     view(-20,30)
%                     hold off
% 
%                 end
%                 saveas(gcf, [outputPath '3D ' 'X' X_set{curInd} ' ' targetpar '.jpg'])
%                 saveas(gcf, [outputPath '3D ' 'X' X_set{curInd} ' ' targetpar '.fig'])
%                 
% 
%             case 'Y'
%                 figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
%                 for i = 1:length(fixedindex)
%                     curInd = fixedindex(i);
%                     pX(:,:,:) = bigX(curInd,:,:,:);
%                     pY(:,:,:) = bigY(curInd,:,:,:);
%                     pvr(:,:,:) = bigvr(curInd,:,:,:);
%                     ptheta(:,:,:) = bigtheta(curInd,:,:,:);
% 
%                     eval(['p' targetpar '(:,:,:) = big' targetpar '(curInd,:,:,:);']);
%                     subplot(subplotSize_1,subplotSize_2,i)
%                      for j = 1:length(isoRange)
%                         eval([ 'figh(i) = patch(isosurface(ptheta,pX,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);  
% 
%                         set(figh(i),'FaceColor',color_set{j});
%                         set(figh(i),'FaceAlpha',0.3)
%                         xlabel('\Theta')
%                         ylabel('Inline Motion')
%                         zlabel('Reduced Velocity')
%                         hold on
%                     end
%                     legend(isoRange_set)
%                     title(Y_set{curInd}) 
%                     grid on
%                     view(-20,30)
%                     hold off
% 
%                 end
% 
%                 suptitle(targetpar)
%                 saveas(gcf, [outputPath '3D ' 'Y' Y_set{curInd} ' ' targetpar '.jpg'])
%                 saveas(gcf, [outputPath '3D ' 'Y' Y_set{curInd} ' ' targetpar '.fig'])
% 
%         end

%% case 2
    case 2
        switch fixedpar
            case 'vr'
                
                for i = 1:length(fixedindex)
                    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                    curInd = fixedindex(i);
                    pX(:,:,:) = bigX(:,:,:,curInd);
                    pY(:,:,:) = bigY(:,:,:,curInd);
                    pvr(:,:,:) = bigvr(:,:,:,curInd);
                    ptheta(:,:,:) = bigtheta(:,:,:,curInd);

                    eval(['p' targetpar '(:,:,:) = big' targetpar '(:,:,:,curInd);']);
                    
                    eval(['minP = min(min(min(p' targetpar ')));']);
                    eval(['maxP = max(max(max(p' targetpar ')));']); 
                    
                    isoRange = linspace(0, maxP-minP,5);
                    isoRange = minP + isoRange(2:4);
                    isoRange_set = cell(1,3);
                    for m = 1:3
                         isoRange_set{m} = num2str(isoRange(m),2);
                    end                    
                  
                    for j = 1:length(isoRange)

                        eval([ 'figh(i,j) = patch(isosurface(pX,pY,ptheta,p' targetpar ',' num2str(isoRange(j)) '));']);

                        set(figh(i,j),'FaceColor',color_set{j});
                        set(figh(i,j),'FaceAlpha',0.3)
                        xlabel('Inline Motion')
                        ylabel('Crossflow Motion')
                        zlabel('\Theta')
                        hold on
                    end
                    legend(isoRange_set)
                    title([targetpar ' at ' 'vr ' vr_set{curInd} ])
                    grid on
                    view(-20,30)
                    hold off
                end
                saveas(gcf, [outputPath '3D ' 'vr' vr_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'vr' vr_set{curInd} ' ' targetpar '.fig'])


            case 'theta'

                for i = 1:length(fixedindex)
                    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                    curInd = fixedindex(i);
                    pX(:,:,:) = bigX(:,:,curInd,:);
                    pY(:,:,:) = bigY(:,:,curInd,:);
                    pvr(:,:,:) = bigvr(:,:,curInd,:);
                    ptheta(:,:,:) = bigtheta(:,:,curInd,:);

                    eval(['p' targetpar '(:,:,:) = big' targetpar '(:,:,curInd,:);']);
                    
                    eval(['minP = min(min(min(p' targetpar ')));']);
                    eval(['maxP = max(max(max(p' targetpar ')));']);                    
                    isoRange = linspace(0, maxP-minP,5);
                    isoRange = minP + isoRange(2:4);
                    isoRange_set = cell(1,3);
                    for m = 1:3
                         isoRange_set{m} = num2str(isoRange(m),2);
                    end                      
                    for j = 1:length(isoRange)
                        eval([ 'figh(i,j) = patch(isosurface(pX,pY,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);           

                        set(figh(i,j),'FaceColor',color_set{j});
                        set(figh(i,j),'FaceAlpha',0.3)
                        xlabel('Inline Motion')
                        ylabel('Crossflow Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                    legend(isoRange_set)
                    title(['fixed ' fixedpar  ' v.s. ' targetpar])
                    grid on
                    view(-20,30)
                    hold off

                end
                saveas(gcf, [outputPath '3D ' 'theta' theta_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'theta' theta_set{curInd} ' ' targetpar '.fig'])


            case 'X'

                for i = 1:length(fixedindex)
                    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                    curInd = fixedindex(i);
                    pX(:,:,:) = bigX(:,curInd,:,:);
                    pY(:,:,:) = bigY(:,curInd,:,:);
                    pvr(:,:,:) = bigvr(:,curInd,:,:);
                    ptheta(:,:,:) = bigtheta(:,curInd,:,:);

                    eval(['p' targetpar '(:,:,:) = big' targetpar '(:,curInd,:,:);']);

                    eval(['minP = min(min(min(p' targetpar ')));']);
                    eval(['maxP = max(max(max(p' targetpar ')));']);                    
                    isoRange = linspace(0, maxP-minP,5);
                    isoRange = minP + isoRange(2:4);
                    isoRange_set = cell(1,3);
                    for m = 1:3
                         isoRange_set{m} = num2str(isoRange(m),2);
                    end                      
                    for j = 1:length(isoRange)
                        eval([ 'figh(i) = patch(isosurface(ptheta,pY,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);

                        set(figh(i),'FaceColor',color_set{j});
                        set(figh(i),'FaceAlpha',0.3)
                        xlabel('\Theta')
                        ylabel('Crossflow Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                    legend(isoRange_set)
                    title(['fixed ' fixedpar  ' v.s. ' targetpar])
                    grid on
                    view(-20,30)
                    hold off

                end
                saveas(gcf, [outputPath '3D ' 'X' X_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'X' X_set{curInd} ' ' targetpar '.fig'])


            case 'Y'

                for i = 1:length(fixedindex)
                    curInd = fixedindex(i);
                    pX(:,:,:) = bigX(curInd,:,:,:);
                    pY(:,:,:) = bigY(curInd,:,:,:);
                    pvr(:,:,:) = bigvr(curInd,:,:,:);
                    ptheta(:,:,:) = bigtheta(curInd,:,:,:);

                    eval(['p' targetpar '(:,:,:) = big' targetpar '(curInd,:,:,:);']);

                    eval(['minP = min(min(min(p' targetpar ')));']);
                    eval(['maxP = max(max(max(p' targetpar ')));']);                    
                    isoRange = linspace(0, maxP-minP,5);
                    isoRange = minP + isoRange(2:4);
                    isoRange_set = cell(1,3);
                    for m = 1:3
                         isoRange_set{m} = num2str(isoRange(m),2);
                    end                      
                     for j = 1:length(isoRange)
                        eval([ 'figh(i) = patch(isosurface(ptheta,pX,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);  

                        set(figh(i),'FaceColor',color_set{j});
                        set(figh(i),'FaceAlpha',0.3)
                        xlabel('\Theta')
                        ylabel('Inline Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                    legend(isoRange_set)
                    title(['fixed ' fixedpar  ' v.s. ' targetpar])
                    grid on
                    view(-20,30)
                    hold off

                end
                saveas(gcf, [outputPath '3D ' 'Y' Y_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'Y' Y_set{curInd} ' ' targetpar '.fig'])


        end
        
    case 3
        subplotSize_1 = 2;
        subplotSize_2 = ceil(length(isoRange)/subplotSize_1);

        switch fixedpar
            case 'vr'
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                for j = 1:length(isoRange)
                    subplot(subplotSize_1,subplotSize_2,j)
                    for i = 1:length(fixedindex)
                        curInd = fixedindex(i);
                        pX(:,:,:) = bigX(:,:,:,curInd);
                        pY(:,:,:) = bigY(:,:,:,curInd);
                        pvr(:,:,:) = bigvr(:,:,:,curInd);
                        ptheta(:,:,:) = bigtheta(:,:,:,curInd);
                        
                        eval(['p' targetpar '(:,:,:) = big' targetpar '(:,:,:,curInd);']);
                        
                        eval(['minP = min(min(min(p' targetpar ')));']);
                        eval(['maxP = max(max(max(p' targetpar ')));']);                    
                        isoRange = linspace(0, maxP-minP,5);
                        isoRange = minP + isoRange(2:4);
                        isoRange_set = cell(1,3);
                        for m = 1:3
                             isoRange_set{m} = num2str(isoRange(m),2);
                        end                     
                        eval([ 'figh(i,j) = patch(isosurface(pX,pY,ptheta,p' targetpar ',' num2str(isoRange(j)) '));']);

                        set(figh(i,j),'FaceColor',color_set{i});
                        set(figh(i,j),'FaceAlpha',0.3)
                        xlabel('Inline Motion')
                        ylabel('Crossflow Motion')
                        zlabel('\Theta')
                        hold on
                    end
                    legend(vr_set)
                    title(num2str(isoRange(j)))
                    grid on
                    view(-20,30)
                    hold off
                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '3D ' 'vr' vr_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'vr' vr_set{curInd} ' ' targetpar '.fig'])

            case 'theta'
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                for j = 1:length(isoRange)
                    subplot(subplotSize_1,subplotSize_2,j)
                    for i = 1:length(fixedindex)
                        curInd = fixedindex(i);
                        pX(:,:,:) = bigX(:,:,curInd,:);
                        pY(:,:,:) = bigY(:,:,curInd,:);
                        pvr(:,:,:) = bigvr(:,:,curInd,:);
                        ptheta(:,:,:) = bigtheta(:,:,curInd,:);

                        eval(['p' targetpar '(:,:,:) = big' targetpar '(:,:,curInd,:);']);

                        eval(['minP = min(min(min(p' targetpar ')));']);
                        eval(['maxP = max(max(max(p' targetpar ')));']);                    
                        isoRange = linspace(0, maxP-minP,5);
                        isoRange = minP + isoRange(2:4);
                        isoRange_set = cell(1,3);
                        for m = 1:3
                             isoRange_set{m} = num2str(isoRange(m),2);
                        end                  
                        eval([ 'figh(i,j) = patch(isosurface(pX,pY,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);           

                        set(figh(i,j),'FaceColor',color_set{i});
                        set(figh(i,j),'FaceAlpha',0.3)
                        xlabel('Inline Motion')
                        ylabel('Crossflow Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                    legend(theta_set)
                    title(num2str(isoRange(j)))
                    grid on
                    view(-20,30)
                    hold off

                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '3D ' 'theta' theta_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'theta' theta_set{curInd} ' ' targetpar '.fig'])


            case 'X'
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)]) 
                for j = 1:length(isoRange)
                        subplot(subplotSize_1,subplotSize_2,j)
                        for i = 1:length(fixedindex)
                        curInd = fixedindex(i);
                        pX(:,:,:) = bigX(:,curInd,:,:);
                        pY(:,:,:) = bigY(:,curInd,:,:);
                        pvr(:,:,:) = bigvr(:,curInd,:,:);
                        ptheta(:,:,:) = bigtheta(:,curInd,:,:);

                        eval(['p' targetpar '(:,:,:) = big' targetpar '(i,:,:,:);']);
                        
                        eval(['minP = min(min(min(p' targetpar ')));']);
                        eval(['maxP = max(max(max(p' targetpar ')));']);                    
                        isoRange = linspace(0, maxP-minP,5);
                        isoRange = minP + isoRange(2:4);
                        isoRange_set = cell(1,3);
                        for m = 1:3
                             isoRange_set{m} = num2str(isoRange(m),2);
                        end  
                    
                        eval([ 'figh(i) = patch(isosurface(ptheta,pY,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);

                        set(figh(i),'FaceColor',color_set{i});
                        set(figh(i),'FaceAlpha',0.3)                        
                        xlabel('\Theta')
                        ylabel('Crossflow Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                    legend(X_set)
                    title(num2str(isoRange(j)))
                    grid on
                    view(-20,30)
                    hold off

                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '3D ' 'X' X_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'X' X_set{curInd} ' ' targetpar '.fig'])

            case 'Y'
                figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
                for j = 1:length(isoRange)
                        subplot(subplotSize_1,subplotSize_2,j)
                        for i = 1:length(fixedindex)
                        curInd = fixedindex(i);
                        pX(:,:,:) = bigX(curInd,:,:,:);
                        pY(:,:,:) = bigY(curInd,:,:,:);
                        pvr(:,:,:) = bigvr(curInd,:,:,:);
                        ptheta(:,:,:) = bigtheta(curInd,:,:,:);

                        eval(['p' targetpar '(:,:,:) = big' targetpar '(curInd,:,:,:);']);
                        
                        eval(['minP = min(min(min(p' targetpar ')));']);
                        eval(['maxP = max(max(max(p' targetpar ')));']);                    
                        isoRange = linspace(0, maxP-minP,5);
                        isoRange = minP + isoRange(2:4);
                        isoRange_set = cell(1,3);
                        for m = 1:3
                             isoRange_set{m} = num2str(isoRange(m),2);
                        end                    

                        eval([ 'figh(i) = patch(isosurface(ptheta,pX,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);  

                        set(figh(i),'FaceColor',color_set{i});
                        set(figh(i),'FaceAlpha',0.3)
                        xlabel('\Theta')
                        ylabel('Inline Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                    legend(Y_set)
                    title(num2str(isoRange(j)))
                    grid on
                    view(-20,30)
                    hold off
                end
                suptitle(targetpar)
                saveas(gcf, [outputPath '3D ' 'Y' Y_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'Y' Y_set{curInd} ' ' targetpar '.fig'])


        end
        
        case 4
        switch fixedpar
            case 'vr'
                for j = 1:length(isoRange)
                    figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])

                   for i = 1:length(fixedindex)
              
                        curInd = fixedindex(i);
                        pX(:,:,:) = bigX(:,:,:,curInd);
                        pY(:,:,:) = bigY(:,:,:,curInd);
                        pvr(:,:,:) = bigvr(:,:,:,curInd);
                        ptheta(:,:,:) = bigtheta(:,:,:,curInd);
                        
                        eval(['p' targetpar '(:,:,:) = big' targetpar '(:,:,:,curInd);']);
                        
                        eval(['minP = min(min(min(p' targetpar ')));']);
                        eval(['maxP = max(max(max(p' targetpar ')));']);                    
                        isoRange = linspace(0, maxP-minP,5);
                        isoRange = minP + isoRange(2:4);
                        isoRange_set = cell(1,3);
                        for m = 1:3
                             isoRange_set{m} = num2str(isoRange(m),2);
                        end  
                    
                        eval([ 'figh(i,j) = patch(isosurface(pX,pY,ptheta,p' targetpar ',' num2str(isoRange(j)) '));']);

                        set(figh(i,j),'FaceColor',color_set{i});
                        set(figh(i,j),'FaceAlpha',0.3)
                        xlabel('Inline Motion')
                        ylabel('Crossflow Motion')
                        zlabel('\Theta')
                        hold on
                   end
                        legend(vr_set)
                        title(num2str(isoRange(j)))
                        grid on
                        view(-20,30)
                        hold off
                end
                saveas(gcf, [outputPath '3D ' 'vr' vr_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'vr' vr_set{curInd} ' ' targetpar '.fig'])
                
            case 'theta'
                for j = 1:length(isoRange)
                      figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
          
                    for i = 1:length(fixedindex)
                      
                        curInd = fixedindex(i);
                        pX(:,:,:) = bigX(:,:,curInd,:);
                        pY(:,:,:) = bigY(:,:,curInd,:);
                        pvr(:,:,:) = bigvr(:,:,curInd,:);
                        ptheta(:,:,:) = bigtheta(:,:,curInd,:);

                        eval(['p' targetpar '(:,:,:) = big' targetpar '(:,:,curInd,:);']);
       
                        eval(['minP = min(min(min(p' targetpar ')));']);
                        eval(['maxP = max(max(max(p' targetpar ')));']);                    
                        isoRange = linspace(0, maxP-minP,5);
                        isoRange = minP + isoRange(2:4);
                        isoRange_set = cell(1,3);
                        for m = 1:3
                             isoRange_set{m} = num2str(isoRange(m),2);
                        end  
                    
                        eval([ 'figh(i,j) = patch(isosurface(pX,pY,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);           

                        set(figh(i,j),'FaceColor',color_set{i});
                        set(figh(i,j),'FaceAlpha',0.3)
                        xlabel('Inline Motion')
                        ylabel('Crossflow Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                        legend(theta_set)
                        title(num2str(isoRange(j)))
                        grid on
                        view(-20,30)
                        hold off

                end
                saveas(gcf, [outputPath '3D ' 'theta' theta_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'theta' theta_set{curInd} ' ' targetpar '.fig'])
              

            case 'X'
                for j = 1:length(isoRange)
                      figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
          
                    for i = 1:length(fixedindex)
                     
                        curInd = fixedindex(i);
                        pX(:,:,:) = bigX(:,curInd,:,:);
                        pY(:,:,:) = bigY(:,curInd,:,:);
                        pvr(:,:,:) = bigvr(:,curInd,:,:);
                        ptheta(:,:,:) = bigtheta(:,curInd,:,:);

                        eval(['p' targetpar '(:,:,:) = big' targetpar '(i,:,:,:);']);
                        
                        eval(['minP = min(min(min(p' targetpar ')));']);
                        eval(['maxP = max(max(max(p' targetpar ')));']);                    
                        isoRange = linspace(0, maxP-minP,5);
                        isoRange = minP + isoRange(2:4);
                        isoRange_set = cell(1,3);
                        for m = 1:3
                             isoRange_set{m} = num2str(isoRange(m),2);
                        end  
                    
                        eval([ 'figh(i) = patch(isosurface(ptheta,pY,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);

                        set(figh(i),'FaceColor',color_set{i});
                        set(figh(i),'FaceAlpha',0.3)
                        xlabel('\Theta')
                        ylabel('Crossflow Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                        legend(X_set)
                        title(num2str(isoRange(j)))
                        grid on
                        view(-20,30)
                        hold off
                end
                saveas(gcf, [outputPath '3D ' 'X' X_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'X' X_set{curInd} ' ' targetpar '.fig'])

            case 'Y'
                for j = 1:length(isoRange)
                      figure('Position',[1 scrsz(4) scrsz(3) scrsz(4)])
          
                    for i = 1:length(fixedindex)
                       
                        curInd = fixedindex(i);
                        pX(:,:,:) = bigX(curInd,:,:,:);
                        pY(:,:,:) = bigY(curInd,:,:,:);
                        pvr(:,:,:) = bigvr(curInd,:,:,:);
                        ptheta(:,:,:) = bigtheta(curInd,:,:,:);

                        eval(['p' targetpar '(:,:,:) = big' targetpar '(curInd,:,:,:);']);
                        
                        eval(['minP = min(min(min(p' targetpar ')));']);
                        eval(['maxP = max(max(max(p' targetpar ')));']);                    
                        isoRange = linspace(0, maxP-minP,5);
                        isoRange = minP + isoRange(2:4);
                        isoRange_set = cell(1,3);
                        for m = 1:3
                             isoRange_set{m} = num2str(isoRange(m),2);
                        end  
                    
                        eval([ 'figh(i) = patch(isosurface(ptheta,pX,pvr,p' targetpar ',' num2str(isoRange(j)) '));']);  

                        set(figh(i),'FaceColor',color_set{i});
                        set(figh(i),'FaceAlpha',0.3)
                        xlabel('\Theta')
                        ylabel('Inline Motion')
                        zlabel('Reduced Velocity')
                        hold on
                    end
                        legend(Y_set)
                        title(num2str(isoRange(j)))
                        grid on
                        view(-20,30)
                        hold off
                   
                end
                saveas(gcf, [outputPath '3D ' 'Y' Y_set{curInd} ' ' targetpar '.jpg'])
                saveas(gcf, [outputPath '3D ' 'Y' Y_set{curInd} ' ' targetpar '.fig'])

            end

end

close all
