function plot_para_2D_single(dataPath, fixedpar_1, fixedindex_1, fixedpar_2, fixedindex_2, targetpar, vr_set, theta_set, X_set, Y_set, all_Y, iterNo)
% clear
% plot_para_2D('theta',[1:2],' vr',[1:4],'CL1',1)

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


for m = 1:length(vr_set)
    eval(['load ' dataPath filesep 'database' filesep 'vr' char(vr_set(m)) 'iterNo' num2str(iterNo) '.mat']);
    
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
    end
    
    bigX(:,:,:,m) = newX;
    bigY(:,:,:,m) = newY;
    bigtheta(:,:,:,m) = newtheta;
    bigvr(:,:,:,m) = newvr;
    
    bigCL1(:,:,:,m) = newCL1;
    bigCL3(:,:,:,m) = newCL3;
    bigCD2(:,:,:,m) = newCD2;
    bigpow(:,:,:,m) = avgPow;
    bigCmy(:,:,:,m) = newCmy;
    bigCmx(:,:,:,m) = newCmx;
    bigCDv(:,:,:,m) = newCDv;
    bigCDa(:,:,:,m) = newCDa;
    bigCLv(:,:,:,m) = newCLv;
    bigCLa(:,:,:,m) = newCLa;
end

switch fixedpar_1
     case 'vr'
        switch fixedpar_2
            case 'theta'
            for i = 1:length(fixedindex_1)
                curInd_1 = fixedindex_1(i);
                
                for j = 1:length(fixedindex_2)
                    curInd_2 = fixedindex_2(j);
                    pX(:,:) = bigX(:,:,curInd_2,curInd_1);
                    pY(:,:) = bigY(:,:,curInd_2,curInd_1);
                    pvr(:,:) = bigvr(:,:,curInd_2,curInd_1);
                    ptheta(:,:) = bigtheta(:,:,curInd_2,curInd_1);

                    eval(['p' targetpar '(:,:) = big' targetpar '(:,:,curInd_2,curInd_1);']);
                    eval(['[C h] = contour( pX, pY, p' targetpar ',[-10:0.5:10]);'])

                    set(h,'ShowText','on')
                    colorbar
                    title(strcat('\theta: ', theta_set{curInd_2}, '    Vr: ', vr_set{curInd_1}, ' iterNo ', num2str(iterNo)))
                end
            end
        end  
end
