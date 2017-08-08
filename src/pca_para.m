function pca_para(dataPath, outputDir, vr_set, theta_set, X_set, Y_set, Y_vec  )
% clear
% plot_para_2D('theta',[1:2],' vr',[1:4],'CL1',1)

% vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
% theta_set = {'-pi','-3/4pi','-1/2pi','-1/4pi','0','1/4pi','1/2pi','3/4pi','pi' };
% X_set = {'0','0.15','0.3','0.45','0.6','0.75'};
% Y_set = {'0.25','0.5','0.75','1','1.25','1.5'};
% color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
% 
TargetVar_set = {'Cmy', 'Cmx','CLv','CLa','CDv','CDa','pow'};
para_set={'AmpX', 'AmpY', 'theta', 'Vr'};
col = distinguishable_colors(20);
col = col([13,15,3,8,1,10,14,7,4],:);
shp = ['o';'s';'p';'v';'*';'d';'^';'>';'<';'.'];

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
       for j = 1:length(Y_vec)
            ind = find(allYad(i,:) == Y_vec(j));

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

% Convert 4D matrix into 1D vector

allX_vec = reshape(bigX, size(bigX,1)* size(bigX,2)* size(bigX,3)* size(bigX,4),1);
allY_vec = reshape(bigY, size(bigX,1)* size(bigX,2)* size(bigX,3)* size(bigX,4),1);
alltheta_vec = reshape(bigtheta, size(bigX,1)* size(bigX,2)* size(bigX,3)* size(bigX,4),1);
allvr_vec = reshape(bigvr, size(bigX,1)* size(bigX,2)* size(bigX,3)* size(bigX,4),1);

for i =1:length(TargetVar_set)
    TargetVarName = TargetVar_set{i};
    eval(['target_vec = reshape(big' TargetVarName ', size(bigX,1)* size(bigX,2)* size(bigX,3)* size(bigX,4),1);']);


    A = [target_vec allX_vec allY_vec alltheta_vec allvr_vec];
    AwoNaN = A(find(sum(isnan(A),2)==0),:);

    figure
    plot(target_vec,'b')
    hold on
    plot(allX_vec-2,'r')
    plot(allY_vec+4,'k')
    plot(alltheta_vec./360,'g')
    plot(allvr_vec./2,'m')
    legend(TargetVarName, 'X-2', 'Y+4', 'theta/360', 'vr/2')
    suptitle(TargetVarName)
    saveas(gcf, [outputDir 'Dependency' TargetVarName '.fig'])
    saveas(gcf, [outputDir 'Dependency' TargetVarName '.jpg'])
    close
    
    % % pca analysis
    % [coeff, score, latent, tsquare] = princomp(AwoNaN);

    % multivariable regression
    y = AwoNaN(:,1);
    X = AwoNaN(:,2:5);
    
    corrTargetVar(i,:) = corr(y, X);
    
    % scatter plot
    figure
    for j =1:4
        subplot(2,2,j)
        scatter(X(:,j), y)
        title(para_set{j})
    end
    suptitle(TargetVarName)
    saveas(gcf, [outputDir 'Scatter' TargetVarName '.fig'])
    saveas(gcf, [outputDir 'Scatter' TargetVarName '.jpg'])
    close
    
    % multivariable linear regression
    X = [ones(size(y)) X];
    % stats are the R2 statistic, the F statistic and its p value, and an estimate of the error variance.
    [b(i,:),bint,r,rint,stats(i,:)] = regress(y,X);


    % square reduced velocity to make it more linear
    A_2 = [target_vec allX_vec allY_vec alltheta_vec allvr_vec.^2];
    AwoNaN_2 = A_2(find(sum(isnan(A_2),2)==0),:);
    y_2 = AwoNaN_2(:,1);
    X_2 = AwoNaN_2(:,2:5);
    X_2 = [ones(size(y_2)) X_2 ];
    [b_2(i,:),bint_2,r_2,rint_2,stats_2(i,:)] = regress(y_2,X_2);

%     % multivariable nonlinear regression
%     beta0 = [1 1 1 1 1 1 1];
%     [beta,R,J,CovB,MSE,ErrorModelInfo ]= nlinfit(X, y,@hougen,beta0);

end

figure
for i = 1:length(TargetVar_set)-1
    plot(abs(corrTargetVar(i,:)),'Color',col(i,:),'Marker',shp(i),'MarkerSize',6, 'LineWidth',3)
    hold on
end
xlim([0 5])
legend(TargetVar_set)
set(gca,'xtick',[ 1 2 3 4])
set(gca,'xticklabel', para_set)
title(['Correlations'],'FontSize',12)
saveas(gcf, [outputDir 'Correlations' '.fig'])
saveas(gcf, [outputDir 'Correlations' '.jpg'])
close

for i =1:length(TargetVar_set)-1
    TargetVar_set_Rsquare{i} = [TargetVar_set{i} ', R^2=' num2str(stats(i,1),2)];
end

figure
for i = 1:length(TargetVar_set)-1
    plot(abs(b(i,2:5)),'Color',col(i,:),'Marker',shp(i),'MarkerSize',6, 'LineWidth',3)
    hold on
end
xlim([0 5])
legend(TargetVar_set_Rsquare)
set(gca,'xtick',[ 1 2 3 4])
set(gca,'xticklabel', para_set)
title(['Linear Regression'],'FontSize',12)
saveas(gcf, [outputDir 'LinearRegression' '.fig'])
saveas(gcf, [outputDir 'LinearRegression' '.jpg'])
close

for i =1:length(TargetVar_set)-1
    TargetVar_set_Rsquare_2{i} = [TargetVar_set{i} ', R^2=' num2str(stats_2(i,1),2)];
end
figure
for i = 1:length(TargetVar_set)-1
    plot(abs(b_2(i,2:5)),'Color',col(i,:),'Marker',shp(i),'MarkerSize',6, 'LineWidth',3)
    hold on
end
xlim([0 5])
legend(TargetVar_set_Rsquare_2)
set(gca,'xtick',[ 1 2 3 4])
set(gca,'xticklabel', para_set)
title(['Linear Regression Improved'],'FontSize',12)
saveas(gcf, [outputDir 'LinearRegressionImproved' '.fig'])
saveas(gcf, [outputDir 'LinearRegressionImproved' '.jpg'])
close


