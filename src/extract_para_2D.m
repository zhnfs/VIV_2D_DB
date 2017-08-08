function [value, vr, theta, X, Y] = extract_para_2D(databaseFileName, thetaIndex, AmpXIndex, AmpYIndex, targetpar, theta_set, X_set, Y_set, all_Y)
% clear
% plot_para_2D('theta',[1:2],' vr',[1:4],'CL1',1)

% vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
% theta_set = {'-pi','-3/4pi','-1/2pi','-1/4pi','0','1/4pi','1/2pi','3/4pi','pi' };
% X_set = {'0','0.15','0.3','0.45','0.6','0.75'};
% Y_set = {'0.25','0.5','0.75','1','1.25','1.5'};
% color_set = {'r','y','g','c','b','m','k',[255/255 140/255 0],'r'};
% targetpar_set = {'CL1','CL3','CD2','Cmy', 'Cmx','pow','CDv'};

newX = NaN(length(Y_set),length(X_set),length(theta_set));
newY = NaN(length(Y_set),length(X_set),length(theta_set));
newtheta = NaN(length(Y_set),length(X_set),length(theta_set));
eval(['new'  targetpar '= NaN(length(Y_set),length(X_set),length(theta_set));']);

load(databaseFileName); 

[ind1,ind2] = find(isnan(Cmx));
Cmx(ind1,ind2) = 1;%% just for ploting reason

for i = 1:length(theta_set)
   for j = 1:length(all_Y)
        ind = find(allYad(i,:) == all_Y(j));
        for k = ind
            newX(j,1:length(ind),i) = allXad(i,ind);
            newY(j,1:length(ind),i) = allYad(i,ind);
            newtheta(j,1:length(ind),i) = alltheta(i,ind);
            newvr(j,1:length(ind),i) = Vr(i,ind);
            eval(['new'  targetpar '(j,1:length(ind),i) = ' targetpar '(i,ind);']);
           
        end
   end    

end

X = newX(AmpYIndex, AmpXIndex, thetaIndex);
Y = newY(AmpYIndex, AmpXIndex, thetaIndex);
vr = newvr(AmpYIndex, AmpXIndex, thetaIndex);
theta = newtheta(AmpYIndex, AmpXIndex, thetaIndex);

eval(['value = new' targetpar '(AmpYIndex, AmpXIndex, thetaIndex);']);
