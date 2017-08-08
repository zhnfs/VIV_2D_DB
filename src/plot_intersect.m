% Plot intersection of added mass surfaces with excitation region

clear all;



names = {'vr6'};
cols = {'m' 'b' 'g' 'y' 'r'};

%eval(['load ' char(names(m))]);
load vr6
    for i = 1:8
        k = 1;
        for j = 1:6
        X(j,:,i) = allXad(i,k:k+5);
        Y(j,:,i) = allYad(i,k:k+5);
        newtheta(j,:,i) = alltheta(i,k:k+5);
        newCL3(j,:,i) = CL3(i,k:k+5);
        newCL5(j,:,i) = CL5(i,k:k+5);
        newCLv(j,:,i) = CLv(i,k:k+5);
        newCmy(j,:,i) = Cmy(i,k:k+5);
        newCDv(j,:,i) = CDv(i,k:k+5);
        newCmx(j,:,i) = Cmx(i,k:k+5);
        newCDmean(j,:,i) = CDmean(i,k:k+5);
        newCD2(j,:,i) = CD2(i,k:k+5);
        newCD4(j,:,i) = CD4(i,k:k+5);
        avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
        k = k+6;
        end
    end

% X = interp3(X,1);
% Y = interp3(Y,1);
% newtheta = interp3(newtheta,1);
% avgPow = interp3(avgPow,1);
% newCmy = interp3(newCmy,1);
% newCmx = interp3(newCmx,1);

figure(1);hold on;
h = patch(isosurface(X,Y,newtheta,avgPow,0));
set(h,'FaceColor','g');
set(h,'FaceAlpha',0.5)
h2 = patch(isosurface(X,Y,newtheta,newCmy,0.33));
set(h2,'FaceColor','r');
set(h2,'FaceAlpha',0.5)
h3 = patch(isosurface(X,Y,newtheta,newCmx,0.2929));
set(h3,'FaceColor','b');
set(h3,'FaceAlpha',0.5)
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('\theta','FontSize',16);
%title('Avg Pow','FontSize',16);
view(-15,20);
grid on
legend('Avg Pow = 0', 'Cmy = 0.33', 'Cmx = 0.2929');

print(gcf,'-depsc','Vr6_intersect');

% Powvert = get(h,'Vertices');
% Phasepow = Powvert(:,3)./180;
% Powvert2 = [Powvert(:,1) Powvert(:,2) Phasepow];
% Cmyvert = get(h2,'Vertices');
% Phasecmy = Cmyvert(:,3)./180;
% Cmyvert2 = [Cmyvert(:,1) Cmyvert(:,2) Phasecmy];
% Cmxvert = get(h3,'Vertices');
% Phasecmx = Cmxvert(:,3)./180;
% Cmxvert2 = [Cmxvert(:,1) Cmxvert(:,2) Phasecmx];
% 
% d_pow_cmy = zeros(length(Powvert2),length(Cmyvert2));
% d_pow_cmx = zeros(length(Powvert2),length(Cmxvert2));
% d_cmy_cmx = zeros(length(Cmyvert2),length(Cmxvert2));
% 
% for i = 1:length(Powvert2)
%     for j = 1:length(Cmyvert2)
%         for k = 1:length(Cmxvert2)
%             d(i,j,k) = sqrt(sum((Powvert2(i,:) - Cmyvert2(j,:)).^2)) + ...
%                 sqrt(sum((Powvert2(i,:) - Cmxvert2(k,:)).^2));
%         end
%     end
% end
% 
% g = min(min(min(d)));
% ind = find(d == g);
% 
% [indPow,indCmy,indCmx] = ind2sub(size(d),ind);
% 
% prediction = (Powvert(indPow,:) + Cmyvert(indCmy,:) + Cmxvert(indCmx,:))/3;
