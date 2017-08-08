% Derivatives of Force Measurements

clear all


vr_set = input('Reduced velocity set to run (e.g. 5, 5p5,...) :','s');
ph = {'180n' '135n' '90n' '45n' '0' '45' '90' '135'};

% load vr6p5.mat
% sam's attempt
eval(['load vr' vr_set '.mat']);

%newCL3 = zeros(8*36,1);

[junk1,junk2] = find(isnan(CDv));
CDv(junk1,junk2) = 0;

for i = 1:8
    k = 1;
    for j = 1:6
    X(j,:,i) = allXad(i,k:k+5);
    Y(j,:,i) = allYad(i,k:k+5);
    newtheta(j,:,i) = alltheta(i,k:k+5);
    newCL3(j,:,i) = CL3(i,k:k+5);%./sqrt(2);
    newCL5(j,:,i) = CL5(i,k:k+5);%./sqrt(2);
    newCLv(j,:,i) = CLv(i,k:k+5);%./sqrt(2);
    newCLvfft(j,:,i) = CLvfft(i,k:k+5);%./sqrt(2);
    newCmy(j,:,i) = Cmy(i,k:k+5);%./sqrt(2);
    newCDv(j,:,i) = CDv(i,k:k+5);%./sqrt(2);
    newCmx(j,:,i) = Cmx(i,k:k+5);%./sqrt(2);
    newCDmean(j,:,i) = CDmean(i,k:k+5);
    newCD2(j,:,i) = CD2(i,k:k+5);%./sqrt(2);
    newCD4(j,:,i) = CD4(i,k:k+5);%./sqrt(2);
    if strcmp(vr_set,'4p5') == 1 | strcmp(vr_set,'5') == 1 | strcmp(vr_set,'5p5') == 1
        avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.18.^3*0.6858*0.0381);
    else
        avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
    end
    k = k+6;
    end
    %newCL3(j:j+35) = CL3(i,:);
    %newYad(j:j+35) = allYad(i,:);
    %newXad(j:j+35) = allXad(i,:);
    %newtheta(j:j+35) = alltheta(i,:);
    %j = j+36;
end

X(:,:,9) = X(:,:,1);
Y(:,:,9) = Y(:,:,1);
newtheta(:,:,9) = 180;
newCLv(:,:,9) = newCLv(:,:,1);
newCDv(:,:,9) = newCDv(:,:,1);
avgPow(:,:,9) = avgPow(:,:,1);

sp_avgPow(:,:,1) = avgPow(:,:,8);
sp_avgPow(:,:,2:10) = avgPow(:,:,1:9);
sp_avgPow(:,:,11) = avgPow(:,:,2);
sp_CLv(:,:,1) = newCLv(:,:,8);
sp_CLv(:,:,2:10) = newCLv(:,:,1:9);
sp_CLv(:,:,11) = newCLv(:,:,2);
sp_CDv(:,:,1) = newCDv(:,:,8);
sp_CDv(:,:,2:10) = newCDv(:,:,1:9);
sp_CDv(:,:,11) = newCDv(:,:,2);

aa = newCLv(:,3,1);
bb = Y(:,3,1);
cc = newtheta(:,3,1)*pi/180;
aaa = newCDv(:,3,1);
aaaa = avgPow(:,3,1);

for i = 1:8
    aa = [aa newCLv(:,3,i+1)];
    bb = [bb Y(:,3,i+1)];
    cc = [cc newtheta(:,3,i+1)*pi/180];
    aaa = [aaa newCDv(:,3,i+1)];
    aaaa = [aaaa avgPow(:,3,i+1)];
end

a5 = sp_avgPow(:,3,1);
a6 = sp_CLv(:,3,1);
a7 = sp_CDv(:,3,1);

for i = 1:10
    a5 = [a5 sp_avgPow(:,3,i+1)];
    a6 = [a6 sp_CLv(:,3,i+1)];
    a7 = [a7 sp_CDv(:,3,i+1)];
end

dd = newCLv(4,:,1);
ee = X(4,:,1);
ff = newtheta(4,:,1)*pi/180;
ddd = newCDv(4,:,1);
dddd = avgPow(4,:,1);

for i = 1:8
    dd = [dd; newCLv(4,:,i+1)];
    ee = [ee; X(4,:,i+1)];
    ff = [ff; newtheta(4,:,i+1)*pi/180];
    ddd = [ddd; newCDv(4,:,i+1)];
    dddd = [dddd; avgPow(4,:,i+1)];
end

d5 = sp_avgPow(4,:,1);
d6 = sp_CLv(4,:,1);
d7 = sp_CDv(4,:,1);

for i = 1:10
    d5 = [d5; sp_avgPow(4,:,i+1)];
    d6 = [d6; sp_CLv(4,:,i+1)];
    d7 = [d7; sp_CDv(4,:,i+1)];
end

[P_clv_ax,P_clv_ay,P_clv_theta] = gradient(newCLv,0.15,0.25,pi/4);
[P1,P2] = gradient(newCLv(:,:,5),0.15,0.25);
[P3,P4] = gradient(aa,pi/4,0.25);
[P5,P6] = gradient(dd,0.15,pi/4);
[PP1,PP2] = gradient(newCDv(:,:,5),0.15,0.25);
[PP3,PP4] = gradient(aaa,pi/4,0.25);
[PP5,PP6] = gradient(ddd,0.15,pi/4);
[PPP1,PPP2] = gradient(avgPow(:,:,5),0.15,0.25);
[PPP3,PPP4] = gradient(aaaa,pi/4,0.25);
[PPP5,PPP6] = gradient(dddd,0.15,pi/4);
[PS3,PS4] = gradient(a5,pi/4,0.25);
[PS5,PS6] = gradient(d5,0.15,pi/4);
[PPS3,PPS4] = gradient(a6,pi/4,0.25);
[PPS5,PPS6] = gradient(d6,0.15,pi/4);
[PPPS3,PPPS4] = gradient(a7,pi/4,0.25);
[PPPS5,PPPS6] = gradient(d7,0.15,pi/4);
%[P_cdv_ax,P_cdv_ay,P_cdv_theta] = gradient(newCDv,0.15,0.25,pi/4);

 conts = [-4 -2 0 2 4];
 cols = {'m' 'b' 'g' 'y' 'r'};

figure(1);clf;hold on;
[f1,h1] = contourf(X(:,:,5),Y(:,:,5),P1);
plot(0.31,0.91,'mp','MarkerSize',10);
%quiver(X(:,:,5),Y(:,:,5),P1./sqrt(P1.^2+P2.^2),P2./sqrt(P1.^2+P2.^2),0.2);
%quiver(X(:,:,5),Y(:,:,5),P1,P2,0.2);
clabel(f1,h1);
xlabel('A_{x}/D');
ylabel('A_{y}/D');
title('Contours of d(CLv)/d(Ax*)')
axis([0 0.75 0.25 1.5])
%axis image

figure(2);clf;hold on;
[f1,h1] = contourf(X(:,:,5),Y(:,:,5),P2);
plot(0.31,0.91,'mp','MarkerSize',10);
%quiver(X(:,:,5),Y(:,:,5),P1./sqrt(P1.^2+P2.^2),P2./sqrt(P1.^2+P2.^2),0.2);
%quiver(X(:,:,5),Y(:,:,5),P1,P2,0.2);
clabel(f1,h1);
xlabel('A_{x}/D');
ylabel('A_{y}/D');
title('Contours of d(CLv)/d(Ay*)')
axis([0 0.75 0.25 1.5])
%axis image

figure(3);clf;hold on;
[f2,h2] = contourf(cc,bb,PPS3(:,2:10));
plot(0,0.91,'mp','MarkerSize',10);
%quiver(cc,bb,P3./sqrt(P3.^2+P4.^2),P4./sqrt(P3.^2+P4.^2),0.2);
%quiver(cc,bb,P3,P4,0.2);
clabel(f2,h2);
xlabel('\theta');
ylabel('A_{y}/D');
title('Contours of d(CLv)/d(\theta)')
axis([-pi pi 0.25 1.5])
%axis image

figure(4);clf;hold on;
[f2,h2] = contourf(cc,bb,PPS4(:,2:10));
plot(0,0.91,'mp','MarkerSize',10);
%quiver(cc,bb,P3./sqrt(P3.^2+P4.^2),P4./sqrt(P3.^2+P4.^2),0.2);
%quiver(cc,bb,P3,P4,0.2);
clabel(f2,h2);
xlabel('\theta');
ylabel('A_{y}/D');
title('Contours of d(CLv)/d(Ay*)')
axis([-pi pi 0.25 1.5])
%axis image

figure(5);clf;hold on;
[f3,h3] = contourf(ee,ff,PPS5(2:10,:));
plot(0.31,0,'mp','MarkerSize',10);
%quiver(ee,ff,P5./sqrt(P5.^2+P6.^2),P6./sqrt(P5.^2+P6.^2),0.2);
%quiver(ee,ff,P5,P6,0.2);
clabel(f3,h3);
xlabel('A_{x}/D');
ylabel('\theta');
title('Contours of d(CLv)/d(Ax*)')
axis([0 0.75 -pi pi])
%axis image

figure(6);clf;hold on;
[f3,h3] = contourf(ee,ff,PPS6(2:10,:));
plot(0.31,0,'mp','MarkerSize',10);
%quiver(ee,ff,P5./sqrt(P5.^2+P6.^2),P6./sqrt(P5.^2+P6.^2),0.2);
%quiver(ee,ff,P5,P6,0.2);
clabel(f3,h3);
xlabel('A_{x}/D');
ylabel('\theta');
title('Contours of d(CLv)/d(\theta)')
axis([0 0.75 -pi pi])
%axis image

figure(7);clf;hold on;
[f1,h1] = contourf(X(:,:,5),Y(:,:,5),PP1);
plot(0.31,0.91,'mp','MarkerSize',10);
%quiver(X(:,:,5),Y(:,:,5),P1./sqrt(P1.^2+P2.^2),P2./sqrt(P1.^2+P2.^2),0.2);
%quiver(X(:,:,5),Y(:,:,5),P1,P2,0.2);
clabel(f1,h1);
xlabel('A_{x}/D');
ylabel('A_{y}/D');
title('Contours of d(CDv)/d(Ax*)')
axis([0 0.75 0.25 1.5])
%axis image

figure(8);clf;hold on;
[f1,h1] = contourf(X(:,:,5),Y(:,:,5),PP2);
plot(0.31,0.91,'mp','MarkerSize',10);
%quiver(X(:,:,5),Y(:,:,5),P1./sqrt(P1.^2+P2.^2),P2./sqrt(P1.^2+P2.^2),0.2);
%quiver(X(:,:,5),Y(:,:,5),P1,P2,0.2);
clabel(f1,h1);
xlabel('A_{x}/D');
ylabel('A_{y}/D');
title('Contours of d(CDv)/d(Ay*)')
axis([0 0.75 0.25 1.5])
%axis image

figure(9);clf;hold on;
[f2,h2] = contourf(cc,bb,PPPS3(:,2:10));
plot(0,0.91,'mp','MarkerSize',10);
%quiver(cc,bb,P3./sqrt(P3.^2+P4.^2),P4./sqrt(P3.^2+P4.^2),0.2);
%quiver(cc,bb,P3,P4,0.2);
clabel(f2,h2);
xlabel('\theta');
ylabel('A_{y}/D');
title('Contours of d(CDv)/d(\theta)')
axis([-pi pi 0.25 1.5])
%axis image

figure(10);clf;hold on;
[f2,h2] = contourf(cc,bb,PPPS4(:,2:10));
plot(0,0.91,'mp','MarkerSize',10);
%quiver(cc,bb,P3./sqrt(P3.^2+P4.^2),P4./sqrt(P3.^2+P4.^2),0.2);
%quiver(cc,bb,P3,P4,0.2);
clabel(f2,h2);
xlabel('\theta');
ylabel('A_{y}/D');
title('Contours of d(CDv)/d(Ay*)')
axis([-pi pi 0.25 1.5])
%axis image

figure(11);clf;hold on;
[f3,h3] = contourf(ee,ff,PPPS5(2:10,:));
plot(0.31,0,'mp','MarkerSize',10);
%quiver(ee,ff,P5./sqrt(P5.^2+P6.^2),P6./sqrt(P5.^2+P6.^2),0.2);
%quiver(ee,ff,P5,P6,0.2);
clabel(f3,h3);
xlabel('A_{x}/D');
ylabel('\theta');
title('Contours of d(CDv)/d(Ax*)')
axis([0 0.75 -pi pi])
%axis image

figure(12);clf;hold on;
[f3,h3] = contourf(ee,ff,PPPS6(2:10,:));
plot(0.31,0,'mp','MarkerSize',10);
%quiver(ee,ff,P5./sqrt(P5.^2+P6.^2),P6./sqrt(P5.^2+P6.^2),0.2);
%quiver(ee,ff,P5,P6,0.2);
clabel(f3,h3);
xlabel('A_{x}/D');
ylabel('\theta');
title('Contours of d(CDv)/d(\theta)')
axis([0 0.75 -pi pi])
%axis image


% figure(1);clf;hold on;
% [f1,h1] = contourf(X(:,:,5),Y(:,:,5),PPP1);
% plot(0.31,0.91,'rp','MarkerSize',10);
% %quiver(X(:,:,5),Y(:,:,5),P1./sqrt(P1.^2+P2.^2),P2./sqrt(P1.^2+P2.^2),0.2);
% %quiver(X(:,:,5),Y(:,:,5),P1,P2,0.2);
% clabel(f1,h1);
% xlabel('A_{x}/D');
% ylabel('A_{y}/D');
% title('Contours of d(Cp)/d(Ax*)')
% axis([0 0.75 0.25 1.5])
% %axis image
% 
% figure(2);clf;hold on;
% [f1,h1] = contourf(X(:,:,5),Y(:,:,5),PPP2);
% plot(0.31,0.91,'rp','MarkerSize',10);
% %quiver(X(:,:,5),Y(:,:,5),P1./sqrt(P1.^2+P2.^2),P2./sqrt(P1.^2+P2.^2),0.2);
% %quiver(X(:,:,5),Y(:,:,5),P1,P2,0.2);
% clabel(f1,h1);
% xlabel('A_{x}/D');
% ylabel('A_{y}/D');
% title('Contours of d(Cp)/d(Ay*)')
% axis([0 0.75 0.25 1.5])
% %axis image
% 
% figure(3);clf;hold on;
% [f2,h2] = contourf(cc,bb,PS3(:,2:10));
% plot(0,0.91,'rp','MarkerSize',10);
% %quiver(cc,bb,P3./sqrt(P3.^2+P4.^2),P4./sqrt(P3.^2+P4.^2),0.2);
% %quiver(cc,bb,P3,P4,0.2);
% clabel(f2,h2);
% xlabel('\theta');
% ylabel('A_{y}/D');
% title('Contours of d(Cp)/d(\theta)')
% axis([-pi pi 0.25 1.5]);
% %axis image
% 
% figure(4);clf;hold on;
% [f2,h2] = contourf(cc,bb,PS4(:,2:10));
% plot(0,0.91,'rp','MarkerSize',10);
% %quiver(cc,bb,P3./sqrt(P3.^2+P4.^2),P4./sqrt(P3.^2+P4.^2),0.2);
% %quiver(cc,bb,P3,P4,0.2);
% clabel(f2,h2);
% xlabel('\theta');
% ylabel('A_{y}/D');
% title('Contours of d(Cp)/d(Ay*)')
% axis([-pi pi 0.25 1.5]);
% %axis image
% 
% figure(5);clf;hold on;
% [f3,h3] = contourf(ee,ff,PS5(2:10,:));
% plot(0.31,0,'rp','MarkerSize',10);
% %quiver(ee,ff,P5./sqrt(P5.^2+P6.^2),P6./sqrt(P5.^2+P6.^2),0.2);
% %quiver(ee,ff,P5,P6,0.2);
% clabel(f3,h3);
% xlabel('A_{x}/D');
% ylabel('\theta');
% title('Contours of d(Cp)/d(Ax)')
% axis([0 0.75 -pi pi]);
% %axis image
% 
% figure(6);clf;hold on;
% [f3,h3] = contourf(ee,ff,PS6(2:10,:));
% plot(0.31,0,'rp','MarkerSize',10);
% %quiver(ee,ff,P5./sqrt(P5.^2+P6.^2),P6./sqrt(P5.^2+P6.^2),0.2);
% %quiver(ee,ff,P5,P6,0.2);
% clabel(f3,h3);
% xlabel('A_{x}/D');
% ylabel('\theta');
% title('Contours of d(Cp)/d(\theta)')
% axis([0 0.75 -pi pi]);
% %axis image

% figure(4);clf;hold on;
% [f4,h4] = contourf(X(:,:,5),Y(:,:,5),PP2);
% plot(0.31,0.91,'rp','MarkerSize',10);
% %quiver(X(:,:,5),Y(:,:,5),PP1./sqrt(PP1.^2+PP2.^2),PP2./sqrt(PP1.^2+PP2.^2),0.2);
% %quiver(X(:,:,5),Y(:,:,5),PP1,PP2,0.2);
% clabel(f4,h4);
% xlabel('A_{x}/D');
% ylabel('A_{y}/D');
% title('Contours of C_{Dv} with gradient')
% %axis image
% 
% figure(5);clf;hold on;
% [f5,h5] = contourf(cc,bb,PP4);
% plot(0,0.91,'rp','MarkerSize',10);
% %quiver(cc,bb,PP3./sqrt(PP3.^2+PP4.^2),PP4./sqrt(PP3.^2+PP4.^2),0.2);
% %quiver(cc,bb,PP3,PP4,0.2);
% clabel(f5,h5);
% xlabel('\theta');
% ylabel('A_{y}/D');
% title('Contours of C_{Dv} with gradient')
% %axis image
% 
% figure(6);clf;hold on;
% [f6,h6] = contourf(ee,ff,PP6);
% plot(0.31,0,'rp','MarkerSize',10);
% %quiver(ee,ff,PP5./sqrt(PP5.^2+PP6.^2),PP6./sqrt(PP5.^2+PP6.^2),0.2);
% %quiver(ee,ff,PP5,PP6,0.2);
% clabel(f6,h6);
% xlabel('A_{x}/D');
% ylabel('\theta');
% title('Contours of C_{Dv} with gradient')
% %axis image

% figure(7);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta*pi/180,P_clv_ax,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
%     
%     %p = sqrt(P_cdv_ax.^2 + P_cdv_ay.^2+P_cdv_theta.^2);
%     %daspect([1,1,5]);
%     %h2 = coneplot(X,Y,newtheta*pi/180,P_clv_ax,P_clv_ay,P_clv_theta,X,Y,newtheta*pi/180,1);
% %     isonormals(X,Y,newtheta,newCL3,h);
% %     vertices = get(h,'Vertices');
% %     vertexnormals = get(h,'VertexNormals');
% %     quiver3(vertices(:,1),vertices(:,2),vertices(:,3),vertexnormals(:,1),vertexnormals(:,2),vertexnormals(:,3),0.005);
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% axis([0 0.75 0 1.5 -pi pi]);
% %title('CL_{3}','FontSize',16);
% view(-15,20);
% grid on
% legend(num2str(conts'));
% 
% %edited by sam
% %print(gcf,'-depsc','CL3_vr6p5');
% eval(['print(gcf,''-depsc'',''CL3_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);
% 
% conts = [0.2 0.6 1];
% cols = {'m' 'b' 'g'};

% figure(2);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCL5,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('CL_{5}','FontSize',16);
% view(-15,20);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','CL5_vr6p5');
% eval(['print(gcf,''-depsc'',''CL5_vr' vr_set ''');'])
% eval(['print(gcf,''-djpeg100'',''CL5_vr' vr_set ''');'])
% 
% conts = [-2 -1 0 1 2];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(3);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCLv,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('CL_{v}','FontSize',16);
% view(-60,30);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','CLv_vr6p5');
% eval(['print(gcf,''-depsc'',''CLv_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''CLv_vr' vr_set ''');']);
% 
% conts = [-2 -1 0 1 2];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(12);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCLvfft,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('CL_{v}','FontSize',16);
% view(-60,30);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','CLv_vr6p5');
% eval(['print(gcf,''-depsc'',''CLvfft_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''CLvfft_vr' vr_set ''');']);
% 
% conts = [-2 -1 0 1 2];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(4);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCmy,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('C_{my}','FontSize',16);
% view(-30,20);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','Cmy_vr6p5');
% eval(['print(gcf,''-depsc'',''Cmy_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''Cmy_vr' vr_set ''');']);
% 
% conts = [2 3 4 5 6];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(5);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCDmean,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('Mean CD','FontSize',16);
% view(-15,20);
% grid on
% legend('2','3','4','5','6');
% 
% % edited by sam
% %print(gcf,'-depsc','CDmean_vr6p5');
% eval(['print(gcf,''-depsc'',''CDmean_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''CDmean_vr' vr_set ''');']);
% 
% conts = [0.5 1 1.5 2 2.5];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(6);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCD2,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('Fluctuating CD_{2}','FontSize',16);
% view(-30,30);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','CD2_vr6p5');
% eval(['print(gcf,''-depsc'',''CD2_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''CD2_vr' vr_set ''');']);
% 
% conts = [0.2 0.6 1];
% cols = {'m' 'b' 'g'};
% 
% figure(7);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCD4,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('Fluctuating CD_{4}','FontSize',16);
% view(-15,20);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','CD4_vr6p5');
% eval(['print(gcf,''-depsc'',''CD4_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''CD4_vr' vr_set ''');']);
% 
% 
% conts = [-3 -2 -1 0 1];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(8);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCDv,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('CD_{v}','FontSize',16);
% view(-15,20);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','CDv_vr6p5');
% eval(['print(gcf,''-depsc'',''CDv_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''CDv_vr' vr_set ''');']);
% 
% conts = [-2 -1 0 1 2];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(9);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCmx,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('C_{mx}','FontSize',16);
% view(-30,30);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','Cmx_vr6p5');
% eval(['print(gcf,''-depsc'',''Cmx_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''Cmx_vr' vr_set ''');']);
% 
% conts = [-2 -1 0 0.1 0.3];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(10);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,avgPow,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('Normalized Total Avg Power','FontSize',16);
% view(-60,30);
% grid on
% legend(num2str(conts'));
% 
% % edited by sam
% %print(gcf,'-depsc','Pow_vr6p5');
% eval(['print(gcf,''-depsc'',''Pow_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''Pow_vr' vr_set ''');']);
% 
% conts = [-2 -1 0 0.1 0.3];
% cols = {'m' 'b' 'g' 'y' 'r'};
% 
% figure(11);hold on;
% h1 = patch(isosurface(X,Y,newtheta,avgPow,0));
% set(h1,'FaceColor','k')
% set(h1,'FaceAlpha',0.5)
% h2 = patch(isosurface(X,Y,newtheta,newCL3,0.5));
% set(h2,'FaceColor','m')
% set(h2,'FaceAlpha',0.5)
% h3 = patch(isosurface(X,Y,newtheta,newCL3,1));
% set(h3,'FaceColor','b')
% set(h3,'FaceAlpha',0.5)
% h5 = patch(isosurface(X,Y,newtheta,newCL3,2));
% set(h5,'FaceColor','r')
% set(h5,'FaceAlpha',0.5)
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% %title('CL_{3} with excitation region','FontSize',16);
% view(-45,50);
% grid on
% legend('Excitation Region','CL_{3} = 0.5','CL_{3} = 1','CL_{3} = 2');
% 
% % edited by sam
% %print(gcf,'-depsc','Exc_and_CL3_vr6p5');
% eval(['print(gcf,''-depsc'',''Exc_and_CL3_vr' vr_set ''');']);
% eval(['print(gcf,''-djpeg100'',''Exc_and_CL3_vr' vr_set ''');']);

% for i = 1:8
% 
% figure(2);clf;hold on;
% [C,h] = contour(X(:,:,i),Y(:,:,i),newCL3(:,:,i));
% xlabel('A_{x}/D','FontSize',18);
% ylabel('A_{y}/D','FontSize',18);
% axis([0 0.75 0.25 1.5]);
% clabel(C,h,'FontSize',16);
% set(gca,'FontSize',16);
% %title('CL_{3}','FontSize',16);
% %view(-15,20);
% %grid on
% %legend(num2str(conts'));
% 
% %edited by sam
% %print(gcf,'-depsc','CL3_vr6p5');
% %eval(['print(gcf,''-depsc'',''contours/CL3_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% %eval(['print(gcf,''-djpeg100'',''contours/CL3_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% 
% figure(3);clf;hold on;
% [C,h] = contour(X(:,:,i),Y(:,:,i),newCLv(:,:,i));
% xlabel('A_{x}/D','FontSize',18);
% ylabel('A_{y}/D','FontSize',18);
% axis([0 0.75 0.25 1.5]);
% clabel(C,h,'FontSize',16);
% set(gca,'FontSize',16);
% 
% %eval(['print(gcf,''-depsc'',''contours/CLv_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% %eval(['print(gcf,''-djpeg100'',''contours/CLv_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% 
% figure(4);clf;hold on;
% [C,h] = contour(X(:,:,i),Y(:,:,i),newCDv(:,:,i));
% xlabel('A_{x}/D','FontSize',18);
% ylabel('A_{y}/D','FontSize',18);
% axis([0 0.75 0.25 1.5]);
% clabel(C,h,'FontSize',16);
% set(gca,'FontSize',16);
% 
% %eval(['print(gcf,''-depsc'',''contours/CDv_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% %eval(['print(gcf,''-djpeg100'',''contours/CDv_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% 
% figure(5);clf;hold on;
% [C,h] = contour(X(:,:,i),Y(:,:,i),avgPow(:,:,i));
% xlabel('A_{x}/D','FontSize',18);
% ylabel('A_{y}/D','FontSize',18);
% axis([0 0.75 0.25 1.5]);
% clabel(C,h,'FontSize',16);
% set(gca,'FontSize',16);
% 
% %eval(['print(gcf,''-depsc'',''contours/Pow_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% %eval(['print(gcf,''-djpeg100'',''contours/Pow_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% 
% figure(6);clf;hold on;
% [C,h] = contour(X(:,:,i),Y(:,:,i),newCmy(:,:,i));
% xlabel('A_{x}/D','FontSize',18);
% ylabel('A_{y}/D','FontSize',18);
% axis([0 0.75 0.25 1.5]);
% clabel(C,h,'FontSize',16);
% set(gca,'FontSize',16);
% 
% %eval(['print(gcf,''-depsc'',''contours/Cmy_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% %eval(['print(gcf,''-djpeg100'',''contours/Cmy_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% 
% figure(7);clf;hold on;
% [C,h] = contour(X(:,:,i),Y(:,:,i),newCmx(:,:,i));
% xlabel('A_{x}/D','FontSize',18);
% ylabel('A_{y}/D','FontSize',18);
% axis([0 0.75 0.25 1.5]);
% clabel(C,h,'FontSize',16);
% set(gca,'FontSize',16);
% 
% %eval(['print(gcf,''-depsc'',''contours/Cmx_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% %eval(['print(gcf,''-djpeg100'',''contours/Cmx_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
% 
% end