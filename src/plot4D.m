% Make 3D scatter plot of data

clear all

% Sam added an input command to this file
vr_set = input('Reduced velocity set to run (e.g. 5, 5p5,...) :','s');
ph = {'180n' '135n' '90n' '45n' '0' '45' '90' '135'};

% load vr6p5.mat
% sam's attempt
eval(['load vr' vr_set '.mat']);

%newCL3 = zeros(8*36,1);

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

% conts = [0.5 1 2 3 4];
% cols = {'m' 'b' 'g' 'y' 'r'};

% figure(1);hold on;
% for i = 1:length(conts)
%     h = patch(isosurface(X,Y,newtheta,newCL3,conts(i)));
%     eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
%     set(h,'FaceAlpha',0.5)
% %     isonormals(X,Y,newtheta,newCL3,h);
% %     vertices = get(h,'Vertices');
% %     vertexnormals = get(h,'VertexNormals');
% %     quiver3(vertices(:,1),vertices(:,2),vertices(:,3),vertexnormals(:,1),vertexnormals(:,2),vertexnormals(:,3),0.005);
% end
% xlabel('A_{x}/D','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% axis([0 0.75 0 1.5 -180 180]);
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

for i = 1:8

figure(2);clf;hold on;
[C,h] = contour(X(:,:,i),Y(:,:,i),newCL3(:,:,i));
xlabel('A_{x}/D','FontSize',18);
ylabel('A_{y}/D','FontSize',18);
axis([0 0.75 0.25 1.5]);
clabel(C,h,'FontSize',16);
set(gca,'FontSize',16);
%title('CL_{3}','FontSize',16);
%view(-15,20);
%grid on
%legend(num2str(conts'));

%edited by sam
%print(gcf,'-depsc','CL3_vr6p5');
eval(['print(gcf,''-depsc'',''contours/CL3_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
%eval(['print(gcf,''-djpeg100'',''contours/CL3_cont_vr' vr_set '_ph' char(ph(i)) ''');']);

figure(3);clf;hold on;
[C,h] = contour(X(:,:,i),Y(:,:,i),newCLv(:,:,i));
xlabel('A_{x}/D','FontSize',18);
ylabel('A_{y}/D','FontSize',18);
axis([0 0.75 0.25 1.5]);
clabel(C,h,'FontSize',16);
set(gca,'FontSize',16);

eval(['print(gcf,''-depsc'',''contours/CLv_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
%eval(['print(gcf,''-djpeg100'',''contours/CLv_cont_vr' vr_set '_ph' char(ph(i)) ''');']);

figure(4);clf;hold on;
[C,h] = contour(X(:,:,i),Y(:,:,i),newCDv(:,:,i));
xlabel('A_{x}/D','FontSize',18);
ylabel('A_{y}/D','FontSize',18);
axis([0 0.75 0.25 1.5]);
clabel(C,h,'FontSize',16);
set(gca,'FontSize',16);

eval(['print(gcf,''-depsc'',''contours/CDv_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
%eval(['print(gcf,''-djpeg100'',''contours/CDv_cont_vr' vr_set '_ph' char(ph(i)) ''');']);

figure(5);clf;hold on;
[C,h] = contour(X(:,:,i),Y(:,:,i),avgPow(:,:,i));
xlabel('A_{x}/D','FontSize',18);
ylabel('A_{y}/D','FontSize',18);
axis([0 0.75 0.25 1.5]);
clabel(C,h,'FontSize',16);
set(gca,'FontSize',16);

eval(['print(gcf,''-depsc'',''contours/Pow_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
%eval(['print(gcf,''-djpeg100'',''contours/Pow_cont_vr' vr_set '_ph' char(ph(i)) ''');']);

figure(6);clf;hold on;
[C,h] = contour(X(:,:,i),Y(:,:,i),newCmy(:,:,i));
xlabel('A_{x}/D','FontSize',18);
ylabel('A_{y}/D','FontSize',18);
axis([0 0.75 0.25 1.5]);
clabel(C,h,'FontSize',16);
set(gca,'FontSize',16);

eval(['print(gcf,''-depsc'',''contours/Cmy_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
%eval(['print(gcf,''-djpeg100'',''contours/Cmy_cont_vr' vr_set '_ph' char(ph(i)) ''');']);

figure(7);clf;hold on;
[C,h] = contour(X(:,:,i),Y(:,:,i),newCmx(:,:,i));
xlabel('A_{x}/D','FontSize',18);
ylabel('A_{y}/D','FontSize',18);
axis([0 0.75 0.25 1.5]);
clabel(C,h,'FontSize',16);
set(gca,'FontSize',16);

eval(['print(gcf,''-depsc'',''contours/Cmx_cont_vr' vr_set '_ph' char(ph(i)) ''');']);
%eval(['print(gcf,''-djpeg100'',''contours/Cmx_cont_vr' vr_set '_ph' char(ph(i)) ''');']);

end