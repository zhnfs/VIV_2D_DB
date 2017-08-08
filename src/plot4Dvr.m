% Plot 4D data as a function of Vr, holding phase constant

clear all;

vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};

for p = 1:length(vr_set)
    eval(['load vr' char(vr_set(p)) '.mat']);

    for i = 1:8
        k = 1;
        for j = 1:6
        X(j,:,i) = allXad(i,k:k+5);
        Y(j,:,i) = allYad(i,k:k+5);
        avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
        vr(j,:,i) = Vr(i,k:k+5);
        k = k+6;
        end
        %newCL3(j:j+35) = CL3(i,:);
        %newYad(j:j+35) = allYad(i,:);
        %newXad(j:j+35) = allXad(i,:);
        %newtheta(j:j+35) = alltheta(i,:);
        %j = j+36;
    end
    vr180n(:,:,p) = vr(:,:,1);
    vr135n(:,:,p) = vr(:,:,2);
    vr90n(:,:,p) = vr(:,:,3);
    vr45n(:,:,p) = vr(:,:,4);
    vr0(:,:,p) = vr(:,:,5);
    vr45(:,:,p) = vr(:,:,6);
    vr90(:,:,p) = vr(:,:,7);
    vr135(:,:,p) = vr(:,:,8);
    avgPow180n(:,:,p) = avgPow(:,:,1);
    avgPow135n(:,:,p) = avgPow(:,:,2);
    avgPow90n(:,:,p) = avgPow(:,:,3);
    avgPow45n(:,:,p) = avgPow(:,:,4);
    avgPow0(:,:,p) = avgPow(:,:,5);
    avgPow45(:,:,p) = avgPow(:,:,6);
    avgPow90(:,:,p) = avgPow(:,:,7);
    avgPow135(:,:,p) = avgPow(:,:,8);
end

%conts = [-2 -1 0 0.1 0.3];
conts = [0];
%cols = {'m' 'b' 'g' 'y' 'r'};
cols = {'g'};

figure(1);hold on;
for i = 1:length(conts)
    h = patch(isosurface(X,Y,vr180n,avgPow180n,conts(i)));
    eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
    set(h,'FaceAlpha',0.5)
end
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('V_{r}','FontSize',16);
title('Normalized Avg Power, \theta = -180','FontSize',16);
axis([0 0.75 0 1.5 4.5 8])
view(-40,20);
grid on
legend(num2str(conts'));

%movieims(1) = getframe(gcf);
% 
% %edited by sam
print(gcf,'-depsc','avgpow_180n');
print(gcf,'-djpeg100','avgpow_180n');
% %eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);

figure(2);hold on;
for i = 1:length(conts)
    h = patch(isosurface(X,Y,vr135n,avgPow135n,conts(i)));
    eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
    set(h,'FaceAlpha',0.5)
end
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('V_{r}','FontSize',16);
title('Normalized Avg Power, \theta = -135','FontSize',16);
axis([0 0.75 0 1.5 4.5 8])
view(-40,20);
grid on
legend(num2str(conts'));

%movieims(2) = getframe(gcf);
% 
% %edited by sam
print(gcf,'-depsc','avgpow_135n');
print(gcf,'-djpeg100','avgpow_135n');
% %eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);

figure(3);hold on;
for i = 1:length(conts)
    h = patch(isosurface(X,Y,vr90n,avgPow90n,conts(i)));
    eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
    set(h,'FaceAlpha',0.5)
end
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('V_{r}','FontSize',16);
title('Normalized Avg Power, \theta = -90','FontSize',16);
axis([0 0.75 0 1.5 4.5 8])
view(-40,20);
grid on
legend(num2str(conts'));

%movieims(3) = getframe(gcf);
% 
% %edited by sam
print(gcf,'-depsc','avgpow_90n');
print(gcf,'-djpeg100','avgpow_90n');
% %eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);

figure(4);hold on;
for i = 1:length(conts)
    h = patch(isosurface(X,Y,vr45n,avgPow45n,conts(i)));
    eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
    set(h,'FaceAlpha',0.5)
end
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('V_{r}','FontSize',16);
title('Normalized Avg Power, \theta = -45','FontSize',16);
axis([0 0.75 0 1.5 4.5 8])
view(-40,20);
grid on
legend(num2str(conts'));

%movieims(4) = getframe(gcf);
% 
% %edited by sam
print(gcf,'-depsc','avgpow_45n');
print(gcf,'-djpeg100','avgpow_45n');
% %eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);

figure(5);hold on;
for i = 1:length(conts)
    h = patch(isosurface(X,Y,vr0,avgPow0,conts(i)));
    eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
    set(h,'FaceAlpha',0.5)
end
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('V_{r}','FontSize',16);
title('Normalized Avg Power, \theta = 0','FontSize',16);
axis([0 0.75 0 1.5 4.5 8])
view(-40,20);
grid on
legend(num2str(conts'));

%movieims(5) = getframe(gcf);
% 
% %edited by sam
print(gcf,'-depsc','avgpow_0');
print(gcf,'-djpeg100','avgpow_0');
% %eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);

figure(6);hold on;
for i = 1:length(conts)
    h = patch(isosurface(X,Y,vr45,avgPow45,conts(i)));
    eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
    set(h,'FaceAlpha',0.5)
end
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('V_{r}','FontSize',16);
title('Normalized Avg Power, \theta = 45','FontSize',16);
axis([0 0.75 0 1.5 4.5 8])
view(-40,20);
grid on
legend(num2str(conts'));

%movieims(6) = getframe(gcf);
% 
% %edited by sam
print(gcf,'-depsc','avgpow_45');
print(gcf,'-djpeg100','avgpow_45');
% %eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);

figure(7);hold on;
for i = 1:length(conts)
    h = patch(isosurface(X,Y,vr90,avgPow90,conts(i)));
    eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
    set(h,'FaceAlpha',0.5)
end
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('V_{r}','FontSize',16);
title('Normalized Avg Power, \theta = 90','FontSize',16);
axis([0 0.75 0 1.5 4.5 8])
view(-40,20);
grid on
legend(num2str(conts'));

%movieims(7) = getframe(gcf);
% 
% %edited by sam
print(gcf,'-depsc','avgpow_90');
print(gcf,'-djpeg100','avgpow_90');
% %eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);

figure(8);hold on;
for i = 1:length(conts)
    h = patch(isosurface(X,Y,vr135,avgPow135,conts(i)));
    eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
    set(h,'FaceAlpha',0.5)
end
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('V_{r}','FontSize',16);
title('Normalized Avg Power, \theta = 135','FontSize',16);
axis([0 0.75 0 1.5 4.5 8])
view(-40,20);
grid on
legend(num2str(conts'));

%movieims(8) = getframe(gcf);
% 
% %edited by sam
print(gcf,'-depsc','avgpow_135');
print(gcf,'-djpeg100','avgpow_135');
% %eval(['print(gcf,''-djpeg100'',''CL3_vr' vr_set ''');']);

%  movname = 'avgpow_phase';
%  mov= avifile(movname,'compression','Cinepak','FPS',3);
%  mov = addframe(mov,movieims);
%  mov= close(mov);
 