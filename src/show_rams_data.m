%clear all;

load rams_data.txt ;
fr = rams_data(:,1);
ad = rams_data(:,2);
cm = rams_data(:,3);
clv = rams_data(:,4);

%cm = -cla / pi^3 ./ad / 2 ./ fr.^2 ;
cla = -cm * pi^3 .* ad * 2 .* fr.^2 ;

vr = 1./fr ;

dvr = .05 ; % for interpolation
dad = .02 ;
dfr = .005 ;

vrax1 = [1.8 8.7] ; % cla all axes
adax1 = [.08 1.04] ;

vrax2 = [1.8 7.6] ; % clv all axes
adax2 = [.08 1.04] ;

vrax3 = [3.5 7.9] ; % cla mit axes
adax3 = [.22 .94] ; 

vrax4 = [3.5 7.4] ; % clv mit axes
adax4 = [.23 .94] ;

vrax5 = [1.8 8.6] ; % cm all axes
adax5 = [.08 1.04] ;

vrax6 = [3.5 7.25] ; % cm mit axes
adax6 = [.22 .94] ;

frax = [min(fr) max(fr)];

dxfig = .4 ;
dyfig = .4 ;
xllfig = .2 ;
yllfig=.2 ;

% figure(1);clf;hold off;
% %subplot('Position',[xllfig yllfig dxfig dyfig]);
% [xif,yif] = meshgrid([min(fr):dfr:max(fr)],[min(ad):dad:max(ad)]);
% zi = griddata(fr,ad,cla,xif,yif);
% [cs,h] = contour(xif,yif,zi,[-10:.5:10],'k') ;
% clabel(cs,h,[-10:1:10]);
% axis('tight');
% xlabel('f_r');
% ylabel('a/d');
% title('C_{la}');
% %print -djpeg rams_cla_vs_fr.jpg;
% %print -deps rams_cla_vs_fr.eps;

% figure(2);clf;hold off;
% %subplot('Position',[xllfig yllfig dxfig dyfig]);
% [xiv,yiv] = meshgrid([min(vr):dvr:max(vr)],[min(ad):dad:max(ad)]);
% zi = griddata(vr,ad,cla,xiv,yiv);
% contour(xiv,yiv,zi) ;
% [cs,h] = contour(xiv,yiv,zi,[-10:.5:10],'k') ;
% clabel(cs,h,[-10:1:10]);
% axis([vrax1 adax1]);
% xlabel('V_r');
% ylabel('a/d');
% title('C_{la}');
% %print -djpeg rams_cla_vs_vr_all_axes.jpg
% %print -deps rams_cla_vs_vr_all_axes.eps

% figure(3);clf;hold off;
% %subplot('Position',[xllfig yllfig dxfig dyfig]);
% [xif,yif] = meshgrid([min(fr):dfr:max(fr)],[min(ad):dad:max(ad)]);
% zi = griddata(fr,ad,clv,xif,yif);
% [cs,h] = contour(xif,yif,zi,[-10:.2:10],'k') ;
% clabel(cs,h,[-10:.4:10]);
% axis('tight');
% xlabel('f_r');
% ylabel('a/d');
% title('C_{lv}');
% %print -djpeg rams_clv_vs_fr.jpg;
% %print -deps rams_clv_vs_fr.eps;

% figure(4);clf;hold off;
% %subplot('Position',[xllfig yllfig dxfig dyfig]);
% [xiv,yiv] = meshgrid([min(vr):dvr:max(vr)],[min(ad):dad:max(ad)]);
% zi = griddata(vr,ad,clv,xiv,yiv);
% contour(xiv,yiv,zi) ;
% [cs,h] = contour(xiv,yiv,zi,[-10:.2:10],'k') ;
% clabel(cs,h,[-10:.4:10]);
% axis([vrax2 adax2]);
% xlabel('V_r');
% ylabel('a/d');
% title('C_{lv}');
% %print -djpeg rams_clv_vs_vr_all_axes.jpg
% %print -deps rams_clv_vs_vr_all_axes.eps

% figure(5);clf;hold off;
% %subplot('Position',[xllfig yllfig dxfig dyfig]);
% [xiv,yiv] = meshgrid([min(vr):dvr:max(vr)],[min(ad):dad:max(ad)]);
% zi = griddata(vr,ad,cla,xiv,yiv);
% contour(xiv,yiv,zi) ;
% [cs,h] = contour(xiv,yiv,zi,[-10:.25:10],'k') ;
% clabel(cs,h,[-10:.5:10]);
% axis([vrax3 adax3]);
% xlabel('V_r');
% ylabel('a/d');
% title('C_{la}');
% %print -djpeg rams_cla_vs_vr_mit_axes.jpg
% %print -deps rams_cla_vs_vr_mit_axes.eps

figure(6);clf;hold off;
%subplot('Position',[xllfig yllfig dxfig dyfig]);
[xiv,yiv] = meshgrid([min(vr):dvr:max(vr)],[min(ad):dad:max(ad)]);
zi = griddata(vr,ad,clv,xiv,yiv);
contour(xiv,yiv,zi) ;
[cs,h] = contour(xiv,yiv,zi,[-10:.5:10],'k') ;
clabel(cs,h,[-10:.5:10]);
%axis([vrax4 adax4]);
xlabel('V_{r}','FontSize',16);
ylabel('A_{y}/D','FontSize',16);
title('C_{Lv} from Gopalkrishnan');
axis([5 8 0.25 1.5])
set(gca,'FontSize',14);
print(gcf,'-djpeg','CLv_gop');
%print -deps rams_clv_vs_vr_mit_axes.eps

% figure(7);clf;hold off;
% %subplot('Position',[xllfig yllfig dxfig dyfig]);
% [xiv,yiv] = meshgrid([min(vr):dvr:max(vr)],[min(ad):dad:max(ad)]);
% zi = griddata(vr,ad,cm,xiv,yiv);
% contour(xiv,yiv,zi) ;
% [cs,h] = contour(xiv,yiv,zi,[-10:.25:10],'k') ;
% clabel(cs,h,[-10:.5:10]);
% axis([vrax5 adax5]);
% xlabel('V_r');
% ylabel('a/d');
% title('C_m');
% %print -djpeg rams_cm_vs_vr_all_axes.jpg
% %print -deps rams_cm_vs_vr_all_axes.eps

figure(8);clf;hold off;
%subplot('Position',[xllfig yllfig dxfig dyfig]);
[xiv,yiv] = meshgrid([min(vr):dvr:max(vr)],[min(ad):dad:max(ad)]);
zi = griddata(vr,ad,cm,xiv,yiv);
contour(xiv,yiv,zi) ;
[cs,h] = contour(xiv,yiv,zi,[-10:0.5:10],'k') ;
clabel(cs,h,[-10:.5:10]);
%axis([vrax6 adax6]);
xlabel('V_{r}','FontSize',16);
ylabel('A_{y}/D','FontSize',16);
title('C_{my} from Gopalkrishnan');
axis([5 8 0.25 1.5])
set(gca,'FontSize',14);
print(gcf,'-djpeg','Cmy_gop');
%print -deps rams_cm_vs_vr_mit_axes.eps
