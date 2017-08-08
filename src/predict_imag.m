% Try to predict motions based on structural properties

clear all

load vr7.mat

ifac = 2;
int = 'cubic';

for i = 1:8
        k = 1;
        for j = 1:6
        newX(j,:,i) = allXad(i,k:k+5);
        newY(j,:,i) = allYad(i,k:k+5);
        newtheta(j,:,i) = alltheta(i,k:k+5);
        newCL3(j,:,i) = CL3(i,k:k+5);%./sqrt(2);
        newCL1(j,:,i) = CL1(i,k:k+5);%./sqrt(2);
        newCLv(j,:,i) = CLv(i,k:k+5);%./sqrt(2);
        newCmy(j,:,i) = Cmy(i,k:k+5);%./sqrt(2);
        newCmy_corr(j,:,i) = Cmy_corr(i,k:k+5);
        newCDv(j,:,i) = CDv(i,k:k+5);%./sqrt(2);
        newCmx(j,:,i) = Cmx(i,k:k+5);%./sqrt(2);
        newCDmean(j,:,i) = CDmean(i,k:k+5);
        newCD2(j,:,i) = CD2(i,k:k+5);%./sqrt(2);
        newCD4(j,:,i) = CD4(i,k:k+5);%./sqrt(2);
        newPsi(j,:,i) = Psi(i,k:k+5);
        avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
        newPhase(j,:,i) = Phase(i,k:k+5);
        newPhi2(j,:,i) = Phi2(i,k:k+5);
        newPsi(j,:,i) = Psi_corr(i,k:k+5);
        newfft_posy(j,:,i) = fft_posy(i,k:k+5);
        newfft_accx(j,:,i) = fft_accx(i,k:k+5);
        newfft_velx(j,:,i) = fft_velx(i,k:k+5);
        newfft_posx(j,:,i) = fft_posx(i,k:k+5);
        newfft_accy(j,:,i) = fft_accy(i,k:k+5);
        newfft_vely(j,:,i) = fft_vely(i,k:k+5);
        newfft_cdf(j,:,i) = fft_cdf(i,k:k+5);
        newfft_cl(j,:,i) = fft_cl(i,k:k+5);
        k = k+6;
        end
end
    
newXi = interp3(newX,ifac,int);
newYi = interp3(newY,ifac,int);
newthetai = interp3(newtheta,ifac,int);
newCD2i = interp3(newCD2,ifac,int);
newCL1i = interp3(newCL1,ifac,int);
newPhasei = interp3(newPhase,ifac,int);
newPhi2i = interp3(newPhi2,ifac,int);
newCL3i = interp3(newCL3,ifac,int);
newPsii = interp3(newPsi,ifac,int);
avgPowi = interp3(avgPow,ifac,int);

newfft_posy_ri = interp3(real(newfft_posy),ifac,int);
newfft_accx_ri = interp3(real(newfft_accx),ifac,int);
newfft_velx_ri = interp3(real(newfft_velx),ifac,int);
newfft_posx_ri = interp3(real(newfft_posx),ifac,int);
newfft_accy_ri = interp3(real(newfft_accy),ifac,int);
newfft_vely_ri = interp3(real(newfft_vely),ifac,int);
newfft_cdf_ri = interp3(real(newfft_cdf),ifac,int);
newfft_cl_ri = interp3(real(newfft_cl),ifac,int);
newfft_posy_ii = interp3(imag(newfft_posy),ifac,int);
newfft_accx_ii = interp3(imag(newfft_accx),ifac,int);
newfft_velx_ii = interp3(imag(newfft_velx),ifac,int);
newfft_posx_ii = interp3(imag(newfft_posx),ifac,int);
newfft_accy_ii = interp3(imag(newfft_accy),ifac,int);
newfft_vely_ii = interp3(imag(newfft_vely),ifac,int);
newfft_cdf_ii = interp3(imag(newfft_cdf),ifac,int);
newfft_cl_ii = interp3(imag(newfft_cl),ifac,int);

U = 0.349;
D = 0.0762;
S = 2;
mx = 45.9;
my = 51.8;
kx = 3235;
ky = 1013;
zetax = 0.0620;
zetay = 0.0250;
fny = 0.704;
fnx = 1.336;
by = zetay*2*my*fny*2*pi;
bx = zetax*2*mx*fnx*2*pi;
by=0;
bx=0;
vr_guess = 7;
omega =  2*pi*U/(vr_guess*D);

%lhsx = newXi.*(-4*mx*omega^2+kx).*(omega/2.*cos(newthetai*pi/180) + sqrt(-1)*omega/2.*sin(newthetai*pi/180))./(1/2*1000*U^2*D*S);
%rhsx = omega/2.*(cos(newPhi2i-newthetai*pi/180).*cos(newthetai).*newCD2i - sin(newPhi2i-newthetai*pi/180).*sin(newthetai).*newCD2i) + ...
%    sqrt(-1)*omega/2.*(sin(newPhi2i-newthetai*pi/180).*cos(newthetai).*newCD2i + cos(newPhi2i-newthetai*pi/180).*sin(newthetai).*newCD2i);
%rhsx = omega/2.*(cos(newPhi2i).*cos(newthetai).*newCD2i - sin(newPhi2i).*sin(newthetai).*newCD2i) + ...
%    sqrt(-1)*omega/2.*(sin(newPhi2i).*cos(newthetai).*newCD2i + cos(newPhi2i).*sin(newthetai).*newCD2i);
%lhsy = newYi.*(-my*omega^2+ky)*omega./(1/2*1000*U^2*D*S);
%rhsy = omega.*(newCL1i.*cos(newPhasei) + newCL3i.*cos(newPsii)/3) + sqrt(-1)*omega.*(newCL1i.*sin(newPhasei) + newCL3i.*sin(newPsii)/3);
%lhsnx = newXi.*bx*omega.*(omega/2.*sqrt(-1).*cos(newthetai*pi/180) + omega/2.*sin(newthetai*pi/180))./(1/2*1000*U^2*D*S);
%lhsny = newYi.*by*omega.*(omega.*sqrt(-1))./(1/2*1000*U^2*D*S);
%xr = (-4*mx*omega^2+kx)*newXi.*cos(newthetai*pi/180+pi/2)/(1/2*1000*U^2*D*S);
%xi = (-4*mx*omega^2+kx)*newXi.*sin(newthetai*pi/180+pi/2)/(1/2*1000*U^2*D*S);
%yr = (-my*omega^2+ky)*newYi/(1/2*1000*U^2*D*S);
%yi = zeros(size(yr));
%cdr = newCD2i.*cos(newPhi2i-newthetai*pi/180+pi/2);
%cdi = newCD2i.*sin(newPhi2i-newthetai*pi/180+pi/2);
%cyr = newCL1i.*cos(newPhasei)+newCL3i.*cos(newPsii);
%cyi = newCL1i.*sin(newPhasei)+newCL3i.*sin(newPsii);
%xr = real(lhsx+lhsnx);
%xi = imag(lhsx+lhsnx);
%yr = real(lhsy+lhsny);
%yi = imag(lhsy+lhsny);
%cdr = real(rhsx);
%cdi = imag(rhsx);
%cyr = real(rhsy);
%cyi = imag(rhsy);

xr = (-mx*4*omega^2*D*newXi.*newfft_accx_ri - bx*omega*D*newXi.*newfft_velx_ri + kx*D*newXi.*newfft_posx_ri)./(1/2*1000*U^2*D*S);
xi = (-mx*4*omega^2*D*newXi.*newfft_accx_ii - bx*omega*D*newXi.*newfft_velx_ii + kx*D*newXi.*newfft_posx_ii)./(1/2*1000*U^2*D*S);
yr = (-my*omega^2*D*newYi.*newfft_accy_ri - by*omega*D*newYi.*newfft_vely_ri + ky*D*newYi.*newfft_posy_ri)./(1/2*1000*U^2*D*S);
yi = (-my*omega^2*D*newYi.*newfft_accy_ii - by*omega*D*newYi.*newfft_vely_ii + ky*D*newYi.*newfft_posy_ii)./(1/2*1000*U^2*D*S);
cdr = newfft_cdf_ri.*newCD2i;
cdi = newfft_cdf_ii.*newCD2i;
cyr = newfft_cl_ri.*newCL1i;
cyi = newfft_cl_ii.*newCL1i;

d = abs(xr-cdr) + abs(xi-cdi) + abs(yr-cyr) + abs(yi-cyi);
%d = abs(yr-cyr) + abs(yi-cyi);
%d1 = abs(xr-cdr);
%d2 = abs(xi-cdi);
%d3 = abs(yr-cyr);
%d4 = abs(yi-cyi);

%num = 100;
index = find(d == min(min(min(d))));
%index = find(d1 < num & d2 < num & d3 < num & d4 < num);
[i1,i2,i3] = ind2sub(size(newXi),index);

Xpred = newXi(i1,i2,i3)
Ypred = newYi(i1,i2,i3)
Thetapred = newthetai(i1,i2,i3)

figure(3);clf;hold on;
h = patch(isosurface(newXi,newYi,newthetai,d,min(min(min(d)))*5));
set(h,'FaceColor','g');
set(h,'FaceAlpha',0.5)
xlabel('X/D','FontSize',16);
ylabel('Y/D','FontSize',16);
zlabel('\theta','FontSize',16);
title('D2','FontSize',16);
view(-15,20);
grid on

% figure(4);clf;hold on;
% h = patch(isosurface(newXi,newYi,newthetai,newPhasei,0));
% set(h,'FaceColor','g');
% set(h,'FaceAlpha',0.5)
% h2 = patch(isosurface(newXi,newYi,newthetai,newPhasei,pi));
% set(h2,'FaceColor','r');
% set(h2,'FaceAlpha',0.5)
% h3 = patch(isosurface(newXi,newYi,newthetai,newPhasei,-pi));
% set(h3,'FaceColor','b');
% set(h3,'FaceAlpha',0.5)
% xlabel('X/D','FontSize',16);
% ylabel('Y/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% title('D2','FontSize',16);
% view(-15,20);
% grid on