% Find optimal 1-dof database from forced vibes based on free vibrations
% Jason Dahl 7/7/08

clear all

vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};

ifac = 2;
int = 'spline';

bigX = zeros(6,6,9,8);
bigY = zeros(6,6,9,8);
bigtheta = zeros(6,6,9,8);
bigCL1 = zeros(6,6,9,8);
bigCL3 = zeros(6,6,9,8);
bigCD2 = zeros(6,6,9,8);
bigCmy = zeros(6,6,9,8);
bigCmx = zeros(6,6,9,8);
bigvr = zeros(6,6,9,8);
bigpow = zeros(6,6,9,8);
bigrat = zeros(6,6,9,8);
bigPhase = zeros(6,6,9,8);


for p = 1:length(vr_set)
    eval(['load vr' char(vr_set(p)) '.mat']);
    
    [ind1,ind2] = find(isnan(Cmx));
    Cmx(ind1,ind2) = 1;
    
    for i = 1:8
            k = 1;
            for j = 1:6
            newX(j,:,i) = allXad(i,k:k+5);
            newY(j,:,i) = allYad(i,k:k+5);
            newtheta(j,:,i) = alltheta(i,k:k+5)*pi/180;
            newCL3(j,:,i) = CL3(i,k:k+5);
            newCL1(j,:,i) = CL1(i,k:k+5);
            newCD2(j,:,i) = CD2(i,k:k+5);
            newCmy(j,:,i) = Cmy(i,k:k+5);
            newCmx(j,:,i) = Cmx(i,k:k+5);
            newCLv(j,:,i) = CLv(i,k:k+5);
            newCLa(j,:,i) = CLa(i,k:k+5);
            newCLafft(j,:,i) = CLafft(i,k:k+5);
            newPhase(j,:,i) = Phase(i,k:k+5);
            vr(j,:,i) = Vr(i,k:k+5);
            
            if p == 1 | p == 2 | p == 3
                avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.18.^3*0.6858*0.0381);
            else
                avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
            end
            k = k+6;
            end
    end
    newX(:,:,9) = newX(:,:,1);
    newY(:,:,9) = newY(:,:,1);
    newtheta(:,:,9) = newtheta(:,:,1);
    newCL3(:,:,9) = newCL3(:,:,1);
    newCL1(:,:,9) = newCL1(:,:,1);
    newCD2(:,:,9) = newCD2(:,:,1);
    newCmy(:,:,9) = newCmy(:,:,1);
    newCmx(:,:,9) = newCmx(:,:,1);
    newCLv(:,:,9) = newCLv(:,:,1);
    newPhase(:,:,9) = newPhase(:,:,1);
    newCLa(:,:,9) = newCLa(:,:,1);
    newCLafft(:,:,9) = newCLafft(:,:,1);
    vr(:,:,9) = vr(:,:,1);
    avgPow(:,:,9) = avgPow(:,:,1);
    
    bigX(:,:,:,p) = newX;
    bigY(:,:,:,p) = newY;
    bigtheta(:,:,:,p) = newtheta;
    bigCL1(:,:,:,p) = newCL1;
    bigCL3(:,:,:,p) = newCL3;
    bigCD2(:,:,:,p) = newCD2;
    bigvr(:,:,:,p) = vr;
    bigpow(:,:,:,p) = avgPow;
    bigCmy(:,:,:,p) = newCmy;
    bigCmx(:,:,:,p) = newCmx;
    bigCLv(:,:,:,p) = newCLv;
    bigPhase(:,:,:,p) = newPhase;
    bigCLa(:,:,:,p) = newCLa;
    bigCLafft(:,:,:,p) = newCLafft;
end

bigXi = interpn(bigX,ifac,int);
bigYi = interpn(bigY,ifac,int);
bigthetai = interpn(bigtheta,ifac,int);
%bigCL1i = interpn(bigCL1,ifac,int);
%bigCL3i = interpn(bigCL3,ifac,int);
%bigCD2i = interpn(bigCD2,ifac,int);
bigvri = interpn(bigvr,ifac,int);
%bigpowi = interpn(bigpow,ifac,int);
bigCmyi = interpn(bigCmy,ifac,int);
%bigCmxi = interpn(bigCmx,ifac,int);
bigCLvi = interpn(bigCLv,ifac,int);
bigPhasei = interpn(bigPhase,ifac,int);
bigCLai = interpn(bigCLa,ifac,int);
%bigCLaffti = interpn(bigCLafft,ifac,int);

Vrinput = [5:0.125:8];
Ayinput = [0.25:0.125:1.5];

for m = 1:length(Vrinput)
    for n = 1:length(Ayinput)
        Vr_out(m,n) = Vrinput(m);
        Ay_out(m,n) = Ayinput(n);
        
%         P = [0.0179 -0.3724 2.5888 -6.2932 2.1594];
%         N = length(P);
%         Arat = 0;
%         for u = 1:N
%             Arat = Arat + P(u)*Vrinput(m).^(N-u);
%         end;
        
        Arat = 0.0873*Vrinput(m).^3 - 1.8451*Vrinput(m).^2 + 12.7064*Vrinput(m) - 28.2932;
        %Theta = 10^3*(0.0317*Vrinput(m).^3 - 0.5664*Vrinput(m).^2 + 3.3146*Vrinput(m) - 6.3152);
        %Theta = -50.7314*Vrinput(m) + 323.2722;
        P2 = [78.3995 -3482.334 66040.8319 -693152.132 4348445.148 -16305237.291 33836700.938 -29979191.49];
        N2 = length(P2);
        Theta = 0;
        for u = 1:N2
            Theta = Theta + P2(u)*Vrinput(m).^(N2-u);
        end;
        
        %Ax_in = Arat.*Ayinput(n);
%         P = [0.1684 -5.3181 66.8639 -418.4617 1304.0536 -1619.0997];
%         N = length(P);
%         Ax_in = 0;
%         if Vrinput(m) < 7.25
%             for u = 1:N
%                 Ax_in = Ax_in + P(u)*Vrinput(m).^(N-u);
%             end;
%         else
%             Ax_in = Ax_in;
%         end
        
        
        Ax_in = 0.25;
        
%         dvr = abs(bigvri-Vrinput(m));
%         indices = find(dvr < 1*0.5/(ifac^2));
%         [p1,p2,p3,p4] = ind2sub(size(bigXi),indices);
%         d = 100*ones(size(bigXi));
%         
%         for z = 1:length(p1)
%             d(p1(z),p2(z),p3(z),p4(z)) = abs((bigXi(p1(z),p2(z),p3(z),p4(z)) - ...
%                 Ax_in)) + abs((bigYi(p1(z),p2(z),p3(z),p4(z)) - ...
%                 Ayinput(n))) + abs(bigthetai(p1(z),p2(z),p3(z),p4(z)) - Theta);
%         end
        
        d = abs(bigXi - Ax_in) + abs(bigYi -Ayinput(n)) + ...
            abs(bigthetai - Theta) + abs(bigvri - Vrinput(m));

%         if min(min(min(min(d)))) == 100
%             vr_force(pp,r) = 0;
%             Yad_force(pp,r) = 0;
%             Xad_force(pp,r) = 0;
%             Theta_force(pp,r) = 0;
%         else
            index(m,n) = find(d == min(min(min(min(d)))));
            [i1,i2,i3,i4] = ind2sub(size(bigXi),index(m,n));
            vr_force(m,n) = bigvri(i1,i2,i3,i4);
            Yad_force(m,n) = bigYi(i1,i2,i3,i4);
            Xad_force(m,n) = bigXi(i1,i2,i3,i4);
            Theta_force(m,n) = bigthetai(i1,i2,i3,i4);
            CLv_force(m,n) = bigCLvi(i1,i2,i3,i4);
            Cmy_force(m,n) = bigCmyi(i1,i2,i3,i4);
            %Cp_force(m,n) = bigpowi(i1,i2,i3,i4);
            Phase_out(m,n) = bigPhasei(i1,i2,i3,i4);
            CLa_out(m,n) = bigCLai(i1,i2,i3,i4);
            %CLafft_out(m,n) = bigCLaffti(i1,i2,i3,i4);
    end
end

%show_rams_data;

figure(6);hold on;
[C1,H1] = contour(Vr_out,Ay_out,CLv_force,[-10:0.5:10]);
xlabel('V_{r}','FontSize',16);
ylabel('A_{y}/D','FontSize',16);
title('C_{Lv} from Dahl');
axis([5 8 0.25 1.5])
clabel(C1,H1);
set(gca,'FontSize',14);
print(gcf,'-djpeg100','CLv_dahl2');

figure(8);hold on;
[C2,H2] = contour(Vr_out,Ay_out,Cmy_force,[-10:0.5:10]);
xlabel('V_{r}','FontSize',16);
ylabel('A_{y}/D','FontSize',16);
title('C_{my} from Dahl');
axis([5 8 0.25 1.5])
clabel(C2,H2);
set(gca,'FontSize',14);
print(gcf,'-djpeg100','Cmy_dahl2');

% figure(3);
% [C3,H3] = contour(Vr_out,Ay_out,Cp_force);
% xlabel('V_{r}');
% ylabel('A_{y}/D');
% clabel(C3,H3);

% figure(9);hold on;
% [C2,H2] = contour(Vr_out,Ay_out,Phase_out);
% xlabel('V_{r}','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% title('C_{my} from Dahl');
% axis([5 8 0.25 1.5])
% clabel(C2,H2);
% set(gca,'FontSize',14);
% 
% phi = atan(-CLv_force./CLa_out);
% 
% figure(10);hold on;
% [C2,H2] = contour(Vr_out,Ay_out,phi);
% xlabel('V_{r}','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% title('C_{my} from Dahl');
% axis([5 8 0.25 1.5])
% clabel(C2,H2);
% set(gca,'FontSize',14);

% 
% figure(10);hold on;
% [C2,H2] = contour(Vr_out,Ay_out,CLafft_out);
% xlabel('V_{r}','FontSize',16);
% ylabel('A_{y}/D','FontSize',16);
% title('C_{my} from Dahl');
% axis([5 8 0.25 1.5])
% clabel(C2,H2);
% set(gca,'FontSize',14);