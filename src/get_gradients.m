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
    
    [ind3,ind4] = find(isnan(CDv));
    CDv(ind3,ind4) = 0;
    
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
            newCDv(j,:,i) = CDv(i,k:k+5);
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
    newtheta(:,:,9) = 180*pi/180;
    newCL3(:,:,9) = newCL3(:,:,1);
    newCL1(:,:,9) = newCL1(:,:,1);
    newCD2(:,:,9) = newCD2(:,:,1);
    newCmy(:,:,9) = newCmy(:,:,1);
    newCmx(:,:,9) = newCmx(:,:,1);
    newCLv(:,:,9) = newCLv(:,:,1);
    newCDv(:,:,9) = newCDv(:,:,1);
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
    bigCDv(:,:,:,p) = newCDv;
    bigCLafft(:,:,:,p) = newCLafft;
end

bigXi = interpn(bigX,ifac,int);
bigYi = interpn(bigY,ifac,int);
bigthetai = interpn(bigtheta,ifac,int);
%bigCL1i = interpn(bigCL1,ifac,int);
%bigCL3i = interpn(bigCL3,ifac,int);
%bigCD2i = interpn(bigCD2,ifac,int);
bigvri = interpn(bigvr,ifac,int);
bigpowi = interpn(bigpow,ifac,int);
bigCmyi = interpn(bigCmy,ifac,int);
bigCmxi = interpn(bigCmx,ifac,int);
bigCLvi = interpn(bigCLv,ifac,int);
bigPhasei = interpn(bigPhase,ifac,int);
bigCLai = interpn(bigCLa,ifac,int);
bigCDvi = interpn(bigCDv,ifac,int);
%bigCLaffti = interpn(bigCLafft,ifac,int);

xspace = 0.0375;
yspace = 0.0625;
thspace = pi/16;
vrspace = 0.125;

[f1,f2,f3,f4] = gradient(bigpowi,0.0375,0.0625,pi/16,0.125);
[dlvdy,dlvdx,dlvdt,dlvdr] = gradient(bigCLvi,yspace,xspace,thspace,vrspace);
[ddvdy,ddvdx,ddvdt,ddvdr] = gradient(bigCDvi,yspace,xspace,thspace,vrspace);
[dmydy,dmydx,dmydt,dmydr] = gradient(bigCmyi,yspace,xspace,thspace,vrspace);
[dmxdy,dmxdx,dmxdt,dmxdr] = gradient(bigCmxi,yspace,xspace,thspace,vrspace);

% Vrinput = [5:0.125:7.5];
Vrinput = [5.5:0.125:7];
%Ayinput = [0.25:0.125:1.5];

for m = 1:length(Vrinput)
    %for n = 1:length(Ayinput)
        %Vr_out(m,n) = Vrinput(m);
        %Ay_out(m,n) = Ayinput(n);
        
%         P = [0.0179 -0.3724 2.5888 -6.2932 2.1594];
%         N = length(P);
%         Arat = 0;
%         for u = 1:N
%             Arat = Arat + P(u)*Vrinput(m).^(N-u);
%         end;
        
        %Arat = 0.0873*Vrinput(m).^3 - 1.8451*Vrinput(m).^2 + 12.7064*Vrinput(m) - 28.2932;
        %Theta = 10^3*(0.0317*Vrinput(m).^3 - 0.5664*Vrinput(m).^2 + 3.3146*Vrinput(m) - 6.3152);
        %Theta = -50.7314*Vrinput(m) + 323.2722;
        P = [0.0547 -0.0010 -0.1649 -0.1109 0.0467 0.3060];
        MU = [6.3141 0.7987];
        N = length(P);
        Ax_in = 0;
        
        if Vrinput(m) < 7.25
            for u = 1:N
                Ax_in = Ax_in + P(u)*((Vrinput(m)-MU(1))./MU(2)).^(N-u);
            end;
        else
            Ax_in = Ax_in;
        end
        
        P2 = [16.2550 -4.4545 -80.6032 13.6259 132.7168 10.6943 -83.9169 3.4497];
        N2 = length(P2);
        Theta = 0;
        for u = 1:N2
            Theta = Theta + P2(u)*((Vrinput(m)-MU(1))./MU(2)).^(N2-u);
        end;
        Theta = Theta*pi/180;
        
        P3 = [0.2088 -0.0465 -0.7944 -0.0068 0.5743 0.9015];
        N3 = length(P3);
        Ay_in = 0;
        for u = 1:N3
            Ay_in = Ay_in + P3(u)*((Vrinput(m)-MU(1))./MU(2)).^(N3-u);
        end;
        
        
        P4 = [0.0317 -0.0891 0.2941 -0.5341 -0.1635];
        N4 = length(P4);
        Cmy_in = 0;
        for u = 1:N4
            Cmy_in = Cmy_in + P4(u)*((Vrinput(m)-MU(1))./MU(2)).^(N4-u);
        end;
        
        P5 = [-0.1217 0.0442 0.4662 -0.1233 -0.4054 0.0839];
        N5 = length(P5);
        CLv_in = 0;
        for u = 1:N5
            CLv_in = CLv_in + P5(u)*((Vrinput(m)-MU(1))./MU(2)).^(N5-u);
        end;
        
        P6 = [-0.0147 0.0405 0.0637 -0.1027 -0.1966 0.2572];
        N6 = length(P6);
        CDv_in = 0;
        for u = 1:N6
            CDv_in = CDv_in + P6(u)*((Vrinput(m)-MU(1))./MU(2)).^(N6-u);
        end;
        
        P7 = [-0.1902 0.1410 0.6892 -0.2864 -1.0265 0.2042];
        N7 = length(P7);
        Cmx_in = 0;
        for u = 1:N7
            Cmx_in = Cmx_in + P7(u)*((Vrinput(m)-MU(1))./MU(2)).^(N7-u);
        end;
        
        d = abs(bigXi - Ax_in) + abs(bigYi - Ay_in) + ...
            abs(bigthetai - Theta) + abs(bigvri - Vrinput(m));

%         if min(min(min(min(d)))) == 100
%             vr_force(pp,r) = 0;
%             Yad_force(pp,r) = 0;
%             Xad_force(pp,r) = 0;
%             Theta_force(pp,r) = 0;
%         else
            index(m) = find(d == min(min(min(min(d)))));
            [i1,i2,i3,i4] = ind2sub(size(bigXi),index(m));
            vr_force(m) = bigvri(i1,i2,i3,i4);
            Yad_force(m) = bigYi(i1,i2,i3,i4);
            Yad_free(m) = Ay_in;
            Xad_force(m) = bigXi(i1,i2,i3,i4);
            Xad_free(m) = Ax_in;
            Theta_free(m) = Theta;
            Theta_force(m) = bigthetai(i1,i2,i3,i4);
            CLv_force(m) = bigCLvi(i1,i2,i3,i4);
            CDv_force(m) = bigCDvi(i1,i2,i3,i4);
            Cmy_force(m) = bigCmyi(i1,i2,i3,i4);
            Cmx_force(m) = bigCmxi(i1,i2,i3,i4);
            Cmy_free(m) = Cmy_in;
            Cmx_free(m) = Cmx_in;
            CLv_free(m) = CLv_in;
            CDv_free(m) = CDv_in;
            Cp_force(m) = bigpowi(i1,i2,i3,i4);
            Dcp_Dax(m) = f1(i1,i2,i3,i4);
            Dcp_Day(m) = f2(i1,i2,i3,i4);
            Dcp_Dth(m) = f3(i1,i2,i3,i4);
            Dcp_Dvr(m) = f4(i1,i2,i3,i4);
            Dclv_Dax(m) = dlvdx(i1,i2,i3,i4);
            Dclv_Day(m) = dlvdy(i1,i2,i3,i4);
            Dclv_Dth(m) = dlvdt(i1,i2,i3,i4);
            Dclv_Dvr(m) = dlvdr(i1,i2,i3,i4);
            Dcdv_Dax(m) = ddvdx(i1,i2,i3,i4);
            Dcdv_Day(m) = ddvdy(i1,i2,i3,i4);
            Dcdv_Dth(m) = ddvdt(i1,i2,i3,i4);
            Dcdv_Dvr(m) = ddvdr(i1,i2,i3,i4);
            Dcmy_Dax(m) = dmydx(i1,i2,i3,i4);
            Dcmy_Day(m) = dmydy(i1,i2,i3,i4);
            Dcmy_Dth(m) = dmydt(i1,i2,i3,i4);
            Dcmy_Dvr(m) = dmydr(i1,i2,i3,i4);
            Dcmx_Dax(m) = dmxdx(i1,i2,i3,i4);
            Dcmx_Day(m) = dmxdy(i1,i2,i3,i4);
            Dcmx_Dth(m) = dmxdt(i1,i2,i3,i4);
            Dcmx_Dvr(m) = dmxdr(i1,i2,i3,i4);
%             Phase_out(m,n) = bigPhasei(i1,i2,i3,i4);
%             CLa_out(m,n) = bigCLai(i1,i2,i3,i4);
            %CLafft_out(m,n) = bigCLaffti(i1,i2,i3,i4);
end

% figure(1);clf;hold on;
% plot(vr_force,Cp_force,'ko');
% plot(vr_force,Yad_force,'rs-');
% plot(vr_force,Xad_force,'b^:');
% plot(Vrinput,Yad_free,'gs-');
% xlabel('Vr');
% ylabel('Y/D,X/D,Cp');
% legend('Cp','Y/D','X/D');
% 
% figure(2);clf;hold on;
% plot(vr_force,Cp_force,'ko-');
% plot(vr_force,Dcp_Dax,'rs-.');
% plot(vr_force,Dcp_Day,'b^:');
% plot(vr_force,Dcp_Dth,'gp--');
% plot(vr_force,Dcp_Dvr,'mv:');
% xlabel('Vr');
% ylabel('Cp, Slopes');
% legend('Cp','Dcp/Dax','Dcp/Day','Dcp/Dth','Dcp/Dvr');

del = zeros(4,length(Dclv_Dax));
del2 = zeros(3,length(Dclv_Dax));

for i = 1:length(Dclv_Dax)
%     A = [Dclv_Dax(i) Dclv_Day(i) Dclv_Dth(i); Dcdv_Dax(i) Dcdv_Day(i) Dcdv_Dth(i);Dcmy_Dax(i) Dcmy_Day(i) Dcmy_Dth(i)];
    A = [Dclv_Dax(i) Dclv_Day(i) Dclv_Dth(i) Dclv_Dvr(i); ...
        Dcdv_Dax(i) Dcdv_Day(i) Dcdv_Dth(i) Dcdv_Dvr(i); ...
        Dcmy_Dax(i) Dcmy_Day(i) Dcmy_Dth(i) Dcmy_Dvr(i); ...
        Dcmx_Dax(i) Dcmx_Day(i) Dcmx_Dth(i) Dcmx_Dvr(i)];
%     B = [-CLv_force(i);-CDv_force(i);Cmy_force(i) - Cmy_free(i)];
    %B = [-CLv_free(i);-CDv_free(i);Cmy_force(i) - Cmy_free(i)];
%     B = [-CLv_force(i)+CLv_free(i);-CDv_force(i)+CDv_free(i);-Cmy_force(i) + Cmy_free(i)];
    B = [-CLv_force(i)+CLv_free(i);-CDv_force(i)+CDv_free(i);-Cmy_force(i) + Cmy_free(i); -Cmx_force(i) + Cmx_free(i)];
    del(:,i) = A^(-1)*B;
%     A2 = [Dcp_Dax(i) Dcp_Day(i) Dcp_Dth(i); Dcmy_Dax(i) Dcmy_Day(i) Dcmy_Dth(i); Dcmx_Dax(i) Dcmx_Day(i) Dcmx_Dth(i)];
%     B2 = [0-Cp_force(i); Cmy_free(i) - Cmy_force(i); Cmx_free(i) - CDv_force(i)];
%     del2(:,i) = A2^(-1)*B2;
%     A2 = [Dclv_Dax(i) Dclv_Day(i); Dcdv_Dax(i) Dcdv_Day(i)];
%     B2 = [CLv_free(i)-CLv_force(i); CDv_free(i) - CDv_force(i)];
%     del2(:,i) = A2^(-1)*B2;


end

% load new1p9.mat
% 
% Vr = Vrn*fy./(ypeakfreq/2/pi);
% 
% figure(1);
% plot(Vr(9:25),CLv2(9:25),'ko-',Vrinput,CLv_force,'rs--');
% xlabel('Vr','FontSize',16);
% ylabel('CLv','FontSize',16);
% legend('Free','Force');
% set(gca,'FontSize',14);
% 
% print(gcf,'-djpeg100','clv_comparison');
% 
% figure(2);
% plot(Vr(9:25),-CDv2(9:25),'ko-',Vrinput,CDv_force,'rs--');
% xlabel('Vr','FontSize',16);
% ylabel('CDv','FontSize',16);
% legend('Free','Force');
% set(gca,'FontSize',14);
% 
% print(gcf,'-djpeg100','cdv_comparison');
% 
% figure(3);
% plot(Vr(9:25),Cmy(9:25),'ko-',Vrinput,Cmy_force,'rs--');
% xlabel('Vr','FontSize',16);
% ylabel('Cmy','FontSize',16);
% legend('Free','Force');
% set(gca,'FontSize',14);
% 
% print(gcf,'-djpeg100','cmy_comparison');
% 
% figure(4);
% plot(Vr(9:24),-Cmx(9:24),'ko-',Vrinput,Cmx_force,'rs--');
% xlabel('Vr','FontSize',16);
% ylabel('Cmx','FontSize',16);
% legend('Free','Force');
% set(gca,'FontSize',14);
% 
% print(gcf,'-djpeg100','cmx_comparison');

figure(1);
plot(Vrinput,Yad_free,'ko-',Vrinput,Yad_force,'rs:',Vrinput - del(4,:),Yad_free+del(1,:),'bp--');
xlabel('Vr','FontSize',16);
ylabel('Y/D','FontSize',16);
set(gca,'FontSize',14);
legend('Free','Forced','Stepped','Location','best')
axis([5 7.5 0 1.5])

% print(gcf,'-djpeg100','stepfigs/yad_step')

figure(2);
plot(Vrinput,Xad_free,'ko-',Vrinput,Xad_force,'rs:',Vrinput - del(4,:),Xad_free+del(2,:),'bp--');
xlabel('Vr','FontSize',16);
ylabel('X/D','FontSize',16);
set(gca,'FontSize',14);
legend('Free','Forced','Stepped','Location','best')
axis([5 7.5 0 1])

% print(gcf,'-djpeg100','stepfigs/xad_step')

figure(3);
plot(Vrinput,Theta_free,'ko-',Vrinput,Theta_force,'rs:',Vrinput - del(4,:),Theta_free+del(3,:),'bp--');
xlabel('Vr','FontSize',16);
ylabel('\theta','FontSize',16);
set(gca,'FontSize',14);
legend('Free','Forced','Stepped','Location','best')
axis([5 7.5 -pi pi])

% print(gcf,'-djpeg100','stepfigs/theta_step')

for m = 1:length(Yad_free)
    
    d = abs(bigXi - Xad_free(m)-del(2,m)) + abs(bigYi - Yad_free(m)-del(1,m)) + ...
            abs(bigthetai - Theta-del(3,m)) + abs(bigvri - Vrinput(m) + del(4,m));
        
%         d = abs(bigXi - Xad_free(m)-del2(2,m)) + abs(bigYi - Yad_free(m)-del2(1,m)) + ...
%             abs(bigthetai - Theta) + abs(bigvri - Vrinput(m));
        
        index(m) = find(d == min(min(min(min(d)))));
            [j1,j2,j3,j4] = ind2sub(size(bigXi),index(m));
        
            Yad_force2(m) = bigYi(j1,j2,j3,j4);
            Xad_force2(m) = bigXi(j1,j2,j3,j4);
            Theta_force2(m) = bigthetai(j1,j2,j3,j4);
            CLv_force2(m) = bigCLvi(j1,j2,j3,j4);
            CDv_force2(m) = bigCDvi(j1,j2,j3,j4);
            Cmy_force2(m) = bigCmyi(j1,j2,j3,j4);
            Cp_force2(m) = bigpowi(j1,j2,j3,j4);
            Cmx_force2(m) = bigCmxi(j1,j2,j3,j4);
            vr_force2(m) = bigvri(j1,j2,j3,j4);
end

figure(4);
plot(Vrinput,Cmy_free,'ko-',Vrinput,Cmy_force,'rs:',vr_force2,Cmy_force2,'bp--');
xlabel('Vr','FontSize',16);
ylabel('Cmy','FontSize',16);
set(gca,'FontSize',14);
legend('Free','Forced','Stepped','Location','best')

% print(gcf,'-djpeg100','stepfigs/cmy_step')

figure(5);
plot(Vrinput,Cmx_free,'ko-',Vrinput,Cmx_force,'rs:',vr_force2,Cmx_force2,'bp--');
xlabel('Vr','FontSize',16);
ylabel('Cmx','FontSize',16);
set(gca,'FontSize',14);
legend('Free','Forced','Stepped','Location','best')

% print(gcf,'-djpeg100','stepfigs/cmx_step')

figure(6);
plot(Vrinput,CLv_free,'ko-',Vrinput,CLv_force,'rs:',vr_force2,CLv_force2,'bp--');
xlabel('Vr','FontSize',16);
ylabel('CLv','FontSize',16);
set(gca,'FontSize',14);
legend('Free','Forced','Stepped','Location','best')

% print(gcf,'-djpeg100','stepfigs/clv_step')

figure(7);
plot(Vrinput,CDv_free,'ko-',Vrinput,CDv_force,'rs:',vr_force2,CDv_force2,'bp--');
xlabel('Vr','FontSize',16);
ylabel('CDv','FontSize',16);
set(gca,'FontSize',14);
legend('Free','Forced','Stepped','Location','best')

% print(gcf,'-djpeg100','stepfigs/cdv_step')

figure(8);
plot(vr_force2,Dclv_Dax,'ko-',vr_force2,Dclv_Day,'rs-',vr_force2,Dclv_Dth,'bp-',vr_force2,Dclv_Dvr,'m^-');
xlabel('Vr','FontSize',16);
ylabel('D(CLv)/D(*)','FontSize',16);
set(gca,'FontSize',14);
legend('D(CLv)/D(Ax*)','D(CLv)/D(Ay*)','D(CLv)/D(\theta)','D(CLv)/D(Vr)','Location','northwest')
axis([5 7.5 -10 10]);

% print(gcf,'-djpeg100','stepfigs/dclv_step')

figure(9);
plot(vr_force2,Dcdv_Dax,'ko-',vr_force2,Dcdv_Day,'rs-',vr_force2,Dcdv_Dth,'bp-',vr_force2,Dcdv_Dvr,'m^-');
xlabel('Vr','FontSize',16);
ylabel('D(CDv)/D(*)','FontSize',16);
set(gca,'FontSize',14);
legend('D(CDv)/D(Ax*)','D(CDv)/D(Ay*)','D(CDv)/D(\theta)','D(CDv)/D(Vr)','Location','northeast')
axis([5 7.5 -4 6]);

% print(gcf,'-djpeg100','stepfigs/dcdv_step')

figure(10);
plot(vr_force2,Dcmy_Dax,'ko-',vr_force2,Dcmy_Day,'rs-',vr_force2,Dcmy_Dth,'bp-',vr_force2,Dcmy_Dvr,'m^-');
xlabel('Vr','FontSize',16);
ylabel('D(Cmy)/D(*)','FontSize',16);
set(gca,'FontSize',14);
legend('D(Cmy)/D(Ax*)','D(Cmy)/D(Ay*)','D(Cmy)/D(\theta)','D(Cmy)/D(Vr)','Location','southwest')
axis([5 7.5 -10 4]);

% print(gcf,'-djpeg100','stepfigs/dcmy_step')

figure(11);
plot(vr_force2,Dcmx_Dax,'ko-',vr_force2,Dcmx_Day,'rs-',vr_force2,Dcmx_Dth,'bp-',vr_force2,Dcmx_Dvr,'m^-');
xlabel('Vr','FontSize',16);
ylabel('D(Cmx)/D(*)','FontSize',16);
set(gca,'FontSize',14);
legend('D(Cmx)/D(Ax*)','D(Cmx)/D(Ay*)','D(Cmx)/D(\theta)','D(Cmx)/D(Vr)','Location','southwest')
%axis([5 7.5 -10 4]);

% print(gcf,'-djpeg100','stepfigs/dcmx_step')

