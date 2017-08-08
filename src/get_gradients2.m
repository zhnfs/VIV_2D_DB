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

[dcpdy,dcpdx,dcpdt,dcpdr] = gradient(bigpowi,0.0375,0.0625,pi/16,0.125);
[dlvdy,dlvdx,dlvdt,dlvdr] = gradient(bigCLvi,yspace,xspace,thspace,vrspace);
[ddvdy,ddvdx,ddvdt,ddvdr] = gradient(bigCDvi,yspace,xspace,thspace,vrspace);
[dmydy,dmydx,dmydt,dmydr] = gradient(bigCmyi,yspace,xspace,thspace,vrspace);
[dmxdy,dmxdx,dmxdt,dmxdr] = gradient(bigCmxi,yspace,xspace,thspace,vrspace);


Vrinput = [5:0.125:7.5];
%Ayinput = [0.25:0.125:1.5];
MU = [6.3141 0.7987];

for m = 1:length(Vrinput)
    RMSE = 1;
    del = zeros(3,1);
    count = 1;
    
    P4 = [0.0317 -0.0891 0.2941 -0.5341 -0.1635];
    N4 = length(P4);
    free_cmy = 0;
    for u = 1:N4
        free_cmy = free_cmy + P4(u)*((Vrinput(m)-MU(1))./MU(2)).^(N4-u);
    end;

    P5 = [-0.1217 0.0442 0.4662 -0.1233 -0.4054 0.0839];
    N5 = length(P5);
    free_clv = 0;
    for u = 1:N5
        free_clv = free_clv + P5(u)*((Vrinput(m)-MU(1))./MU(2)).^(N5-u);
    end;

    P6 = [-0.0147 0.0405 0.0637 -0.1027 -0.1966 0.2572];
    N6 = length(P6);
    free_cdv = 0;
    for u = 1:N6
        free_cdv = free_cdv + P6(u)*((Vrinput(m)-MU(1))./MU(2)).^(N6-u);
    end;
    
    P7 = [-0.1902 0.1410 0.6892 -0.2864 -1.0265 0.2042];
    N7 = length(P7);
    free_cmx = 0;
    for u = 1:N7
        free_cmx = free_cmx + P7(u)*((Vrinput(m)-MU(1))./MU(2)).^(N7-u);
    end;
    
    while RMSE > 0.05 & count < 10 
    
        if count == 1
        
            P = [0.0547 -0.0010 -0.1649 -0.1109 0.0467 0.3060];
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

            d = abs(bigXi - Ax_in) + abs(bigYi - Ay_in) + ...
                abs(bigthetai - Theta) + abs(bigvri - Vrinput(m));
        else
            d = abs(bigXi - Xad_in - del(1)) + abs(bigYi - Yad_in - del(2)) + ...
                abs(bigthetai - Theta_in - del(3)) + abs(bigvri - Vrinput(m));
        end


        index = find(d == min(min(min(min(d)))));
        [i1,i2,i3,i4] = ind2sub(size(bigXi),index);
        
        Yad_in = bigYi(i1,i2,i3,i4);
        Xad_in = bigXi(i1,i2,i3,i4);
        Theta_in = bigthetai(i1,i2,i3,i4);
        CLv_in = bigCLvi(i1,i2,i3,i4);
        CDv_in = bigCDvi(i1,i2,i3,i4);
        Cmy_in = bigCmyi(i1,i2,i3,i4);
        Cmx_in = bigCmxi(i1,i2,i3,i4);
        Cp_in = bigpowi(i1,i2,i3,i4);
        dcp_dax = dcpdx(i1,i2,i3,i4);
        dcp_day = dcpdy(i1,i2,i3,i4);
        dcp_dth = dcpdt(i1,i2,i3,i4);
        dcp_dvr = dcpdr(i1,i2,i3,i4);
        dclv_dax = dlvdx(i1,i2,i3,i4);
        dclv_day = dlvdy(i1,i2,i3,i4);
        dclv_dth = dlvdt(i1,i2,i3,i4);
        dclv_dvr = dlvdr(i1,i2,i3,i4);
        dcdv_dax = ddvdx(i1,i2,i3,i4);
        dcdv_day = ddvdy(i1,i2,i3,i4);
        dcdv_dth = ddvdt(i1,i2,i3,i4);
        dcdv_dvr = ddvdr(i1,i2,i3,i4);
        dcmy_dax = dmydx(i1,i2,i3,i4);
        dcmy_day = dmydy(i1,i2,i3,i4);
        dcmy_dth = dmydt(i1,i2,i3,i4);
        dcmy_dvr = dmydr(i1,i2,i3,i4);
        dcmx_dax = dmxdx(i1,i2,i3,i4);
        dcmx_day = dmxdy(i1,i2,i3,i4);
        dcmx_dth = dmxdt(i1,i2,i3,i4);
        dcmx_dvr = dmxdr(i1,i2,i3,i4);
            
%             A = [dclv_dax dclv_day dclv_dth; ...
%                 dcdv_dax dcdv_day dcdv_dth; ...
%                 dcmy_dax dcmy_day dcmy_dth];
%             B = [free_clv - CLv_in; ...
%                 free_cdv - CDv_in; ...
%                 free_cmy - Cmy_in];
            A = [dcp_dax dcp_day dcp_dth; ...
                dcmx_dax dcmx_day dcmx_dth; ...
                dcmy_dax dcmy_day dcmy_dth];
            B = [0 - Cp_in; ...
                free_cmx - Cmx_in; ...
                free_cmy - Cmy_in];
            del = A^(-1)*B*0.05;
            RMSE = sqrt(sum(B.^2)/length(B));
            count = count+1;
        end
    vr_force(m) = bigvri(i1,i2,i3,i4);
    Yad_force(m) = Yad_in;
    Yad_free(m) = Ay_in;
    Xad_force(m) = Xad_in;
    Xad_free(m) = Ax_in;
    Theta_free(m) = Theta;
    Theta_force(m) = Theta_in;
    CLv_force(m) = CLv_in;
    CDv_force(m) = CDv_in;
    Cmy_force(m) = Cmy_in;
    Cmx_force(m) = Cmx_in;
    Cmy_free(m) = free_cmy;
    CLv_free(m) = free_clv;
    CDv_free(m) = free_cdv;
    Cmx_free(m) = free_cmx;
    Cp_force(m) = Cp_in;
    Dcp_Dax(m) = dcp_dax;
    Dcp_Day(m) = dcp_day;
    Dcp_Dth(m) = dcp_dth;
    Dcp_Dvr(m) = dcp_dvr;
    Dclv_Dax(m) = dclv_dax;
    Dclv_Day(m) = dclv_day;
    Dclv_Dth(m) = dclv_dth;
    Dclv_Dvr(m) = dclv_dvr;
    Dcdv_Dax(m) = dcdv_dax;
    Dcdv_Day(m) = dcdv_day;
    Dcdv_Dth(m) = dcdv_dth;
    Dcdv_Dvr(m) = dcdv_dvr;
    Dcmy_Dax(m) = dcmy_dax;
    Dcmy_Day(m) = dcmy_day;
    Dcmy_Dth(m) = dcmy_dth;
    Dcmy_Dvr(m) = dcmy_dvr;
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

% del = zeros(3,length(Dclv_Dax));
% del2 = zeros(2,length(Dclv_Dax));
% 
% for i = 1:length(Dclv_Dax)
%     A = [Dclv_Dax(i) Dclv_Day(i) Dclv_Dth(i); Dcdv_Dax(i) Dcdv_Day(i) Dcdv_Dth(i);Dcmy_Dax(i) Dcmy_Day(i) Dcmy_Dth(i)];
% %     B = [-CLv_force(i);-CDv_force(i);Cmy_force(i) - Cmy_free(i)];
%     %B = [-CLv_free(i);-CDv_free(i);Cmy_force(i) - Cmy_free(i)];
%     B = [CLv_force(i)-CLv_free(i);CDv_force(i)-CDv_free(i);Cmy_force(i) - Cmy_free(i)];
%     del(:,i) = A^(-1)*B;
%     A2 = [Dclv_Dax(i) Dclv_Day(i); Dcdv_Dax(i) Dcdv_Day(i)];
%     B2 = [CLv_force(i)-CLv_free(i);CDv_force(i)-CDv_free(i)];
%     del2(:,i) = A2^(-1)*B2;
% end

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
plot(Vrinput,Yad_free,'ko-',Vrinput,Yad_force,'rs--');
xlabel('Vr');
ylabel('Y/D');

figure(2);
plot(Vrinput,Xad_free,'ko-',Vrinput,Xad_force,'rs--');
xlabel('Vr');
ylabel('X/D');

figure(3);
plot(Vrinput,Theta_free,'ko-',Vrinput,Theta_force,'rs--');
xlabel('Vr');
ylabel('\theta');

figure(4);
plot(Vrinput,Cmy_free,'ko-',Vrinput,Cmy_force,'rs:');
xlabel('Vr');
ylabel('Cmy');

figure(5);
plot(Vrinput,Cmx_free,'ko-',Vrinput,Cmx_force,'rs:');
xlabel('Vr');
ylabel('Cmx');

figure(6);
plot(Vrinput,0,'ko-',Vrinput,Cp_force,'rs:');
xlabel('Vr');
ylabel('Cp');

% for m = 1:length(Yad_free)
%     
%     d = abs(bigXi - Xad_free(m)-del(2,m)) + abs(bigYi - Yad_free(m)-del(1,m)) + ...
%             abs(bigthetai - Theta-del(3,m)) + abs(bigvri - Vrinput(m));
%         
% %         d = abs(bigXi - Xad_free(m)-del2(2,m)) + abs(bigYi - Yad_free(m)-del2(1,m)) + ...
% %             abs(bigthetai - Theta) + abs(bigvri - Vrinput(m));
%         
%         index(m) = find(d == min(min(min(min(d)))));
%             [j1,j2,j3,j4] = ind2sub(size(bigXi),index(m));
%         
%             Yad_force2(m) = bigYi(j1,j2,j3,j4);
%             Xad_force2(m) = bigXi(j1,j2,j3,j4);
%             Theta_force2(m) = bigthetai(j1,j2,j3,j4);
%             CLv_force2(m) = bigCLvi(j1,j2,j3,j4);
%             CDv_force2(m) = bigCDvi(j1,j2,j3,j4);
%             Cmy_force2(m) = bigCmyi(j1,j2,j3,j4);
%             Cp_force2(m) = bigpowi(j1,j2,j3,j4);
% end
% 
% figure(4);
% plot(Vrinput,Cmy_force,'ko-',Vrinput,Cmy_force2,'rs:',Vrinput,Cmy_free,'b^--');
% xlabel('Vr');
% ylabel('Force Coeff');


