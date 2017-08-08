% Make figures for 'optimal' database

clear all;

load new1p9.mat

Vr = Vrn*fy./(ypeakfreq/2/pi);

figure(1);
plot(Vr(4:27),Cmy(4:27),'ko:');
xlabel('V_{r}', 'FontSize',16);
ylabel('C_{my}', 'FontSize',16);
set(gca,'FontSize',14);

%print(gcf,'-djpeg100','Cmy_free');

figure(2);
plot(Vr(4:27),Yad(4:27),'ko:');
xlabel('V_{r}', 'FontSize',16);
ylabel('A_{y}/D', 'FontSize',16);
set(gca,'FontSize',14);

%print(gcf,'-djpeg100','Yad_free');

figure(3);
plot(Vr(4:27),CLv2(4:27),'ko:');
xlabel('V_{r}', 'FontSize',16);
ylabel('C_{Lv}', 'FontSize',16);
set(gca,'FontSize',14);

%print(gcf,'-djpeg100','Clv_free');

Arat = Xad./Yad;

[P,S,MU] = polyfit(Vr(8:27),Xad(8:27),5);

N = length(P);
Pol = zeros(size(Vr(8:27)));

for i = 1:N
    Pol = Pol + P(i)*((Vr(8:27)-MU(1))./MU(2)).^(N-i);
end

%Pol = P(1)*Vr(4:27).^3 + P(2)*Vr(4:27).^2 + P(3)*Vr(4:27) + P(4); 

figure(4);
plot(Vr(8:27),Xad(8:27),'ko',Vr(8:27),Pol,'r--');
xlabel('V_{r}', 'FontSize',16);
ylabel('A_{x}*/A_{y}*', 'FontSize',16);
set(gca,'FontSize',14);

%print(gcf,'-djpeg100','Arat_free');

load phase1p9

phase_ind = find(phasextoy > 180);
phasextoy(phase_ind) = phasextoy(phase_ind)-360;

[P2,S2,MU2] = polyfit(Vr(8:27),phasextoy(8:27)',7);
%[P2,S2] = polyfit(Vr(8:21),phasextoy(8:21)',1);
%[P2,S2] = polyfit([Vr(8:12); Vr(14:21); Vr(26:27)],[phasextoy(8:12)'; phasextoy(14:21)'; phasextoy(26:27)'],3);
%Pol2 = P2(1)*Vr(8:27).^4 + P2(2)*Vr(8:27).^3 + P2(3)*Vr(8:27).^2 + P2(4)*Vr(8:27) + P2(5); 
%Pol2 = P2(1)*Vr(8:21) + P2(2);
%Pol2 = P2(1)*[Vr(8:12); Vr(14:21); Vr(26:27)].^3 + P2(2)*[Vr(8:12); Vr(14:21); Vr(26:27)].^2 + P2(3)*[Vr(8:12); Vr(14:21); Vr(26:27)] + P2(4);

N2 = length(P2);
Pol2 = zeros(size(Vr(8:27)));

for i = 1:N2
    Pol2 = Pol2 + P2(i)*((Vr(8:27)-MU2(1))./MU2(2)).^(N2-i);
end

figure(5);
plot(Vr(8:27),phasextoy(8:27),'ko',Vr(8:27),Pol2,'r--');
xlabel('V_{r}', 'FontSize',16);
ylabel('\theta', 'FontSize',16);
set(gca,'FontSize',14);

%print(gcf,'-djpeg100','Theta_free');


[P3,S3,MU3] = polyfit(Vr(8:27),Yad(8:27),5);
N3 = length(P3);
Pol3 = zeros(size(Vr(8:27)));

for i = 1:N3
    Pol3 = Pol3 + P3(i)*((Vr(8:27)-MU3(1))./MU3(2)).^(N3-i);
end


Vrinput = [5:0.125:8];

P4 = [0.6424 -20.3963 257.4523 -1615.0123 5035.0178 -6240.6325];
N4 = length(P4);
Ay_in = 0;
for u = 1:N4
   Ay_in = Ay_in + P4(u)*Vrinput.^(N3-u);
end;

figure(6);
plot(Vr(8:27),Yad(8:27),'ko',Vr(8:27),Pol3,'r--',Vrinput,Ay_in,'bs:');
xlabel('V_{r}', 'FontSize',16);
ylabel('A_{y}/D', 'FontSize',16);
set(gca,'FontSize',14);

[P5,S5,MU5] = polyfit(Vr(8:27),Cmy(8:27),4);

N5 = length(P5);
Pol5 = zeros(size(Vr(8:27)));

for i = 1:N5
    Pol5 = Pol5 + P5(i)*((Vr(8:27)-MU5(1))./MU5(2)).^(N5-i);
end

figure(7);
plot(Vr(8:27),Cmy(8:27),'ko',Vr(8:27),Pol5,'r--');
xlabel('V_{r}', 'FontSize',16);
ylabel('Cmy', 'FontSize',16);
set(gca,'FontSize',14);

[P6,S6,MU6] = polyfit(Vr(8:27),CLv2(8:27),5);

N6 = length(P6);
Pol6 = zeros(size(Vr(8:27)));

for i = 1:N6
    Pol6 = Pol6 + P6(i)*((Vr(8:27)-MU6(1))./MU6(2)).^(N6-i);
end

figure(8);
plot(Vr(8:27),CLv2(8:27),'ko',Vr(8:27),Pol6,'r--');
xlabel('V_{r}', 'FontSize',16);
ylabel('CLv', 'FontSize',16);
set(gca,'FontSize',14);

[P7,S7,MU7] = polyfit(Vr(8:27),-CDv2(8:27),5);

N7 = length(P7);
Pol7 = zeros(size(Vr(8:27)));

for i = 1:N7
    Pol7 = Pol7 + P7(i)*((Vr(8:27)-MU7(1))./MU7(2)).^(N7-i);
end

figure(9);
plot(Vr(8:27),-CDv2(8:27),'ko',Vr(8:27),Pol7,'r--');
xlabel('V_{r}', 'FontSize',16);
ylabel('CDv', 'FontSize',16);
set(gca,'FontSize',14);

[P8,S8,MU8] = polyfit(Vr(8:24),-Cmx(8:24),5);

N8 = length(P8);
Pol8 = zeros(size(Vr(8:24)));

for i = 1:N8
    Pol8 = Pol8 + P8(i)*((Vr(8:24)-MU8(1))./MU8(2)).^(N8-i);
end

figure(10);
plot(Vr(8:24),-Cmx(8:24),'ko',Vr(8:24),Pol8,'r--');
xlabel('V_{r}', 'FontSize',16);
ylabel('Cmx', 'FontSize',16);
set(gca,'FontSize',14);