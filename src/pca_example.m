% This program is used to explain how PCA works and the idea behind the
% method. This example is given in the ISA book (Chapter 8, simple example)
% The process is composed of a valve a thermometer and a fermenter. The
% valve is used to control the flow rate of cooling water. The thermometer
% is used to measure the temperature of the fermenter.
% The steady state of the system is: valve=10% max. Temperature T=25 D.

% YANG ZHANG         Email: yz@che.utexas.edu
% 2006-April, UTexas-Austin

clear all
% Give valve set points and fermenter temperature.
% Valve=10+rand(1,10)-0.5; % ten value of Valve around mean.
% Temp=25-(Valve-10)+(rand(1,10)-0.5)*0.5; % 10 value of Temp: T=2.5*valve + noise
Temp  = [2.5000e+001  2.5000e+001  2.5000e+001  2.5300e+001  2.5200e+001  2.5200e+001  2.4800e+001  2.4700e+001  2.5100e+001  2.4700e+001];
Valve = [1.0100e+001  9.9000e+000  1.0000e+001  9.8000e+000  9.9000e+000  9.8000e+000  1.0100e+001  1.0300e+001  1.0000e+001  1.0100e+001];

% Define database
x=[Valve;Temp]';

figure(1)
plot(x(:,1),x(:,2),'ko')
xlabel('Valve Position (%)')
ylabel('Fermenter Temperature (C)')

% Calculate the covariance of dataset x (relationship between Valve and
% Temperature)
covariance=cov(x);

% Singular Value Decomposition x=T*E*P
% E is the eignvalue of x; U is orthnormal; V=U'.
%[U E V]=svd(covariance);
%U=U'; V=V';
[U E] = eig(covariance);
E = diag([E(2,2),E(1,1)]);  % arrange the Eignvalue in a decending order
U = [U(:,2),U(:,1)];  % arrange the Eignvector to make it compatable with Eignvalue.

% Mean center the raw data
Valvebar=x(:,1)-mean(x(:,1));
Tempbar=x(:,2)-mean(x(:,2));
xbar=[Valvebar,Tempbar];

% Principal axis rotation of the covariance matrix.
z=U'*xbar';
z=z';

%covz is the covariance matrix of principal components which is equal to E
%matrix. The variance of transformed variable (z1 & z2) will have variance
%0.2283 and 0.0123 respectively.
covz=U'*covariance*U;
% The first column of U factor is - and + with nearly equal value which
% means the first principal component is related to variability which both
% measurements have difference. The second column is both positive which
% means the 2nd PC is concentrate on the varaiability of the common between
% the two.

% scaling of PCS.
for i=1:2
Vs(:,i)=sqrt(E(i,i))*U(:,i);
Ws(:,i)=U(:,i)/sqrt(E(i,i));
end
% after rescale: Vs'*Vs=E; V'*COV*Vs=E^2; Ws'*Ws=inv(E); Ws'*COV*Ws=I
y=Ws'*xbar';
%define T2
T2=diag(y'*y);

% find the control limit of T2
p=2; %PC number
n=10; % sample number
T2Upper=p*(n-1)/(n-p)*finv(0.95,p,n-p);

% draw the T2 control charts
figure(2)
plot(T2,'ko')
hold on
plot(T2Upper*ones(1,13),'k')
xlabel('Sample Number')
ylabel('T2 Score')
hold off

% Draw the control ellipse
% s1^2*s2^2/(s1^2*s2^2-s12^2)*[(x1-x1mean)^2/s1^2+(x2-x2mean)^2/s2^2-2*s12(
% x1-x1mean)(x2-x2mean)/s1^2/s2^2]=T2^2(upper control limit)
% S=covariance;
% a=S(1)*S(4)/(S(1)*S(4)-S(2)^2);
% b=xbar(:,1).^2/S(1)+xbar(:,2).^2/S(4)-2*S(2)*xbar(:,1).*xbar(:,2)./S(1)/S(4);


% error detection
% two new points is added to be detected which are [10.1,25.2];[9.9;24];

% plot two values with the original 10.
Valve2=[10.2,9.9];Temp2=[25.1,25.5];
figure(1)
hold on
plot(Valve2,Temp2,'k*')
hold off

% plot the control charts for valve and temperature separately.
figure(3)
hold on
plot(1:10,Valve,'ko')
plot(11:12,Valve2,'k*')
plot(0:13, mean(Valve)*ones(1,14),'k')
plot(0:13, (mean(Valve)-std(Valve)*2)*ones(1,14),'k')
plot(0:13, (mean(Valve)+std(Valve)*2)*ones(1,14),'k')
xlabel('Sample Number')
ylabel('Valve Position (%)')
hold off

figure(4)
hold on
plot(1:10,Temp,'ko')
plot(11:12,Temp2,'k*')
plot(0:13, mean(Temp)*ones(1,14),'k')
plot(0:13, (mean(Temp)+std(Temp)*2)*ones(1,14),'k')
plot(0:13, (mean(Temp)-std(Temp)*2)*ones(1,14),'k')
xlabel('Sample Number')
ylabel('Temperature (C)')
hold off

% transform new observations to PC axis
y2=Ws'*([Valve2-mean(Valve);Temp2-mean(Temp)]);

% plot the control charts for principal components
figure(5)
hold on
plot(1:10,y(1,:),'ko')
plot(11:12,y2(1,:),'k*')
plot(0:13, mean(y(1,:))*ones(1,14),'k')
plot(0:13, -2*ones(1,14),'k')
plot(0:13, 2*ones(1,14),'k')
xlabel ('Sample Number')
ylabel ('Principal Component 1')
hold off

figure(6)
hold on
plot(1:10,y(2,:),'ko')
plot(11:12,y2(2,:),'k*')
plot(0:13,mean(y(2,:))*ones(1,14),'k')
plot(0:13, -2*ones(1,14),'k')
plot(0:13, 2*ones(1,14),'k')
xlabel ('Sample Number')
ylabel ('Principal Component 2')
hold off

% plot the T2 chart for the two new points
figure(2)
T22=diag(y2'*y2);
hold on
plot([11,12],T22,'k*')
hold off

% draw the control ellipse
ra = sqrt(T2Upper*E(1,1));
rb = sqrt(T2Upper*E(2,2));
V=[sqrt(E(1,1))*U(:,1),sqrt(E(2,2))*-U(:,2)];
k1 = V(2,1)/V(1,1);
k2 = V(2,2)/V(1,2);
ang = atan(k1); %atan(k1);
x0 = mean(x(:,1));
y0 = mean(x(:,2));
% i = 9.5:0.01:10.5;
% s = 9.72:0.01:10.3;
% j = k1*i + 42.6303;
% l = k2*s + 19.3279;
figure (7)
h=ellipse(ra,rb,ang,x0,y0,'k',1000);
hold on
plot(x(:,1),x(:,2),'ko')
plot(Valve2, Temp2,'k*')
% plot(i,j)
% plot(s,l)
xlabel ('Valve Position')
ylabel ('Temperature')
hold off