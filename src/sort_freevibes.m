% Re-sort free vibration database
% Jason Dahl 6/12/08

clear all;

load new1p9;

Vr = Vrn.*fy./(ypeakfreq/2/pi);

Vrbins = [5 5.5 6 6.5 7 7.5];

Vr5ind = find(Vr < 5.25 & Vr > 4.76);
Vr5p5ind = find(Vr < 5.75 & Vr > 5.26);
Vr6ind = find(Vr < 6.25 & Vr > 5.76);
Vr6p5ind = find(Vr < 6.75 & Vr > 6.26);
Vr7ind = find(Vr < 7.25 & Vr > 6.76);
Vr7p5ind = find(Vr < 7.75 & Vr > 7.26);

figure(1);
plot(CLv2(Vr5ind),Yad(Vr5ind),'ko');
xlabel('CLv');
ylabel('A_{y}/D');
title('Vr = 5');

figure(2);
plot(CLv2(Vr5p5ind),Yad(Vr5p5ind),'ko');
xlabel('CLv');
ylabel('A_{y}/D');
title('Vr = 5.5');

figure(3);
plot(CLv2(Vr6ind),Yad(Vr6ind),'ko');
xlabel('CLv');
ylabel('A_{y}/D');
title('Vr = 6');

figure(4);
plot(CLv2(Vr6p5ind),Yad(Vr6p5ind),'ko');
xlabel('CLv');
ylabel('A_{y}/D');
title('Vr = 6.5');

figure(5);
plot(CLv2(Vr7ind),Yad(Vr7ind),'ko');
xlabel('CLv');
ylabel('A_{y}/D');
title('Vr = 7');

figure(6);
plot(CLv2(Vr7p5ind),Yad(Vr7p5ind),'ko');
xlabel('CLv');
ylabel('A_{y}/D');
title('Vr = 7.5');