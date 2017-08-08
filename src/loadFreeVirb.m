load('new1p0.mat')
Cmy_all(1,:) = Cmy';
Cmx_all(1,:) = Cmx';
Xad_all(1,:) = Xad;
Yad_all(1,:) = Yad;
Vrn_all(1,:) = Vrn;
fyn_all(1,:) = fy;
fxm_all(1,:) = fx;
fy_all(1,:)  = ypeakfreq;
fx_all(1,:)  = xpeakfreq;
Vr_all(1,:)  = Vrn.*fy./(ypeakfreq./2./pi);

load('new1p22.mat')
Cmy_all(2,:) = Cmy';
Cmx_all(2,:) = Cmx';
Xad_all(2,:) = Xad;
Yad_all(2,:) = Yad;
Vrn_all(2,:) = Vrn;
fyn_all(2,:) = fy;
fxn_all(2,:) = fx;
fy_all(2,:)  = ypeakfreq;
fx_all(2,:)  = xpeakfreq;
Vr_all(2,:)  = Vrn.*fy./(ypeakfreq./2./pi);

load('new1p37.mat')
Cmy_all(3,:) = Cmy';
Cmx_all(3,:) = Cmx';
Vrn_all(3,:) = Vrn;
Xad_all(3,:) = Xad;
Yad_all(3,:) = Yad;
fyn_all(3,:) = fy;
fxn_all(3,:) = fx;
fy_all(3,:)  = ypeakfreq;
fx_all(3,:)  = xpeakfreq;
Vr_all(3,:)  = Vrn.*fy./(ypeakfreq./2./pi);

load('new1p52.mat')
Cmy_all(4,:) = Cmy';
Cmx_all(4,:) = Cmx';
Vrn_all(4,:) = Vrn;
Xad_all(4,:) = Xad;
Yad_all(4,:) = Yad;
fyn_all(4,:) = fy;
fxn_all(4,:) = fx;
fy_all(4,:)  = ypeakfreq;
fx_all(4,:)  = xpeakfreq;
Vr_all(4,:)  = Vrn.*fy./(ypeakfreq./2./pi);

load('new1p67.mat')
Cmy_all(5,:) = Cmy';
Cmx_all(5,:) = Cmx';
Vrn_all(5,:) = Vrn;
Xad_all(5,:) = Xad;
Yad_all(5,:) = Yad;
fyn_all(5,:) = fy;
fxn_all(5,:) = fx;
fy_all(5,:)  = ypeakfreq;
fx_all(5,:)  = xpeakfreq;
Vr_all(5,:)  = Vrn.*fy./(ypeakfreq./2./pi);

load('new1p9.mat')
Cmy_all(6,:) = Cmy';
Cmx_all(6,:) = Cmx';
Vrn_all(6,:) = Vrn;
Xad_all(6,:) = Xad;
Yad_all(6,:) = Yad;
fyn_all(6,:) = fy;
fxn_all(6,:) = fx;
fy_all(6,:)  = ypeakfreq;
fx_all(6,:)  = xpeakfreq;
Vr_all(6,:)  = Vrn.*fy./(ypeakfreq./2./pi);


load('phase1p0.mat', 'phasextoy')
phase_ind = find(phasextoy > 180);
phasextoy(phase_ind) = phasextoy(phase_ind)-360;
phase_all(1,:) = phasextoy;

load('phase1p22.mat', 'phasextoy')
phase_ind = find(phasextoy > 180);
phasextoy(phase_ind) = phasextoy(phase_ind)-360;
phase_all(2,:) = phasextoy;

load('phase1p37.mat', 'phasextoy')
phase_ind = find(phasextoy > 180);
phasextoy(phase_ind) = phasextoy(phase_ind)-360;
phase_all(3,:) = phasextoy;

load('phase1p52.mat', 'phasextoy')
phase_ind = find(phasextoy > 180);
phasextoy(phase_ind) = phasextoy(phase_ind)-360;
phase_all(4,:) = phasextoy;

load('phase1p67.mat', 'phasextoy')
phase_ind = find(phasextoy > 180);
phasextoy(phase_ind) = phasextoy(phase_ind)-360;
phase_all(5,:) = phasextoy;

load('phase1p9.mat', 'phasextoy')
phase_ind = find(phasextoy > 180);
phasextoy(phase_ind) = phasextoy(phase_ind)-360;
phase_all(6,:) = phasextoy;

save freeVibForceCoef.mat Cmy_all Cmx_all Vrn_all Vr_all fy_all fx_all ...
    fyn_all fxn_all Yad_all Xad_all phase_all

fnRatio = [1 1.22 1.37 1.52 1.67 1.9];
plotInd = [1,2,3,4,5,6];
compIndx = [9:25];
xlimRange = [4 10];
figure
for i = 1:6
    subplot(2,3,i)
    plot(Vrn_all(plotInd(i),compIndx), Yad_all(plotInd(i),compIndx),'o-')
    xlim(xlimRange)
    ylim([0 1.6])
    title(num2str(fnRatio(plotInd(i))))
end
suptitle('Yad')

figure
for i = 1:6
    subplot(2,3,i)
    plot(Vrn_all(plotInd(i),compIndx), Xad_all(plotInd(i),compIndx),'o-')
    xlim(xlimRange)
    ylim([0 0.75])
    title(num2str(fnRatio(plotInd(i))))
end
suptitle('Xad')

figure
for i = 1:6
    subplot(2,3,i)
    plot(Vrn_all(plotInd(i),compIndx), Cmy_all(plotInd(i),compIndx),'o-')
    xlim(xlimRange)
    title(num2str(fnRatio(plotInd(i))))
end
suptitle('Cmy')

figure
for i = 1:6
    subplot(2,3,i)
    plot(Vrn_all(plotInd(i),compIndx), -Cmx_all(plotInd(i),compIndx),'o-')
    xlim(xlimRange)
    title(num2str(fnRatio(plotInd(i))))
end
suptitle('Cmx')

figure
for i = 1:6
    subplot(2,3,i)
    plot(Vrn_all(plotInd(i),compIndx), phase_all(plotInd(i),compIndx),'o-')
    xlim(xlimRange)
    title(num2str(fnRatio(plotInd(i))))
end
suptitle('Phase')

figure
for i = 1:6
    subplot(2,3,i)
    plot(Vrn_all(plotInd(i),compIndx),fx_all(plotInd(i),compIndx)./fy_all(plotInd(i),compIndx),'o-')
    xlim(xlimRange)
    title(num2str(fnRatio(plotInd(i))))
end
suptitle('Effective Frequency Ratio')