% Plot data as ternary diagrams

clear all;

vr_set = input('Reduced velocity set to run (e.g. 5, 5p5,...) :','s');

eval(['load vr' vr_set '.mat']);

j = 1;
for i = 1:8
    
    X(j:j+35) = allXad(i,:);
    Y(j:j+35) = allYad(i,:);
    newtheta(j:j+35) = alltheta(i,:);
    newCL3(j:j+35) = CL3(i,:);
    newCL5(j:j+35) = CL5(i,:);
    newCLv(j:j+35) = CLv(i,:);
    newCmy(j:j+35) = Cmy(i,:);
    newCDv(j:j+35) = CDv(i,:);
    newCmx(j:j+35) = Cmx(i,:);
    newCDmean(j:j+35) = CDmean(i,:);
    newCD2(j:j+35) = CD2(i,:);
    newCD4(j:j+35) = CD4(i,:);
    avgPow(j:j+35) = (Plift(i,:)+Pdrag(i,:))./(0.5*1000*0.23.^3*0.6858*0.0381);
    %newCL3(j:j+35) = CL3(i,:);
    %newYad(j:j+35) = allYad(i,:);
    %newXad(j:j+35) = allXad(i,:);
    %newtheta(j:j+35) = alltheta(i,:);
    j = j+36;
end

figure(1);
[h,hg,htick] = terplot;
h = ternaryc(X,Y,newtheta,newCL3,'o');
hlabel = terlabel('X/D','Y/D','\theta');
set(htick,[-180 180])