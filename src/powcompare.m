clear all

names = {'vr5'; 'vr5p5'; 'vr6'; 'vr7'; 'vr8'};
cols = {'m' 'b' 'g' 'y' 'r'};

for m = 1:length(names)
    eval(['load ' char(names(m))]);
    for i = 1:8
        k = 1;
        for j = 1:6
        X(j,:,i) = allXad(i,k:k+5);
        Y(j,:,i) = allYad(i,k:k+5);
        newtheta(j,:,i) = alltheta(i,k:k+5);
        newCL3(j,:,i) = CL3(i,k:k+5);
        newCL5(j,:,i) = CL5(i,k:k+5);
        newCLv(j,:,i) = CLv(i,k:k+5);
        newCmy(j,:,i) = Cmy(i,k:k+5);
        newCDv(j,:,i) = CDv(i,k:k+5);
        newCmx(j,:,i) = Cmx(i,k:k+5);
        newCDmean(j,:,i) = CDmean(i,k:k+5);
        newCD2(j,:,i) = CD2(i,k:k+5);
        newCD4(j,:,i) = CD4(i,k:k+5);
        avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
        k = k+6;
        end
    end

    

    figure(1);hold on;
    h = patch(isosurface(X,Y,newtheta,avgPow,0));
    eval(['set(h,''FaceColor'',''' char(cols(m)) ''');']);
    set(h,'FaceAlpha',0.5)
    xlabel('X/D','FontSize',16);
    ylabel('Y/D','FontSize',16);
    zlabel('\theta','FontSize',16);
    title('Avg Pow','FontSize',16);
    view(-15,20);
    grid on
    pause
end
legend(char(names));