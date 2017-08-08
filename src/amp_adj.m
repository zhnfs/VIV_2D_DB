% re-calculate amplitudes from free vibration experiments

clear all;

newmat = {'new1p0' 'new1p22' 'new1p37' 'new1p52' 'new1p67' 'new1p9'};
ampfiles = {'ampfix1p0' 'ampfix1p22' 'ampfix1p37' 'ampfix1p52' 'ampfix1p67' 'ampfix1p9'};

for p = 1:6
    eval(['load ' char(newmat(p)) '.mat']);
    
    for j = 1:37
        if j == 1
            CLcalc = CL2;
            CDcalc = CD2;
            Ycalc = Y./0.0762;
            Xcalc = X./0.0762;
            Tcalc = t;
        else
            eval(['CLcalc = CL2_' num2str(points(j)) ';']);
            eval(['CDcalc = CD2_' num2str(points(j)) ';']);
            eval(['Ycalc = Y_' num2str(points(j)) './0.0762;']);
            eval(['Xcalc = X_' num2str(points(j)) './0.0762;']);
            eval(['Tcalc = t_' num2str(points(j)) ';']);
        end

        [ymax,ymin] = peakdet(Ycalc-mean(Ycalc),0.05);
        [xmax,xmin] = peakdet(Xcalc-mean(Xcalc),0.05);
        
        if isempty(ymax) == 1 %| isempty(ymin) == 1
            Yad_ten(j) = 0;
        elseif isempty(xmax) == 1 %| isempty(xmin) == 1
            Xad_ten(j) = 0;
        else
        maxy = sort(ymax(:,2),'descend');
        num1 = round(1*length(ymax(:,2)));
        Yad_fix(j) = mean(maxy(1:num1));

        maxx = sort(xmax(:,2),'descend');
        num2 = round(1*length(xmax(:,2)));
        Xad_fix(j) = mean(maxx(1:num2));
        end
        
    end
    clear X_* Y_* C* ps* F*
    eval(['save ' char(ampfiles(p)) '.mat']);
    
    figure(p)
    plot(Vrn,Yad,'ko',Vrn,Yad_fix,'ro',Vrn,Xad,'ks',Vrn,Xad_fix,'rs');
    xlabel('Vrn');
    ylabel('Y/D, X/D');
    
end


