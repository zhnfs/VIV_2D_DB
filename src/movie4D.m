% Make movie file for 4D data plots

clear all

vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
vr_value = [4.5 5 5.5 6 6.5 7 7.5 8];

for a = 1:length(vr_set)
    
    % % Sam added an input command to this file
    % vr_set = input('Reduced velocity set to run (e.g. 5, 5p5,...) :','s');
    % 
    % % load vr6p5.mat
    % % sam's attempt
    eval(['load vr' char(vr_set(a)) '.mat']);
    % 
    % %newCL3 = zeros(8*36,1);

    for i = 1:8
        k = 1;
        for j = 1:6
            X(j,:,i) = allXad(i,k:k+5);
            Y(j,:,i) = allYad(i,k:k+5);
            newtheta(j,:,i) = alltheta(i,k:k+5);
            newCL3(j,:,i) = CL3(i,k:k+5);
            newCL5(j,:,i) = CL5(i,k:k+5);
            newCLv(j,:,i) = CLv(i,k:k+5);
            newCLvfft(j,:,i) = CLvfft(i,k:k+5);
            newCmy(j,:,i) = Cmy(i,k:k+5);
            newCDv(j,:,i) = CDv(i,k:k+5);
            newCmx(j,:,i) = Cmx(i,k:k+5);
            newCDmean(j,:,i) = CDmean(i,k:k+5);
            newCD2(j,:,i) = CD2(i,k:k+5);
            newCD4(j,:,i) = CD4(i,k:k+5);
            if strcmp(vr_set,'4p5') == 1 | strcmp(vr_set,'5') == 1 | strcmp(vr_set,'5p5') == 1
                avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.18.^3*0.6858*0.0381);
            else
                avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
            end
        k = k+6;
        end
        %newCL3(j:j+35) = CL3(i,:);
        %newYad(j:j+35) = allYad(i,:);
        %newXad(j:j+35) = allXad(i,:);
        %newtheta(j:j+35) = alltheta(i,:);
        %j = j+36;
    end

    %conts = [-2 -1 0 0.1];
    conts = 0;
    %cols = {'m' 'b' 'g' 'y' 'r'};
    cols = 'g';

    figure(a);hold on;
    for i = 1:length(conts)
        h = patch(isosurface(X,Y,newtheta,avgPow,conts(i)));
        eval(['set(h,''FaceColor'',''' char(cols(i)) ''');']);
        set(h,'FaceAlpha',0.3)
    end
    xlabel('X/D','FontSize',16);
    ylabel('Y/D','FontSize',16);
    zlabel('\theta','FontSize',16);
    title(sprintf('V_{r} = %g',vr_value(a)),'FontSize',16);
    axis([0 0.75 0 1.5 -180 180])
    set(gcf,'Position',[200 200 800 600])
    grid on


    newmatfiles = {'new1p0' 'new1p22' 'new1p37' 'new1p52' 'new1p67' 'new1p9'};
    %newmatfiles = {'new1p37' 'new1p52' 'new1p67' 'new1p9'};
    phasematfiles = {'phase1p0' 'phase1p22' 'phase1p37' 'phase1p52' 'phase1p67' 'phase1p9'};
    ampfixmatfiles = {'ampfix1p0' 'ampfix1p22' 'ampfix1p37' 'ampfix1p52' 'ampfix1p67' 'ampfix1p9'};
    %phasematfiles = {'phase1p37' 'phase1p52' 'phase1p67' 'phase1p9'};
    freevibs = [1.0 1.22 1.37 1.52 1.67 1.9];
    %leg = [conts freevibs];
%     leg = {'-2' '-1' '0' '0.1' 'f_{nx}/f_{ny} = 1.0' 'f_{nx}/f_{ny} = 1.22' ...
%         'f_{nx}/f_{ny} = 1.37' 'f_{nx}/f_{ny} = 1.52' 'f_{nx}/f_{ny} = 1.67' ...
%         'f_{nx}/f_{ny} = 1.9' 'Flex Beam Cyl' 'Jauvtis & Williamson (2004)'};
    leg = {'P_{avg} = 0' 'f_{nx}/f_{ny} = 1.0' 'f_{nx}/f_{ny} = 1.22' ...
        'f_{nx}/f_{ny} = 1.37' 'f_{nx}/f_{ny} = 1.52' 'f_{nx}/f_{ny} = 1.67' ...
        'f_{nx}/f_{ny} = 1.9' 'Flex Beam Cyl' 'Jauvtis & Williamson (2004)'};
    freecols = [1 0 0;0 1 1; 0 0 0;255/255 140/255 0;25/255 25/255 112/255; 255/255 20/255 147/255];

    for p = 1:length(newmatfiles);
        eval(['load ' char(newmatfiles(p)) '.mat;']);
        vr = Vrn.*fy./(ypeakfreq./2./pi);
        fxfy = xpeakfreq./ypeakfreq;
        indfxfy = find(fxfy > 1.95 & fxfy < 2.05);
        ind_vr = find(vr(indfxfy) > vr_value(a)-0.1 & vr(indfxfy) < vr_value(a)+0.1);
        eval(['load ' char(phasematfiles(p)) '.mat;']);
        phase_ind = find(phasextoy > 180);
        phasextoy(phase_ind) = phasextoy(phase_ind)-360;
        eval(['load ' char(ampfixmatfiles(p)) '.mat;']);
        %h2 = scatter3(Xad(indfxfy(ind_vr))./sqrt(2),Yad(indfxfy(ind_vr))./sqrt(2),phasextoy(indfxfy(ind_vr)));
        h2 = scatter3(Xad(indfxfy(ind_vr)),Yad(indfxfy(ind_vr)),phasextoy(indfxfy(ind_vr)));
        set(h2,'MarkerFaceColor',freecols(p,:),'MarkerEdgeColor',freecols(p,:));
    end

    load beamcyl
        ind_vr_beam = find(vr_all > vr_value(a)-0.1 & vr_all < vr_value(a)+0.1);
        h3 = scatter3(Xad_all(ind_vr_beam),Yad_all(ind_vr_beam),phasexy(ind_vr_beam));
        set(h3,'MarkerFaceColor',[131/255 111/255 255/255],'MarkerEdgeColor',[131/255 111/255 255/255]);
        load jauvtis
        load fstar
        for cc = 1:length(theta_jauv(:,2))
            if theta_jauv(cc,2) > 180
                theta_jauv(cc,2) = theta_jauv(cc,2)-360;
            else
                theta_jauv(cc,2) = theta_jauv(cc,2);
            end
        end
        ind_vr_jauv = find(Yad_jauv(:,1)./(fstar(:,2)) > vr_value(a)-0.1 & Yad_jauv(:,1)./(fstar(:,2)) < vr_value(a)+0.1);
        h4 = scatter3(Xad_jauv(ind_vr_jauv,2),Yad_jauv(ind_vr_jauv,2),theta_jauv(ind_vr_jauv,2));
        set(h4,'MarkerFaceColor',[93/255 71/255 139/255],'MarkerEdgeColor',[93/255 71/255 139/255]);
    %legend(num2str(leg'));
    legend(leg', 'Location','EastOutside');
    view(-60,30);
    
    %print(gcf,'-depsc',['comparison_vr' char(vr_set(a))])
    %print(gcf,'-djpeg100',['comparison_vr' char(vr_set(a))])
    
    % q = 1;
    %for i = -130:10:130
    %view(i,30);
    %movieims(q) = getframe(gcf);
    %q = q+1;
    %end
    
    %movname = ['movie_compare_vr' char(vr_set(a))];
    %mov= avifile(movname,'compression','Cinepak','FPS',5);
    %mov = addframe(mov,movieims);
    %mov= close(mov);
    clear C* all* X* Y*
end;
