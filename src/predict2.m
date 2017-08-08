% Calculate the intersection of mass and power surfaces for free vibration
% data.  Uses distance from mass values for minimization.

clear all;

vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};
vr_value = [4.5 5 5.5 6 6.5 7 7.5 8];
newmatfiles = {'new1p0' 'new1p22' 'new1p37' 'new1p52' 'new1p67' 'new1p9'};
phasematfiles = {'phase1p0' 'phase1p22' 'phase1p37' 'phase1p52' 'phase1p67' 'phase1p9'};
corrmatfiles = {'corr1p0' 'corr1p22' 'corr1p37' 'corr1p52' 'corr1p67' 'corr1p9'};
ampfixmatfiles = {'ampfix1p0' 'ampfix1p22' 'ampfix1p37' 'ampfix1p52' 'ampfix1p67' 'ampfix1p9'};
%newmatfiles = {'new1p9'};
%phasematfiles = {'phase1p9'};
%corrmatfiles = {'corr1p9'};
%ampfixmatfiles = {'ampfix1p9'};
freevibs = [1.0 1.22 1.37 1.52 1.67 1.9];

q = 1;

for a = 1:length(vr_set)
    eval(['load vr' char(vr_set(a)) '.mat']);
    
    check = isnan(Cmx);
    [ind1,ind2] = find(check == 1);
    Cmx(ind1,ind2) = 1;
    
    % Put force and position data in proper plotting format
    for i = 1:8
        k = 1;
        for j = 1:6
        newX(j,:,i) = allXad(i,k:k+5);
        newY(j,:,i) = allYad(i,k:k+5);
        newtheta(j,:,i) = alltheta(i,k:k+5);
        newCL3(j,:,i) = CL3(i,k:k+5);%./sqrt(2);
        newCL5(j,:,i) = CL5(i,k:k+5);%./sqrt(2);
        newCLv(j,:,i) = CLv(i,k:k+5);%./sqrt(2);
        newCmy(j,:,i) = Cmy(i,k:k+5);%./sqrt(2);
        newCmy_corr(j,:,i) = Cmy_corr(i,k:k+5);
        newCDv(j,:,i) = CDv(i,k:k+5);%./sqrt(2);
        newCmx(j,:,i) = Cmx(i,k:k+5);%./sqrt(2);
        newCDmean(j,:,i) = CDmean(i,k:k+5);
        newCD2(j,:,i) = CD2(i,k:k+5);%./sqrt(2);
        newCD4(j,:,i) = CD4(i,k:k+5);%./sqrt(2);
        newPsi(j,:,i) = Psi(i,k:k+5);
        avgPow(j,:,i) = (Plift(i,k:k+5)+Pdrag(i,k:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
        k = k+6;
        end
    end
    
    
     newXi = interp3(newX,2,'cubic');
     newYi = interp3(newY,2,'cubic');
     newthetai = interp3(newtheta,2,'cubic');
     newCmyi = interp3(newCmy,2,'cubic');
     newCmy_corri = interp3(newCmy_corr,2,'cubic');
     newCmxi = interp3(newCmx,2,'cubic');
     avgPowi = interp3(avgPow,2,'cubic');
    
    for pp = 1:length(newmatfiles);
        eval(['load ' char(newmatfiles(pp)) '.mat;']);
        vr = Vrn.*fy./(ypeakfreq./2./pi);
        fxfy = xpeakfreq./ypeakfreq;
        indfxfy = find(fxfy > 1.95 & fxfy < 2.05);
        ind_vr = find(vr(indfxfy) > vr_value(a)-0.15 & vr(indfxfy) < vr_value(a)+0.15);
             
        eval(['load ' char(phasematfiles(pp)) '.mat;']);
        phase_ind = find(phasextoy > 180);
        phasextoy(phase_ind) = phasextoy(phase_ind)-360;
        eval(['load ' char(corrmatfiles(pp)) '.mat;']);
        eval(['load ' char(ampfixmatfiles(pp)) '.mat;']);
        
        
        for r = 1:length(ind_vr)
            amy = ky./(2*pi*vels(indfxfy(ind_vr(r)))./(vr(indfxfy(ind_vr(r)))*0.0762)).^2 - my + (pi*1000*0.0762^2/4*2);
            %amx = kx./(4*ky)*(my+amy) - mx + (pi*1000*0.0762^2/4*2);
            amx = kx./(4*pi*vels(indfxfy(ind_vr(r)))./(vr(indfxfy(ind_vr(r)))*0.0762)).^2 - mx + (pi*1000*0.0762^2/4*2);
            %amy = ky./(2*pi*vels(indfxfy())./(vr(indfxfy(ind_vr(r)))*0.0762)).^2 - my + (pi*1000*0.0762^2/4*2)
            %amx = kx./(4*pi*vels(indfxfy(ind_vr(r)))./(vr(indfxfy(ind_vr(r)))*0.0762)).^2 - mx + (pi*1000*0.0762^2/4*2);
            Cmy_pred = amy./(pi*1000*0.0762^2/4*2);
            Cmx_pred = amx./(pi*1000*0.0762^2/4*2);
            mn = pi*1000*0.0762^2/4*2;
            
            %d2 = abs(newCmyi - Cmy3(indfxfy(ind_vr(r)))) + abs(newCmxi - -Cmx(indfxfy(ind_vr(r)))) + abs(avgPowi);
            %d2 = abs(newCmyi - Cmy_pred) + abs(newCmxi - Cmx_pred) + abs(avgPowi);
            
            
%             figure(1);clf;hold on;
%             h = patch(isosurface(newX,newY,newtheta,avgPow,0));
%             set(h,'FaceColor','g');
%             set(h,'FaceAlpha',0.5)
%             %h2 = patch(isosurface(newX,newY,newtheta,newCmy_corr,Cmy3(indfxfy(ind_vr(r)))));
%             h2 = patch(isosurface(newX,newY,newtheta,newCmy_corr,Cmy_pred));
%             %h2 = patch(isosurface(newX,newY,newtheta,kx/(4*ky)*(my/mn+newCmy_corr-1) - mx/mn - newCmx +1,0));
%             set(h2,'FaceColor','r');
%             set(h2,'FaceAlpha',0.5)
%             %h3 = patch(isosurface(newX,newY,newtheta,newCmx,-Cmx(indfxfy(ind_vr(r)))));
%             h3 = patch(isosurface(newX,newY,newtheta,newCmx,Cmx_pred));
%             set(h3,'FaceColor','b');
%             set(h3,'FaceAlpha',0.5)
% %             h4 = patch(isosurface(newX,newY,newtheta,3*newCL3.*cos(newPsi) ,Cmy_pred));
% %             set(h4,'FaceColor','c');
% %             set(h4,'FaceAlpha',0.5)
%             xlabel('X/D','FontSize',16);
%             ylabel('Y/D','FontSize',16);
%             zlabel('\theta','FontSize',16);
%             title('Avg Pow','FontSize',16);
%             view(-15,20);
%             grid on
            
            %print(gcf,'-djpeg100','vr6p5_fxfy1p67');
            %legend('Avg Pow = 0', sprintf('Cmy = %g',Cmy(indfxfy(ind_vr(r)))), sprintf('Cmx = %g',-Cmx(indfxfy(ind_vr(r)))));
            %legend('Avg Pow = 0', sprintf('Cmy = %g',Cmy_pred), sprintf('Cmx = %g',Cmx_pred));
            
%             Powvert = get(h,'Vertices');
%             Phasepow = Powvert(:,3)./180;
%             Powvert2 = [Powvert(:,1) Powvert(:,2) Phasepow];
%             Cmyvert = get(h2,'Vertices');
%             Phasecmy = Cmyvert(:,3)./180;
%             Cmyvert2 = [Cmyvert(:,1) Cmyvert(:,2) Phasecmy];
%             Cmxvert = get(h3,'Vertices');
%             Phasecmx = Cmxvert(:,3)./180;
%             Cmxvert2 = [Cmxvert(:,1) Cmxvert(:,2) Phasecmx];
%             
%             if isempty(Powvert2) == 1 | isempty(Cmyvert2) ==1 | isempty(Cmxvert2) == 1
%                 prediction = [0 0 0];
%             else
% 
%                 d = ones(length(Powvert2),length(Cmyvert2),length(Cmxvert2));
%                 
%                 for i = 1:length(Powvert2)
%                     for j = 1:length(Cmyvert2)
%                         for k = 1:length(Cmxvert2)
%                             if Powvert2(i,2) < 0.4 | Cmyvert2(j,2) < 0.4 | Cmxvert(k,2) < 0.4
%                                 d(i,j,k) = 1;
%                             else
%                                 d(i,j,k) = sqrt(sum((Powvert2(i,:) - Cmyvert2(j,:)).^2)) + ...
%                                     sqrt(sum((Powvert2(i,:) - Cmxvert2(k,:)).^2)) + ...
%                                     sqrt(sum((Cmyvert2(j,:) - Cmxvert2(k,:)).^2));
%                             end
%                         end
%                     end
%                 end
% 
%                 g = min(min(min(d)));
%                 ind = find(d == g);
% 
%                 [indPow,indCmy,indCmx] = ind2sub(size(d),ind(1));
% 
%                 prediction = (Powvert(indPow,:) + Cmyvert(indCmy,:) + Cmxvert(indCmx,:))/3;
%             end
            realrats = 1.5:0.1:2.5;
            realcols = [85 26 139;0 0 255;0 191 255; 84 255 159; ...
               192 255 62; 255 255 0; 255 193 37; 255 165 0; ...
               255 0 0; 255 20 147; 0 0 0]/255;
            
            %figure(2);clf;hold on;
            %h = patch(isosurface(newXi,newYi,newthetai,avgPowi,0));
            %set(h,'FaceColor','g');
            %set(h,'FaceAlpha',0.5)
            
            figure(3);clf;hold on;
            
            
            %truefxfy = sqrt((kx/(mx+ -Cmx(indfxfy(ind_vr(r)))*mn+mn))./(ky/(my+Cmy3(indfxfy(ind_vr(r)))*mn+mn)));
            
            %for co = 1:length(realrats)
                %d2 = abs(newCmy_corri - (realrats(co)^2*ky/kx*(mx + newCmxi*mn)./mn - my/mn)) + abs(avgPowi);
                %d2 = abs(newCmyi - Cmy(indfxfy(ind_vr(r)))) + 3*abs(newCmxi- -Cmx(indfxfy(ind_vr(r)))) + abs(avgPowi);
                d2 = abs(newCmyi - Cmy_pred) + 3*abs(newCmxi- Cmx_pred) + abs(avgPowi);
                %d2 = abs(newCmy_corri - (truefxfy^2*ky/kx*(mx + newCmxi*mn+mn)./mn - my/mn-1)) + abs(avgPowi);
                %    abs(newCmxi - (kx/(truefxfy^2*ky)*(my + newCmy_corri*mn+mn)./mn - mx/mn-1))+ abs(avgPowi);
                g = min(min(min(d2)));
                ind = find(d2 == g);
                %inds = find(d2 < g(co)*1.5);
                [ind1,ind2,ind3] = ind2sub(size(d2),ind(1));
                %[inds1,inds2,inds3] = ind2sub(size(d2),inds);
                %indices(co,:) = [ind1 ind2 ind3];
                prediction = [newXi(ind1,ind2,ind3) newYi(ind1,ind2,ind3) newthetai(ind1,ind2,ind3)];
% 
%             comparison = [Yad_fix(indfxfy(ind_vr(r))) Xad_fix(indfxfy(ind_vr(r))) ...
%                 phasextoy(indfxfy(ind_vr(r))) vr(indfxfy(ind_vr(r))); ...
%                 prediction(2) prediction(1) prediction(3) vr_value(a)]
%             
            ratios(q) = freevibs(pp);
            Yad_free(q) = Yad_fix(indfxfy(ind_vr(r)));
            Xad_free(q) = Xad_fix(indfxfy(ind_vr(r)));
            theta_free(q) = phasextoy(indfxfy(ind_vr(r)));
            vr_free(q) = vr(indfxfy(ind_vr(r)));
            Yad_force(q) = prediction(2);
            Xad_force(q) = prediction(1);
            theta_force(q) = prediction(3);
            vr_force(q) = vr_value(a);
            
%             figure(2);clf;hold on;
%             h = patch(isosurface(newXi,newYi,newthetai,d2,min(min(min(d2)))+0.1));
%             set(h,'FaceColor','g');
%             set(h,'FaceAlpha',0.5)
%             h2 = patch(isosurface(newXi,newYi,newthetai,d2,min(min(min(d2)))+0.5));
%             set(h2,'FaceColor','r');
%             set(h2,'FaceAlpha',0.5)
%             h3 = patch(isosurface(newXi,newYi,newthetai,d2,min(min(min(d2)))+1));
%             set(h3,'FaceColor','b');
%             set(h3,'FaceAlpha',0.5)
%             xlabel('X/D','FontSize',16);
%             ylabel('Y/D','FontSize',16);
%             zlabel('\theta','FontSize',16);
%             title('D2','FontSize',16);
%             view(-15,20);
%             grid on
            
            
            
            %for w = 1:length(inds1)
                %figure(2)
                %h2 = scatter3(newXi(inds1(w),inds2(w),inds3(w)),newYi(inds1(w),inds2(w),inds3(w)),newthetai(inds1(w),inds2(w),inds3(w)),40,realcols(co,:));
                %figure(2);
                %h2 = scatter3(newXi(ind1,ind2,ind3),newYi(ind1,ind2,ind3),newthetai(ind1,ind2,ind3),40,realcols(co,:));
                %set(h2,'MarkerFaceColor',realcols(co,:));
            %end
%             h2 = patch(isosurface(newXi,newYi,newthetai,d2,min(min(min(d2)))+0.5));
%             set(h2,'FaceColor','r');
%             set(h2,'FaceAlpha',0.5)
%             h3 = patch(isosurface(newXi,newYi,newthetai,d2,min(min(min(d2)))+1));
%             set(h3,'FaceColor','b');
%             set(h3,'FaceAlpha',0.5)
            %xlabel('X/D','FontSize',16);
            %ylabel('Y/D','FontSize',16);
            %zlabel('\theta','FontSize',16);
            %title('D2','FontSize',16);
            %view(-15,20);
            %grid on
            
            figure(3);clf;hold on;
            h = patch(isosurface(newXi,newYi,newthetai,avgPowi,0));
            set(h,'FaceColor','g');
            set(h,'FaceAlpha',0.3)
            h = patch(isosurface(newXi,newYi,newthetai,d2,min(min(min(d2)))+0.5));
            set(h,'FaceColor','b');
            set(h,'FaceAlpha',0.5)
            h2 = scatter3(Xad_fix(indfxfy(ind_vr(r))),Yad_fix(indfxfy(ind_vr(r))),phasextoy(indfxfy(ind_vr(r))),40,'r');
            set(h2,'MarkerFaceColor','r');
            %h2 = patch(isosurface(newXi,newYi,newthetai,d2,min(min(min(d2)))+0.5));
            %set(h2,'FaceColor','r');
            %set(h2,'FaceAlpha',0.5)
            %h3 = patch(isosurface(newXi,newYi,newthetai,d2,min(min(min(d2)))+1));
            %set(h3,'FaceColor','b');
            %set(h3,'FaceAlpha',0.5)
            xlabel('X/D','FontSize',16);
            ylabel('Y/D','FontSize',16);
            zlabel('\theta','FontSize',16);
            title('D2','FontSize',16);
            view(-15,20);
            grid on
            %pause
            %end
            
            truefxfy = sqrt((kx/(mx+ -Cmx(indfxfy(ind_vr(r)))*mn))./(ky/(my+Cmy(indfxfy(ind_vr(r)))*mn)));
            %num = min(g(1));
            %num2 = find(g == num);
            %num2 = find(abs(realrats-truefxfy) < 0.05);
            %index = indices(num2(1),:);
            %prediction = [newXi(index(1),index(2),index(3)) ...
            %    newYi(index(1),index(2),index(3)) newthetai(index(1),index(2),index(3))];
            %comparison = [Yad_fix(indfxfy(ind_vr(r))) Xad_fix(indfxfy(ind_vr(r))) ...
            %    phasextoy(indfxfy(ind_vr(r))) vr(indfxfy(ind_vr(r))) truefxfy; ...
            %    prediction(2) prediction(1) prediction(3) vr_value(a) realrats(num2(1))]
            
            clear Powvert Cmyvert Cmxvert indPow indCmy indCmx d d2 g
            %pause
            q = q+1;
        end
        %clear vr fxfy indfxfy ind_vr phase_ind
    end
        
     clear C* all* newX newY newtheta newC* X Y avgPow  
end

ind1p9 = find(ratios == 1.9);
ind1p67 = find(ratios == 1.67);
ind1p52 = find(ratios == 1.52);
ind1p37 = find(ratios == 1.37);
ind1p22 = find(ratios == 1.22);
ind1p0 = find(ratios == 1);

indices = {'ind1p0' 'ind1p22' 'ind1p37' 'ind1p52' 'ind1p67' 'ind1p9'};
cols_free = {'ko' 'ko' 'ko' 'ko' 'ko' 'ko'};
cols_force = {'ks' 'ks' 'ks' 'ks' 'ks' 'ks'};
tit = {'f_{x}/f_{y} = 1.0' 'f_{x}/f_{y} = 1.22' 'f_{x}/f_{y} = 1.37' ...
    'f_{x}/f_{y} = 1.52' 'f_{x}/f_{y} = 1.67' 'f_{x}/f_{y} = 1.9'};

fxfyname = [0 22 37 51 67 9];

figure(1); clf;

for pp = 1:length(newmatfiles)
%     eval(['load ' char(newmatfiles(p)) '.mat;']);
%     vr = Vrn.*fy./(ypeakfreq./2./pi);
%     eval(['load ' char(phasematfiles(p)) '.mat;']);
%     phase_ind = find(phasextoy > 180);
%     phasextoy(phase_ind) = phasextoy(phase_ind)-360;
    
    figure(1); clf; hold on;
    %eval(['plot(vr(7:27),Yad(7:27), char(cols_free(p)));']);
    eval(['plot(vr_free(' char(indices(pp)) '),Yad_free(' char(indices(pp)) '), char(cols_free(pp)));']);
    eval(['plot(vr_force(' char(indices(pp)) '),Yad_force(' char(indices(pp)) '), char(cols_force(pp)));']);
    xlabel('Vr');
    ylabel('Y/D');
    axis([4 9 0 1.5]);
    title(char(tit(pp)));
    legend('Observed','Predicted','Location', 'Northwest');
    
    eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(pp)) '_Yad'')']);
    
    figure(2); clf; hold on;
    %eval(['plot(vr(7:27),Yad(7:27), char(cols_free(p)));']);
    eval(['plot(vr_free(' char(indices(pp)) '),Xad_free(' char(indices(pp)) '), char(cols_free(pp)));']);
    eval(['plot(vr_force(' char(indices(pp)) '),Xad_force(' char(indices(pp)) '), char(cols_force(pp)));']);
    xlabel('Vr');
    ylabel('X/D');
    axis([4 9 0 0.6]);
    title(char(tit(pp)));
    legend('Observed','Predicted','Location', 'Northwest');
    
    eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(pp)) '_Xad'')']);
    
    figure(3); clf; hold on;
    %eval(['plot(vr(7:27),Yad(7:27), char(cols_free(p)));']);
    eval(['plot(vr_free(' char(indices(pp)) '),theta_free(' char(indices(pp)) '), char(cols_free(pp)));']);
    eval(['plot(vr_force(' char(indices(pp)) '),theta_force(' char(indices(pp)) '), char(cols_force(pp)));']);
    xlabel('Vr');
    ylabel('\theta');
    axis([4 9 -180 180]);
    title(char(tit(pp)));
    legend('Observed','Predicted','Location', 'Northwest');
    
    eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(pp)) '_Theta'')']);
    
    pause
end
