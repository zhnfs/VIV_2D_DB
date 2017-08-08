% Try to predict motions based on structural properties

clear all

cd('/Users/haining/VIV/src/ILCFprediction/databaseJMD');

vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};

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


for p = 1:length(vr_set)
    eval(['load vr' char(vr_set(p)) '.mat']);
    
    [ind1,ind2] = find(isnan(Cmx));
    Cmx(ind1,ind2) = 1;
    
    for i = 1:8
            k = 1;
            for j = 1:6
            newX(j,:,i) = allXad(i,k:k+5);
            newY(j,:,i) = allYad(i,k:k+5);
            newtheta(j,:,i) = alltheta(i,k:k+5);
            vr(j,:,i) = Vr(i,k:k+5);
            newCL3(j,:,i) = CL3(i,k:k+5);
            newCL1(j,:,i) = CL1(i,k:k+5);
            newCD2(j,:,i) = CD2(i,k:k+5);
            newCmy(j,:,i) = Cmy(i,k:k+5);
            newCmx(j,:,i) = Cmx(i,k:k+5);
            newCLv(j,:,i) = CLv(i,k:k+5);
            newCDv(j,:,i) = CDv(i,k:k+5);
            
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
    newtheta(:,:,9) = 180;
    vr(:,:,9) = vr(:,:,1);
    avgPow(:,:,9) = avgPow(:,:,1);
    newCL3(:,:,9) = newCL3(:,:,1);
    newCL1(:,:,9) = newCL1(:,:,1);
    newCD2(:,:,9) = newCD2(:,:,1);
    newCmy(:,:,9) = newCmy(:,:,1);
    newCmx(:,:,9) = newCmx(:,:,1);
    newCLv(:,:,9) = newCLv(:,:,1);
    newCDv(:,:,9) = newCDv(:,:,1);
        
    bigX(:,:,:,p) = newX;
    bigY(:,:,:,p) = newY;
    bigtheta(:,:,:,p) = newtheta;
    bigvr(:,:,:,p) = vr;        
    bigCL1(:,:,:,p) = newCL1;
    bigCL3(:,:,:,p) = newCL3;
    bigCD2(:,:,:,p) = newCD2;
    bigpow(:,:,:,p) = avgPow;
    bigCmy(:,:,:,p) = newCmy;
    bigCmx(:,:,:,p) = newCmx;
    bigCLv(:,:,:,p) = newCLv;
    bigCDv(:,:,:,p) = newCDv;
end

Y_vec = [0.25:0.25:1.5];
X_vec = [0:0.15:0.75];
phase_vec = [-180 -135 -90 -45 0 45 90 135 180];
vr_vec = [4.5:0.5:8];

save dbJMDori.mat bigX bigY bigtheta bigCL1 bigCD2 bigpow bigCmy bigCmx bigCLv bigCDv Y_vec X_vec phase_vec vr_vec

% interpolation 4D matrix from 6*6*9*8 to  21*21*33*29 if ifac=2
%   VI = INTERPN(V,NTIMES) expands V by interleaving interpolates between
%   every element, working recursively for NTIMES.  
%   INTERPN requires that X1,X2,X3,etc. be monotonic and plaid (as if
%   they were created using NDGRID).  X1,X2,X3,etc. can be non-uniformly
%   spaced.

ifac = 2;
int = 'spline';

int2 = 'linear';
bigXi = interpn(bigX,ifac,int2);
bigYi = interpn(bigY,ifac,int2);
bigvri = interpn(bigvr,ifac,int2);
bigthetai = interpn(bigtheta,ifac,int2);

bigCL1i = interpn(bigCL1,ifac,int);
bigCL3i = interpn(bigCL3,ifac,int);
bigCD2i = interpn(bigCD2,ifac,int);
bigpowi = interpn(bigpow,ifac,int);
bigCmyi = interpn(bigCmy,ifac,int);
bigCmxi = interpn(bigCmx,ifac,int);


% 
% Y_vec = [0.25:0.25:1.5];
% X_vec = [0:0.15:0.75];
% phase_vec = [-180 -135 -90 -45 0 45 90 135 180];
% vr_vec = [4.5:0.5:8];
% 
% bigCL1i_ind = interpn(Y_vec,X_vec, phase_vec, vr_vec, bigCL1,0.5,0.15,-90,6,int);


cd('/Users/haining/VIV/src/ILCFprediction/freeVibJMD');

newmatfiles = {'new1p0' 'new1p22' 'new1p37' 'new1p52' 'new1p67' 'new1p9'};
phasematfiles = {'phase1p0' 'phase1p22' 'phase1p37' 'phase1p52' 'phase1p67' 'phase1p9'};
corrmatfiles = {'corr1p0' 'corr1p22' 'corr1p37' 'corr1p52' 'corr1p67' 'corr1p9'};
ampfixmatfiles = {'ampfix1p0' 'ampfix1p22' 'ampfix1p37' 'ampfix1p52' 'ampfix1p67' 'ampfix1p9'};

for pp = 6:length(newmatfiles)

    eval(['load ' char(newmatfiles(pp)) '.mat;']);
    vr_free = Vrn.*fy./(ypeakfreq./2./pi);
    amy = ky./(ypeakfreq).^2 - my + (pi*1000*0.0762^2/4*2);
    amx = kx./(2*ypeakfreq).^2 - mx + (pi*1000*0.0762^2/4*2);
    Cmy_pred = amy./(pi*1000*0.0762^2/4*2);
    Cmx_pred = amx./(pi*1000*0.0762^2/4*2);
    mn = pi*1000*0.0762^2/4*2;
    
    eval(['load ' char(phasematfiles(pp)) '.mat;']);
    phase_ind = find(phasextoy > 180);
    phasextoy(phase_ind) = phasextoy(phase_ind)-360;
    eval(['load ' char(corrmatfiles(pp)) '.mat;']);
    eval(['load ' char(ampfixmatfiles(pp)) '.mat;']);
    
    cp = (0.5*Yad_fix.^2*0.0762^2.*ypeakfreq'.^2*2*my*2*pi*fy*zetay + ...
        0.5*Xad_fix.^2*0.0762^2.*xpeakfreq'.^2*2*mx*2*pi*fx*zetax)./(0.5*1000*0.0762*2*vels.^3);
    
    for r = 9:25
        dvr = abs(bigvri-vr_free(r));
        indices = find(dvr < 1*0.5/(ifac^2));
        [p1,p2,p3,p4] = ind2sub(size(bigXi),indices);
        d = 100*ones(size(bigXi));
        
        for z = 1:length(p1)
            d(p1(z),p2(z),p3(z),p4(z)) = abs((bigCmyi(p1(z),p2(z),p3(z),p4(z)) - ...
                Cmy_pred(r))) + 5*abs((bigCmxi(p1(z),p2(z),p3(z),p4(z)) - ...
                Cmx_pred(r))) + abs(bigpowi(p1(z),p2(z),p3(z),p4(z)));
        end
        if min(min(min(min(d)))) == 100
            vr_force(pp,r) = 0;
            Yad_force(pp,r) = 0;
            Xad_force(pp,r) = 0;
            Theta_force(pp,r) = 0;
            vrn_force(pp,r) = 0;
        else
            mins = find(d == min(min(min(min(d)))));
%             index(pp,r) = find(d == min(min(min(min(d)))));
            index(pp,r) = mins(1);
            [i1,i2,i3,i4] = ind2sub(size(bigXi),index(pp,r));
            vr_force(pp,r) = bigvri(i1,i2,i3,i4);
            Yad_force(pp,r) = bigYi(i1,i2,i3,i4);
            Xad_force(pp,r) = bigXi(i1,i2,i3,i4);
            Theta_force(pp,r) = bigthetai(i1,i2,i3,i4);
            vr_fr(pp,r) = vr_free(r);
            Yad_fr(pp,r) = Yad_fix(r);
            Xad_fr(pp,r) = Xad_fix(r);
            Theta_fr(pp,r) = phasextoy(r);
            vrn_force(pp,r) = vr_force(pp,r).*(ypeakfreq(r)/2/pi)/fy;
            vrn_free(pp,r) = Vrn(r);
        end
    end
end

fxfyrat = [1 1.22 1.37 1.52 1.67 1.9];
fxfyname = [0 22 37 52 67 9];

Yad_err = abs(Yad_force - Yad_fr)./abs(Yad_fr)
Xad_err = abs(Xad_force - Xad_fr)./abs(Xad_fr)
Theta_err = abs(Theta_force - Theta_fr)./abs(Theta_fr)

for g = 6:6

%     figure(1); clf;
%     plot(vr_fr(g,9:end),Yad_fr(g,9:end),'ko-',vr_force(g,9:end),Yad_force(g,9:end),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(g)),'FontSize',16);
%     xlabel('V_{r}','FontSize',16);
%     ylabel('A_{y}/D','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib','FontSize',14,'Location','Southwest')
%     axis([4 9 0 1.5]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
%     grid
% 
%     %eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(g)) '_Yad_1p9'')']);
%     %eval(['print(gcf,''-djpeg100'', ''fxfy1p' num2str(fxfyname(g)) '_Yad_1p9'')']);
%     
%     figure(2); clf;
%     plot(vr_fr(g,9:end),Xad_fr(g,9:end),'ko-',vr_force(g,9:end),Xad_force(g,9:end),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(g)),'FontSize',16);
%     xlabel('V_{r}','FontSize',16);
%     ylabel('A_{x}/D','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib', 'Predicted from Forced Vib','FontSize',14,'Location','Northwest')
%     axis([4 9 0 0.6]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9],'Ytick',[0 0.15 0.3 0.45 0.6]);
%     grid
%     
%     %eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(g)) '_Xad_1p9'')']);
%     %eval(['print(gcf,''-djpeg100'', ''fxfy1p' num2str(fxfyname(g)) '_Xad_1p9'')']);
% 
% 
%     figure(3); clf;
%     plot(vr_fr(g,9:end),Theta_fr(g,9:end),'ko-',vr_force(g,9:end),Theta_force(g,9:end),'rs','MarkerSize',9);
%     %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(g)),'FontSize',16);
%     xlabel('V_{r}','FontSize',16);
%     ylabel('\theta','FontSize',16);
%     set(gca,'FontSize',14);
%     legend('Observed Free Vib','Predicted from Forced Vib','FontSize',14,'Location','Northwest')
%     axis([4 9 -180 180]);
%     set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
%     grid
%     
%     %eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(g)) '_Theta_1p9'')']);
%     %eval(['print(gcf,''-djpeg100'', ''fxfy1p' num2str(fxfyname(g)) '_Theta_1p9'')']);
%     
%     %pause;
    
    figure(1); clf;
    plot(vrn_free(g,9:end),Yad_fr(g,9:end),'ko-',vrn_force(g,9:end),Yad_force(g,9:end),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(g)),'FontSize',16);
    xlabel('V_{r}','FontSize',16);
    ylabel('A_{y}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','FontSize',14,'Location','Southwest')
    axis([4 9 0 1.5]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9],'Ytick',[0 0.25 0.5 0.75 1 1.25 1.5]);
    grid

    %eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(g)) '_Yad_1p9'')']);
    %eval(['print(gcf,''-djpeg100'', ''fxfy1p' num2str(fxfyname(g)) '_Yad_1p9'')']);
    
    figure(2); clf;
    plot(vrn_free(g,9:end),Xad_fr(g,9:end),'ko-',vrn_force(g,9:end),Xad_force(g,9:end),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(g)),'FontSize',16);
    xlabel('V_{r}','FontSize',16);
    ylabel('A_{x}/D','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib', 'Predicted from Forced Vib','FontSize',14,'Location','Northwest')
    axis([4 9 0 0.6]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9],'Ytick',[0 0.15 0.3 0.45 0.6]);
    grid
    
    %eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(g)) '_Xad_1p9'')']);
    %eval(['print(gcf,''-djpeg100'', ''fxfy1p' num2str(fxfyname(g)) '_Xad_1p9'')']);


    figure(3); clf;
    plot(vrn_free(g,9:end),Theta_fr(g,9:end),'ko-',vrn_force(g,9:end),Theta_force(g,9:end),'rs','MarkerSize',9);
    %title(sprintf('Prediction Comparison for f_{x}/f_{y} = %2.3g', fxfyrat(g)),'FontSize',16);
    xlabel('V_{r}','FontSize',16);
    ylabel('\theta','FontSize',16);
    set(gca,'FontSize',14);
    legend('Observed Free Vib','Predicted from Forced Vib','FontSize',14,'Location','Northwest')
    axis([4 9 -180 180]);
    set(gca,'Xtick',[4 4.5 5 5.5 6 6.5 7 7.5 8 8.5 9],'Ytick',[-180 -135 -90 -45 0 45 90 135 180]);
    grid
    
    %eval(['print(gcf,''-depsc'', ''fxfy1p' num2str(fxfyname(g)) '_Theta_1p9'')']);
    %eval(['print(gcf,''-djpeg100'', ''fxfy1p' num2str(fxfyname(g))
    %'_Theta_1p9'')']);
    
end

cd('/Users/haining/VIV/src/ILCFprediction/src')
% pX(:,:,:) = bigX(:,:,1,:);
% pY(:,:,:) = bigY(:,:,1,:);
% pvr(:,:,:) = bigvr(:,:,1,:);
% ppow(:,:,:) = bigpow(:,:,1,:);
% 
%  figure(4);
% % h1 = patch(isosurface(pX,pY,pvr,ppow,0));
% % set(h1,'FaceColor','r');
% % set(h1,'FaceAlpha',0.3)
% 
% pX(:,:,:) = bigX(:,:,2,:);
% pY(:,:,:) = bigY(:,:,2,:);
% pvr(:,:,:) = bigvr(:,:,2,:);
% ppow(:,:,:) = bigpow(:,:,2,:);
% 
% h2 = patch(isosurface(pX,pY,pvr,ppow,0));
% set(h2,'FaceColor',[255/255 140/255 0]);
% set(h2,'FaceAlpha',0.3)
% 
% pX(:,:,:) = bigX(:,:,3,:);
% pY(:,:,:) = bigY(:,:,3,:);
% pvr(:,:,:) = bigvr(:,:,3,:);
% ppow(:,:,:) = bigpow(:,:,3,:);
% 
% h3 = patch(isosurface(pX,pY,pvr,ppow,0));
% set(h3,'FaceColor','y');
% set(h3,'FaceAlpha',0.3)
% 
% 
% pX(:,:,:) = bigX(:,:,4,:);
% pY(:,:,:) = bigY(:,:,4,:);
% pvr(:,:,:) = bigvr(:,:,4,:);
% ppow(:,:,:) = bigpow(:,:,4,:);
% 
% h4 = patch(isosurface(pX,pY,pvr,ppow,0));
% set(h4,'FaceColor','g');
% set(h4,'FaceAlpha',0.3)

% pX(:,:,:) = bigX(:,:,5,:);
% pY(:,:,:) = bigY(:,:,5,:);
% pvr(:,:,:) = bigvr(:,:,5,:);
% ppow(:,:,:) = bigpow(:,:,5,:);
% 
% h5 = patch(isosurface(pX,pY,pvr,ppow,0));
% set(h5,'FaceColor','c');
% set(h5,'FaceAlpha',0.3)


% pX(:,:,:) = bigX(:,:,6,:);
% pY(:,:,:) = bigY(:,:,6,:);
% pvr(:,:,:) = bigvr(:,:,6,:);
% ppow(:,:,:) = bigpow(:,:,6,:);
% 
% h6 = patch(isosurface(pX,pY,pvr,ppow,0));
% set(h6,'FaceColor','b');
% set(h6,'FaceAlpha',0.3)

% pX(:,:,:) = bigX(:,:,7,:);
% pY(:,:,:) = bigY(:,:,7,:);
% pvr(:,:,:) = bigvr(:,:,7,:);
% ppow(:,:,:) = bigpow(:,:,7,:);
% 
% h7 = patch(isosurface(pX,pY,pvr,ppow,0));
% set(h7,'FaceColor','m');
% set(h7,'FaceAlpha',0.3)

% pX(:,:,:) = bigX(:,:,8,:);
% pY(:,:,:) = bigY(:,:,8,:);
% pvr(:,:,:) = bigvr(:,:,8,:);
% ppow(:,:,:) = bigpow(:,:,8,:);
% 
% h8 = patch(isosurface(pX,pY,pvr,ppow,0));
% set(h8,'FaceColor',[125/255 38/255 205/255]);
% set(h8,'FaceAlpha',0.3)

% for i = 1:length(Xad_force(6,9:end))
%     if Theta_force(6,8+i) > -180 & Theta_force(6,8+i) < -157.4
%         theta_for_col(i,:) = [1 0 0];
%     elseif Theta_force(6,8+i) > -157.4 & Theta_force(6,8+i) < -112.4
%         theta_for_col(i,:) = [255/255 140/255 0];
%     elseif Theta_force(6,8+i) > -112.4 & Theta_force(6,8+i) < -67.4
%         theta_for_col(i,:) = [1 1 0];
%     elseif Theta_force(6,8+i) > -67.4 & Theta_force(6,8+i) < -22.4
%         theta_for_col(i,:) = [0 1 0];
%     elseif Theta_force(6,8+i) > -22.4 & Theta_force(6,8+i) < 22.4
%         theta_for_col(i,:) = [0 1 1];
%     elseif Theta_force(6,8+i) > 22.4 & Theta_force(6,8+i) < 67.4
%         theta_for_col(i,:) = [0 0 1];
%     elseif Theta_force(6,8+i) > 67.4 & Theta_force(6,8+i) < 112.4
%         theta_for_col(i,:) = [1 0 1];
%     elseif Theta_force(6,8+i) > 112.4 & Theta_force(6,8+i) < 135.4
%         theta_for_col(i,:) = [125/255 38/255 205/255];
%     end
%     
%     %hold on
%     %h9 = plot3(Xad_force(6,8+i),Yad_force(6,8+i),vr_force(6,8+i),'o');
%     %set(h9,'MarkerFaceColor',theta_for_col(i,:),'MarkerEdgeColor',theta_for_col(i,:));
% end

% for i = 1:length(Xad_fr(6,11:end))
%     if Theta_fr(6,10+i) > -180 & Theta_fr(6,10+i) < -157.5
%         theta_fr_col(i,:) = [1 0 0];
%     elseif Theta_fr(6,10+i) > -157.4 & Theta_fr(6,10+i) < -112.4
%         theta_fr_col(i,:) = [255/255 140/255 0];
%     elseif Theta_fr(6,10+i) > -112.4 & Theta_fr(6,10+i) < -67.4
%         theta_fr_col(i,:) = [1 1 0];
%     elseif Theta_fr(6,10+i) > -67.4 & Theta_fr(6,10+i) < -22.4
%         theta_fr_col(i,:) = [0 1 0];
%     elseif Theta_fr(6,10+i) > -22.4 & Theta_fr(6,10+i) < 22.4
%         theta_fr_col(i,:) = [0 1 1];
%     elseif Theta_fr(6,10+i) > 22.4 & Theta_fr(6,10+i) < 67.4
%         theta_fr_col(i,:) = [0 0 1];
%     elseif Theta_fr(6,10+i) > 67.4 & Theta_fr(6,10+i) < 112.4
%         theta_fr_col(i,:) = [1 0 1];
%     elseif Theta_fr(6,10+i) > 112.4 & Theta_fr(6,10+i) < 135.4
%         theta_fr_col(i,:) = [125/255 38/255 205/255];
%     end
%     
%     hold on
%     h10 = plot3(Xad_fr(6,10+i),Yad_fr(6,10+i),vr_fr(6,10+i),'sq');
%     set(h10,'MarkerFaceColor',theta_fr_col(i,:),'MarkerEdgeColor',theta_fr_col(i,:));
%     
% end

%hold on
%h9 = plot3(Xad_force(6,9:end),Yad_force(6,9:end),vr_force(6,9:end),'k-');
%set(h9,'MarkerFaceColor','k','MarkerEdgeColor','k');

%h10 = plot3(Xad_fr(6,11:end),Yad_fr(6,11:end),vr_fr(6,11:end),'r-','LineWidth',2);
%set(h10,'MarkerFaceColor','r','MarkerEdgeColor','r');

% load jauvtis
% load fstar
% 
% for i = 1:length(theta_jauv(:,2))
%     if theta_jauv(i,2) > 180
%         theta_jauv(i,2) = theta_jauv(i,2) - 360
%     else
%         theta_jauv(i,2) = theta_jauv(i,2);
%     end
% end
% 
% for i = 1:length(theta_jauv(:,2))
%     if theta_jauv(i,2) > 157.4 | theta_jauv(i,2) < -157.4
%         theta_jauv_col(i,:) = [1 0 0];
%     elseif theta_jauv(i,2) > -157.4 & theta_jauv(i,2) < -112.4
%         theta_jauv_col(i,:) = [255/255 140/255 0];
%     elseif theta_jauv(i,2) > -112.4 & theta_jauv(i,2) < -67.4
%         theta_jauv_col(i,:) = [1 1 0];
%     elseif theta_jauv(i,2) > -67.4 & theta_jauv(i,2) < -22.4
%         theta_jauv_col(i,:) = [0 1 0];
%     elseif theta_jauv(i,2) > -22.4 & theta_jauv(i,2) < 22.4
%         theta_jauv_col(i,:) = [0 1 1];
%     elseif theta_jauv(i,2) > 22.4 & theta_jauv(i,2) < 67.4
%         theta_jauv_col(i,:) = [0 0 1];
%     elseif theta_jauv(i,2) > 67.4 & theta_jauv(i,2) < 112.4
%         theta_jauv_col(i,:) = [1 0 1];
%     elseif theta_jauv(i,2) > 112.4 & theta_jauv(i,2) < 157.4
%         theta_jauv_col(i,:) = [125/255 38/255 205/255];
%     end
%     
%     hold on
%     h10 = plot3(Xad_jauv(i,2),Yad_jauv(i,2),fstar(i,1)./fstar(i,2),'sq');
%     set(h10,'MarkerFaceColor',theta_jauv_col(i,:),'MarkerEdgeColor',theta_jauv_col(i,:));
%     
% end
% 
% 
% h11 = plot3(Xad_jauv(:,2),Yad_jauv(:,2),fstar(:,1)./fstar(:,2),'r-','LineWidth',2);
% %set(h11,'MarkerFaceColor','r','MarkerEdgeColor','r');
% 
% axis([0 0.75 0 1.5 4 8.5])
% xlabel('A_{x}/D','FontSize',14);
% ylabel('A_{y}/D','FontSize',14);
% zlabel('V_{r}','FontSize',14);
% legend('-135 deg','-90 deg','-45 deg','135 deg')
% grid on
% view(40,15)
% 
% 
% 
% print(gcf,'-depsc','3djauv_view1');
% print(gcf,'-djpeg100','3djauv_view1');
% 
% view(-170,10)
% print(gcf,'-depsc','3djauv_view2');
% print(gcf,'-djpeg100','3djauv_view2');
% 
% view(0,90)
% print(gcf,'-depsc','3djauv_view3');
% print(gcf,'-djpeg100','3djauv_view3');
% 
% view(140,-25)
% print(gcf,'-depsc','3djauv_view4');
% print(gcf,'-djpeg100','3djauv_view4');

% h11 = plot3(Xad_fr(5,9:25),Yad_fr(5,9:25),vr_fr(5,9:25),'bo-');
% set(h11,'MarkerFaceColor','b','MarkerEdgeColor','b');
% axis([0 0.75 0 1.5 4 8.5])
% 
% h12 = plot3(Xad_fr(4,9:25),Yad_fr(4,9:25),vr_fr(4,9:25),'go-');
% set(h12,'MarkerFaceColor','g','MarkerEdgeColor','g');
% axis([0 0.75 0 1.5 4 8.5])

% h13 = plot3(Xad_fr(3,9:25),Yad_fr(3,9:25),vr_fr(3,9:25),'ko-');
% set(h13,'MarkerFaceColor','k','MarkerEdgeColor','k');
% axis([0 0.75 0 1.5 4 8.5])
% 
% h14 = plot3(Xad_fr(2,9:25),Yad_fr(2,9:25),vr_fr(2,9:25),'mo-');
% set(h14,'MarkerFaceColor','m','MarkerEdgeColor','m');
% axis([0 0.75 0 1.5 4 8.5])
% 
% h15 = plot3(Xad_fr(1,9:25),Yad_fr(1,9:25),vr_fr(1,9:25),'co-');
% set(h15,'MarkerFaceColor','c','MarkerEdgeColor','c');
% axis([0 0.75 0 1.5 4 8.5])