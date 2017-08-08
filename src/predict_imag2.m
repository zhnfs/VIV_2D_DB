% Try to predict motions based on structural properties

clear all

vr_set = {'4p5' '5' '5p5' '6' '6p5' '7' '7p5' '8'};

ifac = 2;
int = 'spline';

bigX = zeros(6,5,9,8);
bigY = zeros(6,5,9,8);
bigtheta = zeros(6,5,9,8);
bigCL1 = zeros(6,5,9,8);
bigCL3 = zeros(6,5,9,8);
bigCD2 = zeros(6,5,9,8);
bigvr = zeros(6,5,9,8);
bigfft_posy = zeros(6,5,9,8);
bigfft_vely = zeros(6,5,9,8);
bigfft_accy = zeros(6,5,9,8);
bigfft_posx = zeros(6,5,9,8);
bigfft_velx = zeros(6,5,9,8);
bigfft_accx = zeros(6,5,9,8);
bigfft_cdf = zeros(6,5,9,8);
bigfft_cl1 = zeros(6,5,9,8);
bigfft_cl3 = zeros(6,5,9,8);
bigpow = zeros(6,5,9,8);


for p = 1:length(vr_set)
    eval(['load vr' char(vr_set(p)) '.mat']);
    
    %[ind1,ind2] = find(isnan(fft_xpos));
    %fft_xpos(ind1,ind2) = 1;
    %fft_xvel(ind1,ind2) = 1;
    %fft_xacc(ind1,ind2) = 1;
    %allXad(ind1,ind2) = 0.15;
    
    for i = 1:8
            k = 1;
            for j = 1:6
            newX(j,:,i) = allXad(i,k+1:k+5);
            newY(j,:,i) = allYad(i,k+1:k+5);
            newtheta(j,:,i) = alltheta(i,k+1:k+5);
            newCL3(j,:,i) = CL3(i,k+1:k+5);%./sqrt(2);
            newCL1(j,:,i) = CL1(i,k+1:k+5);%./sqrt(2);
            newCD2(j,:,i) = CD2(i,k+1:k+5);%./sqrt(2);
            newfft_posy(j,:,i) = fft_ypos(i,k+1:k+5);
%             if isnan(fft_xacc(i,k)) == 1
%                 newfft_accx(j,:,i) = [0;0;0;0;0;0];
%                 newfft_velx(j,:,i) = [0;0;0;0;0;0];
%                 newfft_posx(j,:,i) = [0;0;0;0;0;0];
%             else
                newfft_accx(j,:,i) = fft_xacc(i,k+1:k+5);
                newfft_velx(j,:,i) = fft_xvel(i,k+1:k+5);
                newfft_posx(j,:,i) = fft_xpos(i,k+1:k+5);
%             end
            newfft_accy(j,:,i) = fft_yacc(i,k+1:k+5);
            newfft_vely(j,:,i) = fft_yvel(i,k+1:k+5);
            newfft_cdf(j,:,i) = fft_cdf(i,k+1:k+5);
            newfft_cl1(j,:,i) = fft_cl(i,k+1:k+5);
            newfft_cl3(j,:,i) = fft_cl3(i,k+1:k+5);
            vr(j,:,i) = Vr(i,k+1:k+5);
            
            if p == 1 | p == 2 | p == 3
                avgPow(j,:,i) = (Plift(i,k+1:k+5)+Pdrag(i,k+1:k+5))./(0.5*1000*0.18.^3*0.6858*0.0381);
            else
                avgPow(j,:,i) = (Plift(i,k+1:k+5)+Pdrag(i,k+1:k+5))./(0.5*1000*0.23.^3*0.6858*0.0381);
            end
            k = k+6;
            end
    end
    newX(:,:,9) = newX(:,:,1);
    newY(:,:,9) = newY(:,:,1);
    newtheta(:,:,9) = 180;
    newCL3(:,:,9) = newCL3(:,:,1);
    newCL1(:,:,9) = newCL1(:,:,1);
    newCD2(:,:,9) = newCD2(:,:,1);
    vr(:,:,9) = vr(:,:,1);
    newfft_posy(:,:,9) = newfft_posy(:,:,1);
    newfft_accx(:,:,9) = newfft_accx(:,:,1);
    newfft_velx(:,:,9) = newfft_velx(:,:,1);
    newfft_posx(:,:,9) = newfft_posx(:,:,1);
    newfft_accy(:,:,9) = newfft_accy(:,:,1);
    newfft_vely(:,:,9) = newfft_vely(:,:,1);
    newfft_cdf(:,:,9) = newfft_cdf(:,:,1);
    newfft_cl1(:,:,9) = newfft_cl1(:,:,1);
    newfft_cl3(:,:,9) = newfft_cl3(:,:,1);
    avgPow(:,:,9) = avgPow(:,:,1);
    
    bigX(:,:,:,p) = newX;
    bigY(:,:,:,p) = newY;
    bigtheta(:,:,:,p) = newtheta;
    bigCL1(:,:,:,p) = newCL1;
    bigCL3(:,:,:,p) = newCL3;
    bigCD2(:,:,:,p) = newCD2;
    bigvr(:,:,:,p) = vr;
    bigfft_posy(:,:,:,p) = newfft_posy;
    bigfft_vely(:,:,:,p) = newfft_vely;
    bigfft_accy(:,:,:,p) = newfft_accy;
    bigfft_posx(:,:,:,p) = newfft_posx;
    bigfft_velx(:,:,:,p) = newfft_velx;
    bigfft_accx(:,:,:,p) = newfft_accx;
    bigfft_cdf(:,:,:,p) = newfft_cdf;
    bigfft_cl1(:,:,:,p) = newfft_cl1;
    bigfft_cl3(:,:,:,p) = newfft_cl3;
    bigpow(:,:,:,p) = avgPow;
end

bigXi = interpn(bigX,ifac,int);
bigYi = interpn(bigY,ifac,int);
bigthetai = interpn(bigtheta,ifac,int);
bigCL1i = interpn(bigCL1,ifac,int);
bigCL3i = interpn(bigCL3,ifac,int);
bigCD2i = interpn(bigCD2,ifac,int);
bigvri = interpn(bigvr,ifac,int);
bigfft_posy_ri = interpn(real(bigfft_posy),ifac,int);
bigfft_vely_ri = interpn(real(bigfft_vely),ifac,int);
bigfft_accy_ri = interpn(real(bigfft_accy),ifac,int);
bigfft_posx_ri = interpn(real(bigfft_posx),ifac,int);
bigfft_velx_ri = interpn(real(bigfft_velx),ifac,int);
bigfft_accx_ri = interpn(real(bigfft_accx),ifac,int);
bigfft_cdf_ri = interpn(real(bigfft_cdf),ifac,int);
bigfft_cl1_ri = interpn(real(bigfft_cl1),ifac,int);
bigfft_cl3_ri = interpn(real(bigfft_cl3),ifac,int);
bigfft_posy_ii = interpn(imag(bigfft_posy),ifac,int);
bigfft_vely_ii = interpn(imag(bigfft_vely),ifac,int);
bigfft_accy_ii = interpn(imag(bigfft_accy),ifac,int);
bigfft_posx_ii = interpn(imag(bigfft_posx),ifac,int);
bigfft_velx_ii = interpn(imag(bigfft_velx),ifac,int);
bigfft_accx_ii = interpn(imag(bigfft_accx),ifac,int);
bigfft_cdf_ii = interpn(imag(bigfft_cdf),ifac,int);
bigfft_cl1_ii = interpn(imag(bigfft_cl1),ifac,int);
bigfft_cl3_ii = interpn(imag(bigfft_cl3),ifac,int);
bigpowi = interpn(bigpow,ifac,int);

bigfft_posy_ri = bigfft_posy_ri./hypot(bigfft_posy_ri,bigfft_posy_ii);
bigfft_vely_ri = bigfft_vely_ri./hypot(bigfft_vely_ri,bigfft_vely_ii);
bigfft_accy_ri = bigfft_accy_ri./hypot(bigfft_accy_ri,bigfft_accy_ii);
bigfft_posx_ri = bigfft_posx_ri./hypot(bigfft_posx_ri,bigfft_posx_ii);
bigfft_velx_ri = bigfft_velx_ri./hypot(bigfft_velx_ri,bigfft_velx_ii);
bigfft_accx_ri = bigfft_accx_ri./hypot(bigfft_accx_ri,bigfft_accx_ii);
bigfft_cdf_ri = bigfft_cdf_ri./hypot(bigfft_cdf_ri,bigfft_cdf_ii);
bigfft_cl1_ri = bigfft_cl1_ri./hypot(bigfft_cl1_ri,bigfft_cl1_ii);
bigfft_cl3_ri = bigfft_cl3_ri./hypot(bigfft_cl3_ri,bigfft_cl3_ii);
bigfft_posy_ii = bigfft_posy_ii./hypot(bigfft_posy_ri,bigfft_posy_ii);
bigfft_vely_ii = bigfft_vely_ii./hypot(bigfft_vely_ri,bigfft_vely_ii);
bigfft_accy_ii = bigfft_accy_ii./hypot(bigfft_accy_ri,bigfft_accy_ii);
bigfft_posx_ii = bigfft_posx_ii./hypot(bigfft_posx_ri,bigfft_posx_ii);
bigfft_velx_ii = bigfft_velx_ii./hypot(bigfft_velx_ri,bigfft_velx_ii);
bigfft_accx_ii = bigfft_accx_ii./hypot(bigfft_accx_ri,bigfft_accx_ii);
bigfft_cdf_ii = bigfft_cdf_ii./hypot(bigfft_cdf_ri,bigfft_cdf_ii);
bigfft_cl1_ii = bigfft_cl1_ii./hypot(bigfft_cl1_ri,bigfft_cl1_ii);
bigfft_cl3_ii = bigfft_cl3_ii./hypot(bigfft_cl3_ri,bigfft_cl3_ii);

vr_steps = 4.5:0.5/(2^ifac):8;
weights = fliplr(vr_steps);

freefile = '1p9';

eval(['load new' freefile]);

Vr1p9 = Vrn*fy./(ypeakfreq/2/pi);
linx = Vr1p9(9:24);
liny = vels(9:24)';
slope = (linx'*linx)^-1*linx'*liny;
lino = ypeakfreq(9:24)';
slope2 = (liny'*liny)^-1*liny'*lino';

num=1;



for i = 1:length(vr_steps)

    %U = 0.349;
    %U = 0.5;
    D = 0.0762;
    S = 2;
    %mx = 45.9;
    %my = 51.8;
    %kx = 3235;
    %kx = 4000;
    %ky = 1013;
    %zetax = 0.0620;
    %zetay = 0.0250;
    %fny = 0.704;
    %fnx = 1.336;
    by = zetay*2*my*fy*2*pi;
    bx = zetax*2*mx*fx*2*pi;
    %by =0;
    %bx=0;
    %vr_guess = 7;
    U(i) = vr_steps(i)*slope;
    omega(i) = U(i)*slope2;
    %omega(i) =  2*pi*U(i)/(vr_steps(i)*D);
%     omega = 2*pi*fny;
%     U = vr_steps(i)*D*fny;


%     xr = (mx*4*omega.^2*D*bigXi.*bigfft_accx_ri + bx*omega*D*bigXi.*bigfft_velx_ii + kx*D*bigXi.*bigfft_posx_ri)./(1/2*1000*U^2*D*S);
%     xi = (mx*4*omega.^2*D*bigXi.*bigfft_accx_ii + bx*omega*D*bigXi.*bigfft_velx_ri + kx*D*bigXi.*bigfft_posx_ii)./(1/2*1000*U^2*D*S);
%     yr = (my*omega.^2*D*bigYi.*bigfft_accy_ri + by*omega*D*bigYi.*bigfft_vely_ii + ky*D*bigYi.*bigfft_posy_ri)./(1/2*1000*U^2*D*S);
%     yi = (my*omega.^2*D*bigYi.*bigfft_accy_ii + by*omega*D*bigYi.*bigfft_vely_ri + ky*D*bigYi.*bigfft_posy_ii)./(1/2*1000*U^2*D*S);
    xr = (mx*4*omega(i).^2*D*bigXi.*bigfft_accx_ri + bx*omega(i)*D*bigXi.*bigfft_velx_ri + kx*D*bigXi.*bigfft_posx_ri);%./(1/2*1000*U(i)^2*D*S);
    xi = (mx*4*omega(i).^2*D*bigXi.*bigfft_accx_ii + bx*omega(i)*D*bigXi.*bigfft_velx_ii + kx*D*bigXi.*bigfft_posx_ii);%./(1/2*1000*U(i)^2*D*S);
    yr = (my*omega(i).^2*D*bigYi.*bigfft_accy_ri + by*omega(i)*D*bigYi.*bigfft_vely_ri + ky*D*bigYi.*bigfft_posy_ri);%./(1/2*1000*U(i)^2*D*S);
    yi = (my*omega(i).^2*D*bigYi.*bigfft_accy_ii + by*omega(i)*D*bigYi.*bigfft_vely_ii + ky*D*bigYi.*bigfft_posy_ii);%./(1/2*1000*U(i)^2*D*S);
%     xr = real(mx*4*omega^2*D*bigXi.*(bigfft_posx_ri+sqrt(-1)*bigfft_posx_ii) + ...
%         sqrt(-1)*bx*omega*D*bigXi.*(bigfft_posx_ri+sqrt(-1)*bigfft_posx_ii) + ...
%         kx*D*bigXi.*bigfft_posx_ri)./(1/2*1000*U^2*D*S);
%     xi = imag(mx*4*omega^2*D*bigXi.*(bigfft_posx_ri+sqrt(-1)*bigfft_posx_ii) + ...
%         sqrt(-1)*bx*omega*D*bigXi.*(bigfft_posx_ri+sqrt(-1)*bigfft_posx_ii) + ...
%         kx*D*bigXi.*bigfft_posx_ri)./(1/2*1000*U^2*D*S);
%     yr = real(my*omega^2*D*bigYi.*(bigfft_posy_ri+sqrt(-1)*bigfft_posy_ii) + ...
%         sqrt(-1)*by*omega*D*bigYi.*(bigfft_posy_ri+sqrt(-1)*bigfft_posy_ii) + ...
%         ky*D*bigYi.*bigfft_posy_ri)./(1/2*1000*U^2*D*S);
%     yi = imag(my*omega^2*D*bigYi.*(bigfft_posy_ri+sqrt(-1)*bigfft_posy_ii) + ...
%         sqrt(-1)*by*omega*D*bigYi.*(bigfft_posy_ri+sqrt(-1)*bigfft_posy_ii) + ...
%         ky*D*bigYi.*bigfft_posy_ri)./(1/2*1000*U^2*D*S);
%     xr = (-mx*4*omega^2*D*bigXi.*bigfft_accx_ri - bx*D*omega*bigXi.*bigfft_velx_ii + kx*D*bigXi.*bigfft_posx_ri)./(1/2*1000*U^2*D*S);
%     xi = (-mx*4*omega^2*D*bigXi.*bigfft_accx_ii + bx*D*omega*bigXi.*bigfft_velx_ri + kx*D*bigXi.*bigfft_posx_ii)./(1/2*1000*U^2*D*S);
%     yr = (-my*omega^2*D*bigYi.*bigfft_accy_ri - by*D*omega*bigYi.*bigfft_vely_ii + ky*D*bigYi.*bigfft_posy_ri)./(1/2*1000*U^2*D*S);
%     yi = (-my*omega^2*D*bigYi.*bigfft_accy_ii + by*D*omega*bigYi.*bigfft_vely_ri + ky*D*bigYi.*bigfft_posy_ii)./(1/2*1000*U^2*D*S);
    cdr = bigfft_cdf_ri.*bigCD2i*1/2*1000*U(i)^2*D*S;
    cdi = bigfft_cdf_ii.*bigCD2i*1/2*1000*U(i)^2*D*S;
    cyr = bigfft_cl1_ri.*bigCL1i*1/2*1000*U(i)^2*D*S;% + bigfft_cl3_ri.*bigCL3i;
    cyi = bigfft_cl1_ii.*bigCL1i*1/2*1000*U(i)^2*D*S;% + bigfft_cl3_ii.*bigCL3i;

    dvr = abs(bigvri-vr_steps(i));
    indices = find(dvr < 1*0.5/(ifac^2));
    [p1,p2,p3,p4] = ind2sub(size(bigXi),indices);
    d = 10000*ones(size(bigXi));
%     d1 = 100*ones(size(bigXi));
%     d2 = 100*ones(size(bigXi));
%     d3 = 100*ones(size(bigXi));
%     d4 = 100*ones(size(bigXi));
    
    maxxr = max([abs(max(max(max(max(abs(xr)))))) abs(min(min(min(min(xr)))))]);
    maxxi = max([abs(max(max(max(max(abs(xi)))))) abs(min(min(min(min(xi)))))]);
    maxyr = max([abs(max(max(max(max(abs(yr)))))) abs(min(min(min(min(yr)))))]);
    maxyi = max([abs(max(max(max(max(abs(yi)))))) abs(min(min(min(min(yi)))))]);
    
    for z = 1:length(p1)
         %if vr_steps(i) < 6.51
            d(p1(z),p2(z),p3(z),p4(z)) = abs(xr(p1(z),p2(z),p3(z),p4(z))-cdr(p1(z),p2(z),p3(z),p4(z))) + ...%/abs(sqrt(xr(p1(z),p2(z),p3(z),p4(z)).^2 + xi(p1(z),p2(z),p3(z),p4(z)).^2)) + ...
                abs(xi(p1(z),p2(z),p3(z),p4(z))-cdi(p1(z),p2(z),p3(z),p4(z))) + .../abs(sqrt(xr(p1(z),p2(z),p3(z),p4(z)).^2 + xi(p1(z),p2(z),p3(z),p4(z)).^2)) + ...
                abs(yr(p1(z),p2(z),p3(z),p4(z))-cyr(p1(z),p2(z),p3(z),p4(z))) + .../abs(sqrt(yr(p1(z),p2(z),p3(z),p4(z)).^2 + yi(p1(z),p2(z),p3(z),p4(z)).^2)) + ...
                abs(yi(p1(z),p2(z),p3(z),p4(z))-cyi(p1(z),p2(z),p3(z),p4(z)));%/abs(sqrt(yr(p1(z),p2(z),p3(z),p4(z)).^2 + yi(p1(z),p2(z),p3(z),p4(z)).^2));%+ abs(bigpowi(p1(z),p2(z),p3(z),p4(z)));
         %elseif vr_steps(i) > 6.51
         %   d(p1(z),p2(z),p3(z),p4(z)) = 2*abs(xr(p1(z),p2(z),p3(z),p4(z))-cdr(p1(z),p2(z),p3(z),p4(z)))/abs(sqrt(xr(p1(z),p2(z),p3(z),p4(z)).^2 + xi(p1(z),p2(z),p3(z),p4(z)).^2)) + ...
         %       2*abs(xi(p1(z),p2(z),p3(z),p4(z))-cdi(p1(z),p2(z),p3(z),p4(z)))/abs(sqrt(xr(p1(z),p2(z),p3(z),p4(z)).^2 + xi(p1(z),p2(z),p3(z),p4(z)).^2)) + ...
         %       abs(yr(p1(z),p2(z),p3(z),p4(z))-cyr(p1(z),p2(z),p3(z),p4(z)))/abs(sqrt(yr(p1(z),p2(z),p3(z),p4(z)).^2 + yi(p1(z),p2(z),p3(z),p4(z)).^2)) + ...
         %       3*abs(yi(p1(z),p2(z),p3(z),p4(z))-cyi(p1(z),p2(z),p3(z),p4(z)))/abs(sqrt(yr(p1(z),p2(z),p3(z),p4(z)).^2 + yi(p1(z),p2(z),p3(z),p4(z)).^2));% + abs(bigpowi(p1(z),p2(z),p3(z),p4(z)));
         %end
            %         elseif vr_steps(i) > 5.49 & vr_steps(i) < 6.01
%             d(p1(z),p2(z),p3(z),p4(z)) = 2*abs(xr(p1(z),p2(z),p3(z),p4(z))-cdr(p1(z),p2(z),p3(z),p4(z))) + ...
%                 5*abs(xi(p1(z),p2(z),p3(z),p4(z))-cdi(p1(z),p2(z),p3(z),p4(z))) + ...
%                 2*abs(yr(p1(z),p2(z),p3(z),p4(z))-cyr(p1(z),p2(z),p3(z),p4(z))) + ...
%                 abs(yi(p1(z),p2(z),p3(z),p4(z))-cyi(p1(z),p2(z),p3(z),p4(z)));
%         elseif vr_steps(i) > 6.01
%             d(p1(z),p2(z),p3(z),p4(z)) = abs(xr(p1(z),p2(z),p3(z),p4(z))-cdr(p1(z),p2(z),p3(z),p4(z))) + ...
%                 abs(xi(p1(z),p2(z),p3(z),p4(z))-cdi(p1(z),p2(z),p3(z),p4(z))) + ...
%                 abs(yr(p1(z),p2(z),p3(z),p4(z))-cyr(p1(z),p2(z),p3(z),p4(z))) + ...
%                 abs(yi(p1(z),p2(z),p3(z),p4(z))-cyi(p1(z),p2(z),p3(z),p4(z)));
%         end
%             d(p1(z),p2(z),p3(z),p4(z)) = 5*abs((xr(p1(z),p2(z),p3(z),p4(z))-cdr(p1(z),p2(z),p3(z),p4(z)))./...
%                 (abs(complex(xr(p1(z),p2(z),p3(z),p4(z))-cdr(p1(z),p2(z),p3(z),p4(z)),...
%                 xi(p1(z),p2(z),p3(z),p4(z))-cdi(p1(z),p2(z),p3(z),p4(z)))))) + ...
%             5*abs((xi(p1(z),p2(z),p3(z),p4(z))-cdi(p1(z),p2(z),p3(z),p4(z)))./(abs(complex(xr(p1(z),p2(z),p3(z),p4(z))-cdr(p1(z),p2(z),p3(z),p4(z)),xi(p1(z),p2(z),p3(z),p4(z))-cdi(p1(z),p2(z),p3(z),p4(z)))))) + ...
%             abs((yr(p1(z),p2(z),p3(z),p4(z))-cyr(p1(z),p2(z),p3(z),p4(z)))./(abs(complex(yr(p1(z),p2(z),p3(z),p4(z))-cyr(p1(z),p2(z),p3(z),p4(z)),yi(p1(z),p2(z),p3(z),p4(z))-cyi(p1(z),p2(z),p3(z),p4(z)))))) + ...
%             abs((yi(p1(z),p2(z),p3(z),p4(z))-cyi(p1(z),p2(z),p3(z),p4(z)))./(abs(complex(yr(p1(z),p2(z),p3(z),p4(z))-cyr(p1(z),p2(z),p3(z),p4(z)),yi(p1(z),p2(z),p3(z),p4(z))-cyi(p1(z),p2(z),p3(z),p4(z))))));
    
    end
    
%     for z = 1:length(p1)
%         d1(p1(z),p2(z),p3(z),p4(z)) = abs(yr(p1(z),p2(z),p3(z),p4(z))-cyr(p1(z),p2(z),p3(z),p4(z)))/abs(yr(p1(z),p2(z),p3(z),p4(z)));
%         d2(p1(z),p2(z),p3(z),p4(z)) = abs(xr(p1(z),p2(z),p3(z),p4(z))-cdr(p1(z),p2(z),p3(z),p4(z)))/abs(xr(p1(z),p2(z),p3(z),p4(z)));
%         d3(p1(z),p2(z),p3(z),p4(z)) = abs(xi(p1(z),p2(z),p3(z),p4(z))-cdi(p1(z),p2(z),p3(z),p4(z)))/abs(xi(p1(z),p2(z),p3(z),p4(z)));
%         d4(p1(z),p2(z),p3(z),p4(z)) = abs(yi(p1(z),p2(z),p3(z),p4(z))-cyi(p1(z),p2(z),p3(z),p4(z)))/abs(yi(p1(z),p2(z),p3(z),p4(z)));
%     end
%     
%     index = find(d1 < 2 & d2 < 0.5 & d3 < 0.5 & d4 < 2);
%d = abs(xr-cdr) + abs(xi-cdi) + abs(yr-cyr) + abs(yi-cyi);
%d = abs(yr-cyr) + abs(yi-cyi);
%d1 = abs(xr-cdr);
%d2 = abs(xi-cdi);
%d3 = abs(yr-cyr);
%d4 = abs(yi-cyi);

%num = 100;
 index = find(d == min(min(min(min(d)))));
%index = find(d < min(min(min(min(d))))*1.05);
%index = find(d1 < num & d2 < num & d3 < num & d4 < num);
[i1,i2,i3,i4] = ind2sub(size(bigXi),index);

for h = 1:length(index)
%     Xpred(i) = bigXi(i1,i2,i3,i4);
%     Ypred(i) = bigYi(i1,i2,i3,i4);
%     Thetapred(i) = bigthetai(i1,i2,i3,i4);
%     Vrpred(i) = bigvri(i1,i2,i3,i4);
%     Powpred(i) = bigpowi(i1,i2,i3,i4);
    Xpred(num) = bigXi(i1(h),i2(h),i3(h),i4(h));
    Ypred(num) = bigYi(i1(h),i2(h),i3(h),i4(h));
    Thetapred(num) = bigthetai(i1(h),i2(h),i3(h),i4(h));
    Vrpred(num) = bigvri(i1(h),i2(h),i3(h),i4(h));
    Powpred(num) = bigpowi(i1(h),i2(h),i3(h),i4(h));
    num=num+1;
end



end

eval(['load new' freefile]);

Vr1p9 = Vrn*fy./(ypeakfreq/2/pi);

figure(1);
plot(Vr1p9(9:24),Yad(9:24),'ko',vr_steps,Ypred,'rs');
xlabel('Vr');
ylabel('Y/D');
legend('Observed','Predicted');

%print(gcf,'-depsc','Y_amp_prediction_1p9');

figure(2);
plot(Vr1p9(9:24),Xad(9:24),'ko',vr_steps,Xpred,'rs');
xlabel('Vr');
ylabel('X/D');
legend('Observed','Predicted');

%print(gcf,'-depsc','X_amp_prediction_1p9');

eval(['load phase' freefile]);

for i = 1:length(phasextoy)
    if phasextoy(i) > 180
        phasextoy(i) = phasextoy(i)-360;
    else
        phasextoy(i) = phasextoy(i);
    end
end

figure(3);
plot(Vr1p9(9:24),phasextoy(9:24),'ko',vr_steps,Thetapred,'rs');
xlabel('Vr');
ylabel('\theta');
legend('Observed','Predicted');

%print(gcf,'-depsc','Theta_prediction_1p9');

% for i = 1:29
% figure(4);clf;hold on;
% % h = patch(isosurface(bigXi(:,:,:,i),bigYi(:,:,:,i),bigthetai(:,:,:,i),d(:,:,:,i),min(min(min(d(:,:,:,i))))*1.1));
% % set(h,'FaceColor','g');
% % set(h,'FaceAlpha',0.5)
% xlabel('X/D','FontSize',16);
% ylabel('Y/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% title('D2','FontSize',16);
% view(-15,20);
% axis([0 1 0 1.5 -180 180]);
% grid on
% pause
% end


% figure(4);clf;hold on;
% h = patch(isosurface(newXi,newYi,newthetai,newPhasei,0));
% set(h,'FaceColor','g');
% set(h,'FaceAlpha',0.5)
% h2 = patch(isosurface(newXi,newYi,newthetai,newPhasei,pi));
% set(h2,'FaceColor','r');
% set(h2,'FaceAlpha',0.5)
% h3 = patch(isosurface(newXi,newYi,newthetai,newPhasei,-pi));
% set(h3,'FaceColor','b');
% set(h3,'FaceAlpha',0.5)
% xlabel('X/D','FontSize',16);
% ylabel('Y/D','FontSize',16);
% zlabel('\theta','FontSize',16);
% title('D2','FontSize',16);
% view(-15,20);
% grid on