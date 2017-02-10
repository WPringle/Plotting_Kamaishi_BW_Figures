% Script to estimate the bottom shear stress in the 2CLOWNS calculation
clearvars;
clc;
close all;

% Coefficients:
rho = 1025; % density of salt water
ks10 = 0.2; %0.60;  %D50 on rubble mound (10 kg or 300 kg)
ks0  = 1d-3;  %D50 on seabed;
man = 0.025; %mannings
rho_s = 2650; %density of sediment
karman = 0.41;
g = 9.807; %gravity
nu = 1d-6; % viscosity of water

%% 3D part
direc = '../Satakev8.0_3D_Real_BW/'; 
fname = 'Satake_BW_Real_ke.'; 
load([direc 'MAT_files/' fname 'xyznF.mat'])

% Get the nn cell numbers for bottom cell
fb_b = zeros(ie(2)-is(2),ie(1)-is(1)); 
nnip_b = zeros(ie(2)-is(2),ie(1)-is(1)); 
nnjp_b = zeros(ie(2)-is(2),ie(1)-is(1)); 
nnim_b = zeros(ie(2)-is(2),ie(1)-is(1)); 
nnjm_b = zeros(ie(2)-is(2),ie(1)-is(1)); 
nfb_b = zeros(ie(2)-is(2),ie(1)-is(1)); 
ii = 0; jj = 0;
for i = is(1):ie(1)-1
    ii = ii + 1; jj = 0;
    for j = is(2):ie(2)-1
        jj = jj + 1;
        for k = is(3):ie(3)-1
            nf_n = nf(j,k,i);
            if nf_n > -1
                nfb_b(jj,ii) = nfb(j,k-1,i);
                fb_b(jj,ii) = 0;
                nfim_n = nf(j,k,i-1); nfbim_n = nfb(j,k,i-1);
                if nfim_n > -1 || (nfbim_n > 0 && nfbim_n < 10)
                    nnim_b(jj,ii) = mn(j,k,i);
                    fb_b(jj,ii) = fb_b(jj,ii) + 0.25*a(mn(j,k,i),1);
                else
                    nnim_b(jj,ii) = mn(j,k+1,i);
                    fb_b(jj,ii) = fb_b(jj,ii) + 0.25*a(mn(j,k+1,i),1);
                end
                nfjm_n = nf(j-1,k,i); nfbjm_n = nfb(j-1,k,i);
                if nfjm_n > -1 || (nfbjm_n > 0 && nfbjm_n < 10)
                    nnjm_b(jj,ii) = mn(j,k,i);
                    fb_b(jj,ii) = fb_b(jj,ii) + 0.25*a(mn(j,k,i),2);
                else
                    nnjm_b(jj,ii) = mn(j,k+1,i);
                    fb_b(jj,ii) = fb_b(jj,ii) + 0.25*a(mn(j,k+1,i),2);
                end
                nfip_n = nf(j,k,i+1); nfbip_n = nfb(j,k,i+1);
                if nfip_n > -1 || (nfbip_n > 0 && nfbip_n < 10)
                    nnip_b(jj,ii) = mn(j,k,i+1);
                    fb_b(jj,ii) = fb_b(jj,ii) + 0.25*a(mn(j,k,i+1),1);
                else
                    nnip_b(jj,ii) = mn(j,k+1,i+1);
                    fb_b(jj,ii) = fb_b(jj,ii) + 0.25*a(mn(j,k+1,i+1),1);
                end
                nfjp_n = nf(j+1,k,i); nfbjp_n = nfb(j+1,k,i);
                if nfjp_n > -1 || (nfbjp_n > 0 && nfbjp_n < 10)
                    nnjp_b(jj,ii) = mn(j+1,k,i);
                    fb_b(jj,ii) = fb_b(jj,ii) + 0.25*a(mn(j+1,k,i),2);
                else
                    nnjp_b(jj,ii) = mn(j+1,k+1,i);
                    fb_b(jj,ii) = fb_b(jj,ii) + 0.25*a(mn(j+1,k+1,i),2);
                end
                break
            end
        end
    end
end

% load([direc 'MAT_files/' fname 'VelF.mat'])
% 
% Velim = Vel3(:,nnim_b,1);
% Veljm = Vel3(:,nnjm_b,2);
% Velip = Vel3(:,nnip_b,1);
% Veljp = Vel3(:,nnjp_b,2);
% Velim = reshape(Velim,length(time),size(nfb_b,1),size(nfb_b,2));
% Veljm = reshape(Veljm,length(time),size(nfb_b,1),size(nfb_b,2));
% Velip = reshape(Velip,length(time),size(nfb_b,1),size(nfb_b,2));
% Veljp = reshape(Veljp,length(time),size(nfb_b,1),size(nfb_b,2));
% U = 0.5*(Velim(:,:,:) + Velip(:,:,:));
% V = 0.5*(Veljm(:,:,:) + Veljp(:,:,:));
% BVel = sqrt(U.^2 + V.^2);
% % 
% save([direc 'MAT_files/' fname '_BVel.mat'],'BVel')
load([direc 'MAT_files/' fname '_BVel.mat'])

% 
BVelm = squeeze(max(BVel));
%BStd  = squeeze(std(BVel,0,1));
%BVel98  = squeeze(prctile(BVel,98));
% yp = 0.5*2.2*fb_b;
% yk(nfb_b == 0) = ks0*exp(-8*karman);
% yk(nfb_b == 10) = ks10*exp(-8*karman);
% yk = reshape(yk,size(nfb_b,1),size(nfb_b,2));
% Ush = karman * BVel ./ log(yp./yk); 
% Tb = rho * Ush.^2;    
Cf = 0.005;
Tb = 0.5 * rho * Cf .* BVelm.^2;
Tbs(nfb_b == 0) = Tb(nfb_b == 0) / ((rho_s - rho)*g*ks0); 
Tbs(nfb_b == 10) = Tb(nfb_b == 10) / ((rho_s - rho)*g*ks10);  
Tbs = reshape(Tbs,size(nfb_b,1),size(nfb_b,2));

% 
% plot the color plot
figure;
pcolor(x(is(1):ie(1)-1),y(is(2):ie(2)-1),Tbs); 
shading interp
colormap jet
caxis([0 10])

%
xx = 0.5*(x(is(1):ie(1)-1) + x(is(1)+1:ie(1)));
yy = 0.5*(y(is(2):ie(2)-1) + y(is(1)+1:ie(2)));
xx = xx * 1d-6; yy = yy * 1d-6 + 3.5;
grdwrite2(xx,yy,Tbs,[direc fname 'Tbs3D.nc'])

%% 2DH part
for ii = 1:2
    if ii == 1
        direc = '../Satakev8.0_3D_Real_BW/'; 
        fname = 'Satakev8.0_3D_Real_BW_';  
    elseif ii == 2
        direc = '../Satakev8.0_2DH_BW_FD/'; 
        fname = 'Satakev8.0_2DH_BW_FD'; 
    end
    load([direc 'MAT_files/' fname 'xyzF.mat'])
    % 
    is = find ( x >= 406080 & x < 409430);
    js = find ( y >= 844700 & y < 845630);
    x = x(is); y = y(js);
    % 
    % load([direc 'MAT_files/' fname 'VelEta.mat'])
    % 
    % % reorder Eta, Vel and get Depth  
    % [time2D,tI] = sort(time); time2D = time2D';
    % Eta = Eta(tI,:); Vel = Vel(tI,:,:); 
    % H = Eta - repmat(Zk,1,size(Eta,1))';
    % clear Eta
    % 
    % % find the nn values inside domain
    % nn = squeeze(mn(js,is)); nnu = squeeze(mn(js(2):js(end)+1,is(2):is(end)+1));
    % 
    % % get depth
    % H = H(:,nn); H(H <= 0) = NaN; 
    % H = reshape(H,length(time2D),size(nn,1),size(nn,2));
    % 
    % % get cell center velocities
    % Vel1 = Vel(:,nn,:); Vel2 = Vel(:,nnu,:);
    % Vel1 = reshape(Vel1,length(time2D),size(nn,1),size(nn,2),2);
    % Vel2 = reshape(Vel2,length(time2D),size(nnu,1),size(nnu,2),2);
    % U = 0.5*(Vel1(:,:,:,1) + Vel2(:,:,:,1));
    % V = 0.5*(Vel1(:,:,:,2) + Vel2(:,:,:,2));
    % Vel2 = U.^2 + V.^2;
    % 
    % save([direc 'MAT_files/' fname '_HVel2.mat'],'Vel2','H')
    load([direc 'MAT_files/' fname '_HVel2.mat'])

    % Get friction coefficient
%     fS(:,nfb_b == 0)  = 2*log10(H(:,nfb_b == 0)/ks0)+1.14; %30/(sqrt(8)*0.4)*(log(30*H(:,nfb_b == 0)/ks0)-1);  % 
%     fS(:,nfb_b == 10) = 2*log10(H(:,nfb_b == 10)/ks10)+1.14;
%     fS = reshape(fS,size(fS,1),size(nfb_b,1),size(nfb_b,2));
%     f  = (1./fS).^2; 
%     Cf = 0.25 * f;
    %
    %Vel = squeeze(prctile(sqrt(Vel2),95));
    Vel2 = squeeze(max(Vel2));
    
    %Cf = 2*g*man^2./H.^(1/3);
    Cf = 0.005;

    % Now convert to bottom shear stresses...
    Tb2D = 0.5 * rho * Cf .* Vel2;

    % get maximum stress
    % dimensionless
    Tb2Ds(nfb_b == 0) = Tb2D(nfb_b == 0) / ((rho_s - rho)*g*ks0); 
    Tb2Ds(nfb_b == 10) = Tb2D(nfb_b == 10) / ((rho_s - rho)*g*ks10); 
    Tb2Ds = reshape(Tb2Ds,size(nfb_b,1),size(nfb_b,2));

    % plot the color plot
    figure;
    pcolor(x,y,Tb2Ds); 
    shading interp
    colormap jet
    caxis([0 10])
    
    %grdwrite2(xx,yy,Tb2Ds,[direc fname 'Tbs2D.nc'])
end
