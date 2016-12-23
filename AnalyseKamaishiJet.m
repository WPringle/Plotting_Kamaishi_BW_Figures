% Script to analyse the jet through the submerged section of the breakwater
% as computed by 2CLOWNS-3D simulation. 
% By: William Pringle Sep 20 2016

clearvars; clc; close all;
%% Set and load the data
direc = 'Satakev8.0_Real_ke_HR/';
fname = 'SatakeTNorm_HYB_realizable_ke_HR.';
oname = 'Real_ke';
% direc = 'Satakev8.0_ke_HR_DeepOn/';
% fname = 'Tohoku_SatakeT_HYB_ke_k1mm_DeepOn_HR.';
% oname = 'DeepOn';

load([direc 'MAT_files/' fname 'xyzn.mat']);  
load([direc 'MAT_files/' fname 'time_av2.mat'])
%load([direc 'MAT_files/' fname 'pro.mat'])
%turb_pro = mean(turb_pro)';
%load([direc fname time '.mat'])

%% Choose plane and wanted cross-section
plane = 'x-z';
%z_want = -8;
z_want  = 845140;
%z_want  = 407100;

%% Make desired target into grid format for desired plane
switch plane
    case  'x-y'
        % Make desired target into grid format for x-y plane
        xx = 0.5*(x(is(1):ie(1)-1) + x(is(1)+1:ie(1)));
        yy = 0.5*(y(is(2):ie(2)-1) + y(is(2)+1:ie(2))) + 3.5d6;
        xx = xx/1d6;
        yy = yy/1d6;
        dz = 2.2; %z(is(3)+1)-z(is(3));
        js = is(2); je = ie(2)-1;
        iss = is(1) ; iee = ie(1)-1;
        % Get cross section
        k  = find( z > z_want - dz/2 & z < z_want + dz/2);
        % For getting correct indices
        iy = 0; ix = 1; iz = 0;
        jy = 1; jx = 0; jz = 0;
        u1i = 1; u2i = 2;
    case  'x-z'    
        % Make desired target into grid format for x-y plane
        xx = 0.5*(x(is(1):ie(1)-1) + x(is(1)+1:ie(1)));
        yy = 0.5*(z(is(3):ie(3)-1) + z(is(3)+1:ie(3)));
        xx = xx/1d6;
        dy = y(is(2)+1)-y(is(2));
        js = is(3); je = ie(3)-1;
        iss = is(1) ; iee = ie(1)-1;
        % Get cross section
        j  = find( y > z_want - dy/2 & y < z_want + dy/2);
        % For getting correct indices
        iy = 0; ix = 1; iz = 0;
        jy = 0; jx = 0; jz = 1;
        u1i = 1; u2i = 3;
    case  'y-z'
        % Make desired target into grid format for x-y plane
        xx = 0.5*(y(is(2):ie(2)-1) + y(is(2)+1:ie(2))) + 3.5d6;
        yy = 0.5*(z(is(3):ie(3)-1) + z(is(3)+1:ie(3)));
        xx = xx/1d6;
        dx = x(is(1)+1)-x(is(1));
        js = is(3); je = ie(3)-1;
        iss = is(2) ; iee = ie(2)-1;
        % Get cross section
        i  = find( x > z_want - dx/2 & x < z_want + dx/2);
        % For getting correct indices
        iy = 1; ix = 0; iz = 0;
        jy = 0; jx = 0; jz = 1;
        u1i = 2; u2i = 3;
end
%
uv = zeros(length(yy),length(xx));
ux = zeros(length(yy),length(xx));
vt = zeros(length(yy),length(xx));
pro = zeros(length(yy),length(xx));
eps = zeros(length(yy),length(xx));
%
for jj = js:je
    switch plane
        case  'x-y'
            j = jj; 
        case  {'x-z','y-z'}
            k = jj;
    end
    for ii = iss:iee
        switch plane
            case  {'x-y','x-z'}
                i = ii;
            case  'y-z'
                j = ii;
        end
        if nf(j,k,i) == - 1; continue; end
        nn = mn(j,k,i); 
        if nn ~=0
            % Make the array
            %if strcmp(output,'u')
            nnpi = mn(j+iy,k+iz,i+ix); 
            nnpj = mn(j+jy,k+jz,i+jx);
            if nnpi ~=0
                u1 = 0.5 * (u(nn,u1i) + u(nnpi,u1i)); 
            else
                u1 = 0.5 * u(nn,u1i); 
            end
            if nnpj ~=0
                u2 = 0.5 * (u(nn,u2i) + u(nnpj,u2i));
            else
                u2 = 0.5 * u(nn,u2i); 
            end
            ux(jj-2,ii-2)  = ux(jj-2,ii-2) + u1;
            uv(jj-2,ii-2)  = uv(jj-2,ii-2) + sqrt(u1^2 + u2^2);
        %elseif strcmp(output,'nu_t')
            vt(jj-2,ii-2)  = vt(jj-2,ii-2) + nu_t(nn);
            pro(jj-2,ii-2) = pro(jj-2,ii-2) + turb_pro(nn);      
            eps(jj-2,ii-2) = eps(jj-2,ii-2) + turb_eps(nn);   
        %end
        end
    end
end
ux(ux == 0) = NaN;
uv(uv == 0) = NaN;
vt(vt == 0) = NaN;
pro(pro == 0) = NaN;
eps(eps == 0) = NaN;

figure(1);
pcolor(xx,yy,uv)
shading interp
colormap(jet)
caxis([0 10])
title('Velocity')
colorbar
print('-r600','-dpng',['Presentation/' oname '_' plane 'Vel.png'])

figure(2);
pcolor(xx,yy,vt)
shading interp
colormap(jet)
caxis([0 2])
title('Turbulent Viscosity')
colorbar
print('-r600','-dpng',['Presentation/' oname '_' plane 'Vis.png'])

figure(3);
pcolor(xx,yy,pro)
shading interp
colormap(jet)
caxis([0 0.25])
title('Turbulent Production')
colorbar
print('-r600','-dpng',['Presentation/' oname '_' plane 'Pro.png'])

figure(4);
pcolor(xx,yy,eps)
shading interp
colormap(jet)
caxis([0 0.25])
title('Turbulent Dissipation')
colorbar
print('-r600','-dpng',['Presentation/' oname '_' plane 'Dis.png'])

% %% Make desired target into grid format for x-z plane
% xx = 0.5*(x(is(1):ie(1)-1) + x(is(1)+1:ie(1)));
% zz = 0.5*(z(is(3):ie(3)-1) + z(is(3)+1:ie(3)));
% xx = xx/1d6;
% dy = y(is(2)+1)-y(is(2));
% 
% y_want = 845140; % choose cross-section
% j  = find( y > y_want - dy/2 & y < y_want + dy/2);
% uz = zeros(length(zz),length(xx));
% uxz = zeros(length(zz),length(xx));
% vtz = zeros(length(zz),length(xx));
% proz = zeros(length(zz),length(xx));
% %
% for k = is(3):ie(3)-1
%     for i = is(1):ie(1)-1
%         if nf(j,k,i) == - 1; continue; end
%         nn = mn(j,k,i); 
%         if nn ~=0
%             % Make the array
%             %if strcmp(output,'u')
%                 nnpi = mn(j,k,i+1); nnpk = mn(j,k+1,i);
%                 if nnpi ~=0
%                     u1 = 0.5*(u(nn,1) + u(nnpi,1)); 
%                 else
%                     u1 = 0.5*u(nn,1); 
%                 end
%                 if nnpj ~=0
%                     u2 = 0.5*(u(nn,2) + u(nnpk,1));
%                 else
%                     u2 = 0.5*u(nn,2); 
%                 end
%                 uxz(k-is(2)+1,i-is(1)+1)= uxz(k-is(3)+1,i-is(1)+1) + ...
%                                           u(nn,1);
%                 uz(k-is(2)+1,i-is(1)+1) = uz(k-is(3)+1,i-is(1)+1) + ...
%                                           sqrt(u1^2 + u2^2);
%             %elseif strcmp(output,'nu_t')
%                 vtz(k-is(2)+1,i-is(1)+1) = vtz(k-is(3)+1,i-is(1)+1) + ...
%                                            nu_t(nn);
%                 proz(k-is(2)+1,i-is(1)+1) = proz(k-is(3)+1,i-is(1)+1) + ...
%                                            turb_pro(nn);
%             %end
%         end
%     end
% end
% uxz(uxz == 0) = NaN;
% uz(uz == 0) = NaN;
% vtz(vtz == 0) = NaN;
% proz(proz == 0) = NaN;
% 
% figure(4);
% pcolor(xx,zz,uz)
% shading interp
% colormap(jet)
% caxis([0 10])
% colorbar
% 
% figure(5);
% pcolor(xx,zz,vtz)
% shading interp
% colormap(jet)
% caxis([0 14])
% colorbar
% 
% figure(6);
% pcolor(xx,zz,proz)
% shading interp
% colormap(jet)
% caxis([0 1])
% colorbar
% 
%% Plot Um/U0 versus x/h
switch plane
    case  'x-z'
        z_want = -8; dz = 2.2;
        k  = find( z > z_want - dz/2 & z < z_want + dz/2);
        x0 = 408050; i = find(x == x0);
        Um = squeeze(ux(k,1:i));
        U0 = ux(k,i);
        h  = 20;
        xh = (x(i)-x(1:i))/h;

        figure;
        loglog(xh,Um/U0,'x')
        xlim([10 100])
        ylim([0 1.1])

        hold on
        
        if strcmp(oname,'DeepOn')
            I = find(xh >= 20 & xh <= 35);
            p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-')
            text(15,0.9 ,['[x/h]^{' num2str(p1(1)) '}'])

            I = find(xh >= 37 & xh <= 57);
            p2 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p2(2))*xh(I).^p2(1),'-')
            text(60,0.8 ,['[x/h]^{' num2str(p2(1)) '}'])

            I = find(xh >= 65 & xh <= 82);
            p3 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p3(2))*xh(I).^p3(1),'-')
            text(48,0.3 ,['[x/h]^{' num2str(p3(1)) '}'])
        elseif strcmp(oname,'Normal')
            I = find(xh >= 20 & xh <= 40);
            p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-')
            text(15,0.8 ,['[x/h]^{' num2str(p1(1)) '}'])
            
            I = find(xh >= 40 & xh <= 50);
            p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-')
            text(50,0.8 ,['[x/h]^{' num2str(p1(1)) '}'])
            
            I = find(xh >= 50 & xh <= 80);
            p2 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p2(2))*xh(I).^p2(1),'-')
            text(63,0.3 ,['[x/h]^{' num2str(p2(1)) '}'])

            I = find(xh >= 80 & xh <= 100);
            p3 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p3(2))*xh(I).^p3(1),'-')
            text(61,0.08 ,['[x/h]^{' num2str(p3(1)) '}'])
        elseif strcmp(oname,'Real_ke')
            I = find(xh >= 10 & xh <= 30);
            p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-')
            text(15,0.9 ,['[x/h]^{' num2str(p1(1)) '}'])
            
            I = find(xh >= 30 & xh <= 60);
            p2 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p2(2))*xh(I).^p2(1),'-')
            text(40,0.8 ,['[x/h]^{' num2str(p2(1)) '}'])

            I = find(xh >= 60 & xh <= 80);
            p3 = polyfit(log(xh(I)),log(Um(I)/U0),1);
            loglog(xh(I),exp(p3(2))*xh(I).^p3(1),'-')
            text(70,0.8 ,['[x/h]^{' num2str(p3(1)) '}'])            
            
        end
          
        xlabel('x/h')
        ylabel('U_{m0}/U_0')
        
        print('-r600','-dpng',['Presentation/' oname '_Umxh.png'])
end
% % [~,IDX] = max(abs(ux),[],1);
% % u_cross = zeros(length(IDX),1); vt_cross = zeros(length(IDX),1);
% % for i = 1:length(IDX)
% %     u_cross(i) = ux(IDX(i),i);
% %     vt_cross(i) = vt(IDX(i),i);
% % end
% % %u_cross = squeeze(uv(int64(end/2),:));
% % %vt_cross = squeeze(vt(int64(end/2),:));
% % 
% % figure(2);
% % [AX,~,~] = plotyy(xx,-u_cross,xx,vt_cross);
% % 
% % xm = 408150; i = find(x == xm); j = find(y >= 4.345d6 - 3.5d6 & y <= 4.345335d6 - 3.5d6);
% % um = squeeze(ux(j,i));
% % ym = y(j);
% % Qm = trapz(ym,um.^2);
% % alpha1 = 8;
% % CurlyX = sqrt(Qm*(3/64)/alpha1^3);
% % CurlyX = 0.05;
% % 
% % 
% % j = find(y >= 4.3449d6 - 3.5d6 & y <= 4.3455d6 - 3.5d6);
% % xend = 95;
% % P = polyfit((xx(xend+1)-xx(1:xend)).^-0.5,u_cross(1:xend),1);
% % 
% % axes(AX(1))
% % hold on
% % plot(xx(1:xend),P(1)*(xx(xend+1)-xx(1:xend)).^-0.5+P(2),'-b')
% % axes(AX(2))
% % hold on
% % plot(xx(1:xend),CurlyX*(1d6*xx(xend+1)-1d6*xx(1:xend)).^0.5,'-r')
% 
% figure(3);
% xe = 407300; %407600:-10:407200; 
% alpha = 0.116;
% x0    = xm + 1.5d3; % imaginary zero location
% for xee = xe
%     i = find(x == xee);
%     ue = squeeze(ux(j,i));
%     [umm , umI] = min(ue);
%     b = alpha*abs(xee - x0)/2;
%     %vte = squeeze(nu_t(mn(j,k,i)));
%     subplot(2,1,1)
%     %plot((y(j)-y(j(umI)))/b,ue/umm,'.')
%     plot(y(j)-y(j(umI)),-ue,'.')
%     hold on
%     %ueg = umm * exp(-((y(j)-y(j(umI)))/b).^2);
%     %ueg = umm * sech((y(j)-y(j(umI)))/b).^2;
%     ueg = 4*alpha1^2*CurlyX * sech(y(j)-y(j(umI))).^2;
%     plot(y(j)-y(j(umI)),ueg)
%     %plot((y(j)-y(j(umI)))/b,ueg/umm)
%     ue(isnan(ue)) = 0;
%     Qmg = trapz(y(j),ue.^2);
%     subplot(2,1,2)
%     plot(xee,Qmg,'o')
%     hold on
% end
% plot(xm,Qm,'o')
% %plot(y(j),vte)
% 
% % gaussian profile
% 
% %
% %plot(y(j),ueg)
% %Qmg = trapz(y(j),ueg.^2);
% 
% figure(1);
% hold on
% plot(x,alpha*abs(x - x0)/2)
% 
% figure(2);
% axes(AX(2))
% hold on
% b0 = y(end)-y(1);
% um(um == 0) = [];
% uc = 2.41*mean(abs(um))*sqrt(b0./(x0-x));
% plot(x,uc,'-b')