% Script to analyse the jet through the submerged section of the breakwater
% as computed by 2CLOWNS-3D simulation with time. 
% By: William Pringle Oct 14 2016

clearvars; clc; close all;
%% Set and load the data
% direc = 'Satakev8.0_ke/';
% fname = 'Tohoku_SatakeT_HYB_ke_k1mm_wall.';
% oname = 'ke_wall';

direc = 'Satakev8.0_3D_Real_BW/';
fname = 'Satake_BW_Real_ke.';
oname = 'Real_ke';
% direc = 'Satakev8.0_ke_HR_DeepOn/';
% fname = 'Tohoku_SatakeT_HYB_ke_k1mm_DeepOn_HR.';
% oname = 'DeepOn';

fPath = cd;
List = dir(fullfile(fPath, [direc fname '*']));
List([List.isdir]) = [];
L = length(List)/4;
% order list by datenum
Listc = struct2cell(List);
Ordered = Listc'; %sortrows(Listc',5);

load([direc 'MAT_files/' fname 'xyzn.mat']);  

%% Choose wanted points 
if strcmp(direc,'Satakev8.0_3D_Real_BW/')
    L = 210;
end
dz = 2.2;  % vertical cell size
dy = 10;   % horizontal cell size

z_want  = -20:dz:7; % vertical contour
y_want  = 845100:dy:845160; % horizontal contour  845140
x0      = 408050;  % 0 point  of jet
xh_want = 0:25:75; % points to measure velocity

h  = 25;   % height of nozzle

k  = knnsearch(z',z_want');
j  = knnsearch(y',y_want'); 
i  = find(x == x0);
xh = (x(i)-x)/h;
I = knnsearch(xh',xh_want');
I(I == 1) = 2;
% Initialise the ux_want and time
ux_want = zeros(L,length(I));
ux_all = zeros(L,ie(1)); c_mu_all = zeros(L,ie(1)); 
uz_all = zeros(L,ie(1),k(end)); kz_all = zeros(L,ie(1),k(end)); 
uy_all = zeros(L,ie(1),ie(2)-1); ky_all = zeros(L,ie(1),ie(2)-1);
vz_all = zeros(L,ie(1),k(end)); vy_all = zeros(L,ie(1),ie(2)-1);
tz_all = zeros(L,ie(1),k(end)); ty_all = zeros(L,ie(1),ie(2)-1);
time  = zeros(L,1);
%% Loop over the time
for t = 1:L
    
    Ordern = Ordered(t*4-3:t*4,:);
    Ordern = sortrows(Ordern,-3);
    % Get time
    File = Ordern{1}; FL = length(File);
    tstrng = File(FL-5:FL);
    if strcmp(tstrng(1),'t')
        tstrng(1) = [];
    end
    timen = str2double(tstrng);
    time(t) = timen;
    
    load([direc 'MAT_files/' File '.mat'])
    
    % Get the velocities at wanted points
    % Make the Ux vector for all x points
    ux = zeros(ie(1),1); c_mu = zeros(ie(1),1);
    uz05 = zeros(ie(1),k(end));  kz05 = zeros(ie(1),k(end)); 
    vz05 = zeros(ie(1),k(end));  tz05 = zeros(ie(1),k(end));
    uy05 = zeros(ie(1),ie(2)-1); ky05 = zeros(ie(1),ie(2)-1); 
    vy05 = zeros(ie(1),ie(2)-1); 
    for ii = is(1):ie(1)
        nn = mn(j,k,ii);
        ux(ii) = max(-u(nn,1));
        c_mu_n = nu_t(nn).*eps(nn)./rk(nn).^2; 
        c_mu(ii) = mean(mean(c_mu_n(~isnan(c_mu_n))));
        for kk = is(3):k(end)
            nnk = mn(j,kk,ii);
            if nnk ~= 0
                uz05(ii,kk) = max(-u(nnk,1));
                if nnk <= length(rk); 
                    kz05(ii,kk) = max(rk(nnk));
                    vz05(ii,kk) = max(nu_t(nnk));
                end
            end
        end
        for jj = is(2):ie(2)-1
            nnj = mn(jj,k,ii);
            if nnj ~= 0
                uy05(ii,jj) = max(-u(nnj,1));
                if nnj <= length(rk); 
                    ky05(ii,jj) = max(rk(nnj));
                    vy05(ii,jj) = max(nu_t(nnj));           
                end
            end
        end
    end
    ux_want(t,:) = ux(I);
    ux_all(t,:)  = ux; c_mu_all(t,:) = c_mu;
    uz_all(t,:,:) = uz05; kz_all(t,:,:) = kz05; vz_all(t,:,:) = vz05; 
    uy_all(t,:,:) = uy05; ky_all(t,:,:) = ky05; vy_all(t,:,:) = vy05; 
end

[time,tI] = sort(time); ux_want = ux_want(tI,:);
%% Plot points versus time 
figure;
plot(time,ux_want,'-o')

legend(['x/h = ' num2str(xh_want(1))],['x/h = ' num2str(xh_want(2))],...
       ['x/h = ' num2str(xh_want(3))],['x/h = ' num2str(xh_want(4))],...
       'location','best') %['x/h = ' num2str(xh_want(5))]

xlim([1300 2500])
ylim([0 14])
xlabel('time [s]')
ylabel('U_{m} [m/s]')
title('Transient nature of maximum velocities onshore of breakwater')

%print('-r600','-dpng',['../Presentation/' oname '_Umt.png'])

%% Plot min c_mu versus location downstream
figure;
c_mu_min = min(c_mu_all); %?;prctile(c_mu_all,10,1); %;
plot(xh(1:i),c_mu_min(1:i),'-')
xlabel('\itx/h')
xlim([0 78])
ylim([0.03 0.09])
ylabel('\itC_{\mu}')
%title('10th percentile values of C_\mu with distance from breakwater')
%print('-r600','-dpng',['../papers/' oname '_c_mu.png'])
set(gca,'fontsize',7,'Yticklabel',0.03:0.01:0.09)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 9 6],'PaperPositionMode','manual');
print('-r600','-depsc','../Paper/C_mu_min.eps'); %'-fillpage',
%% Time average ux_all for decay analysis
U90 = prctile(ux_all,90,1);
Um = zeros(1,i);
for ii = 1:i
    Um(ii)  = mean(ux_all(ux_all(:,ii) > U90(ii),ii));
end
U0 = Um(i);

figure; 
loglog(xh(1:i),Um/U0,'.')
% 
hold on
%loglog(101,0,'.')

%loglog(xh(1:i),U85(1:i)/U0)
if strcmp(direc,'Satakev8.0_Norm_Real_ke/') || strcmp(direc,'Satakev8.0_3D_Real_BW/')
    %I = find(xh >= 1 & xh <= 15);
    I = find(xh >= 1 & xh <= 20);
    p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
    loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
    text(5,0.83 ,['[x/h]^{' num2str(p1(1)) '}'])

    %I = find(xh >= 15 & xh <= 60);
    I = find(xh >= 18 & xh <= 50);
    p2 = polyfit(log(xh(I)),log(Um(I)/U0),1);
    loglog(xh(I),exp(p2(2))*xh(I).^p2(1),'-','LineWidth',1.25)
    text(26,0.73 ,['[x/h]^{' num2str(p2(1)) '}'])

    %I = find(xh >= 60 & xh <= 85);
    I = find(xh >= 48 & xh <= 66);
    p3 = polyfit(log(xh(I)),log(Um(I)/U0),1);
    loglog(xh(I),exp(p3(2))*xh(I).^p3(1),'-','LineWidth',1.25)
    text(39,0.54 ,['[x/h]^{' num2str(p3(1)) '}'])  
elseif strcmp(direc,'Satakev8.0_ke/')
    I = find(xh >= 1 & xh <= 17);
    p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
    loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
    text(6,0.80 ,['[x/h]^{' num2str(p1(1)) '}'])

    I = find(xh >= 15 & xh <= 35);
    p2 = polyfit(log(xh(I)),log(Um(I)/U0),1);
    loglog(xh(I),exp(p2(2))*xh(I).^p2(1),'-','LineWidth',1.25)
    text(18,0.75 ,['[x/h]^{' num2str(p2(1)) '}'])

    I = find(xh >= 32 & xh <= 40);
    p3 = polyfit(log(xh(I)),log(Um(I)/U0),1);
    loglog(xh(I),exp(p3(2))*xh(I).^p3(1),'-','LineWidth',1.25)
    text(40,0.75 ,['[x/h]^{' num2str(p3(1)) '}'])    
end

grid on 
xlabel('x/h')
ylabel('U_{m}/U_0')
xlim([4 80])
ylim([0.1 1.0])
ax = gca;
ax.XTick = [4 6 8 10 20 40 60];
ax.YTick = [0.1 0.2 0.4 0.6 0.8 1.0];
title('Decay of maximum centreline velocity from breakwater')

%print('-r600','-dpng',['../Presentation/' oname '_Umxh.png'])
    
%% Look at z and y growth rates 
Uz90 = squeeze(prctile(uz_all,90,1)); Kz90 = squeeze(prctile(kz_all,90,1));
Vz90 = squeeze(prctile(vz_all,90,1)); Tz90 = squeeze(prctile(tz_all,90,1));
Uz05 = zeros(i,k(end)); Kz05 = zeros(i,k(end)); 
Vz05 = zeros(i,k(end)); Tz05 = zeros(i,k(end));
zhalf = zeros(i,1);
Uy90 = squeeze(prctile(uy_all,90,1)); Ky90 = squeeze(prctile(ky_all,90,1));
Vy90 = squeeze(prctile(vy_all,90,1)); Ty90 = squeeze(prctile(ty_all,90,1));
Uy05 = zeros(i,ie(2)-1); Ky05 = zeros(i,ie(2)-1); 
Vy05 = zeros(i,ie(2)-1); Ty05 = zeros(i,ie(2)-1);
yhalf1 = zeros(i,1);
yhalf2 = zeros(i,1);
for ii = i:-1:1
    % z growth rate
    for kk = 1:k(end)
        Uz05(ii,kk)  = mean(uz_all(uz_all(:,ii,kk) > Uz90(ii,kk),ii,kk));
        Kz05(ii,kk)  = mean(kz_all(kz_all(:,ii,kk) > Kz90(ii,kk),ii,kk));
        Vz05(ii,kk)  = mean(vz_all(vz_all(:,ii,kk) > Vz90(ii,kk),ii,kk));
        Tz05(ii,kk)  = mean(tz_all(vz_all(:,ii,kk) > Tz90(ii,kk),ii,kk));
    end
    I = find(Uz05(ii,:) <= 0.5*Um(ii)); 
    if ~isempty(I)
        c1 = 0.5*Um(ii) - Uz05(ii,I(end));
        c2 = Uz05(ii,I(end)+1) - 0.5*Um(ii);
        c3 = Uz05(ii,I(end)+1) - Uz05(ii,I(end));
        zhalf(ii) = z(I(end)) * c2/c3 + z(I(end)+1) * c1/c3;
    end
    % y growth rate
    for jj = is(2):ie(2)-1
        Uy05(ii,jj)  = mean(uy_all(uy_all(:,ii,jj) > Uy90(ii,jj),ii,jj));
        Ky05(ii,jj)  = mean(ky_all(ky_all(:,ii,jj) > Ky90(ii,jj),ii,jj));
        Vy05(ii,jj)  = mean(vy_all(vy_all(:,ii,jj) > Vy90(ii,jj),ii,jj));
        Ty05(ii,jj)  = mean(ty_all(ty_all(:,ii,jj) > Ty90(ii,jj),ii,jj));
    end
    I = find(Uy05(ii,:) <= 0.5*Um(ii)); 
    if ~isempty(I) && length(I) > 1
        Idiff = diff(I); [~,ix] = max(Idiff);
        c1 = 0.5*Um(ii) - Uy05(ii,I(ix));
        c2 = Uy05(ii,I(ix)+1) - 0.5*Um(ii);
        c3 = Uy05(ii,I(ix)+1) - Uy05(ii,I(ix));
        yhalf1(ii) = y(I(ix)) * c2/c3 + y(I(ix)+1) * c1/c3;
        c1 = 0.5*Um(ii) - Uy05(ii,I(ix+1));
        c2 = Uy05(ii,I(ix+1)-1) - 0.5*Um(ii);
        c3 = Uy05(ii,I(ix+1)-1) - Uy05(ii,I(ix+1));
        yhalf2(ii) = y(I(ix+1)) * c2/c3 + y(I(ix+1)-1) * c1/c3;
    end   
end
y0 = 845130;
yhalf1 = y0 - yhalf1;
yhalf2 = yhalf2 - y0;
zhalf  = z(k(end)) - zhalf;

figure;
loglog(xh(1:i),zhalf/h,'.')
hold on
loglog(xh(1:i),yhalf1/h,'.')
loglog(xh(1:i),yhalf2/h,'.')
if strcmp(direc,'Satakev8.0_Norm_Real_ke/') || strcmp(direc,'Satakev8.0_3D_Real_BW/')
    I = find(xh >= 4 & xh <= 6.5); I(isnan(zhalf(I))) = [];
    p1 = polyfit(log(xh(I)),log(zhalf(I)'/h),1);
    loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
    text(2.7,2.6 ,['[x/h]^{' num2str(p1(1)) '}'])
    
%     I = find(xh >= 20 & xh <= 70); I(isnan(yhalf2(I))) = [];
%     p1 = polyfit(log(xh(I)),log(yhalf2(I)'/h),1);
%     loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
%     text(30,8 ,['[x/h]^{' num2str(p1(1)) '}'])
elseif strcmp(direc,'Satakev8.0_ke/')
%     I = find(xh >= 1 & xh <= 4);
%     p1 = polyfit(log(xh(I)),log(zhalf(I)'/h),1);
%     loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
%     text(3,1.5 ,['[x/h]^{' num2str(p1(1)) '}'])
end
legend('z_{1/2}/h','yS_{1/2}/h','yN_{1/2}/h','Location','best')

grid on
xlim([1 80])
ylim([1 10])
xlabel('x/h')
ylabel('z_{1/2}/h,y_{1/2}/h')
title('Half width growth rates')
ax = gca;
ax.XTick = [1 2 3 4 6 8 10 20 40 60];
ax.YTick = [1 2 4 6 8 10];
%print('-r600','-dpng',['../Presentation/' oname '_Growthxh.png'])

%% x-y plane velocities
figure;
pcolor(x(is(1):i),y(is(2):ie(2)-1)+3.5d6,Uy05(is(1):i,is(2):ie(2)-1)')
shading interp
colormap jet
caxis([0 10])
colorbar
xlabel('Easting [m]')
ylabel('Northing [m]')
title('Mean of 90th percentile of velocities in x-y plane')
%print('-r600','-dpng',['../Presentation/' oname '_x-yVel.png'])

%% x-z plane velocities 
%Uz05(Uz05 == 0) = NaN;
figure;
pcolor(x(is(1):i),z(is(3):k(end)),Uz05(is(1):i,is(3):k(end))')
shading interp
colormap jet
caxis([0 10])
colorbar
xlabel('Easting [m]')
ylabel('z [m]')
title('Mean of 90th percentile of velocities in x-z plane')
%print('-r600','-dpng',['../Presentation/' oname '_x-zVel.png'])

%% x-y plane TKE
figure;
pcolor(x(is(1):i),y(is(2):ie(2)-1)+3.5d6,Ky05(is(1):i,is(2):ie(2)-1)')
shading interp
colormap jet
caxis([0 5])
colorbar
xlabel('Easting [m]')
ylabel('Northing [m]')
title('Mean of 90th percentile of TKE in x-y plane')
%print('-r600','-dpng',['../Presentation/' oname '_x-yTKE.png'])

%% x-z plane TKE
figure;
pcolor(x(is(1):i),z(is(3):k(end)),Kz05(is(1):i,is(3):k(end))')
shading interp
colormap jet
caxis([0 8])
colorbar
xlabel('Easting [m]')
ylabel('z [m]')
title('Mean of 90th percentile of TKE in x-z plane')
%print('-r600','-dpng',['../Presentation/' oname '_x-zTKE.png'])

%% x-y plane Turb Vis
figure;
pcolor(x(is(1):i),y(is(2):ie(2)-1)+3.5d6,Vy05(is(1):i,is(2):ie(2)-1)')
shading interp
colormap jet
caxis([0 15])
colorbar
xlabel('Easting [m]')
ylabel('Northing [m]')
title('Mean of 90th percentile of \nu_t in x-y plane')
%print('-r600','-dpng',['../Presentation/' oname '_x-yVis.png'])

%% x-z plane Turb Vis
figure;
pcolor(x(is(1):i),z(is(3):k(end)),Vz05(is(1):i,is(3):k(end))')
shading interp
colormap jet
caxis([0 15])
colorbar
xlabel('Easting [m]')
ylabel('z [m]')
title('Mean of 90th percentile of \nu_t in x-z plane')
%print('-r600','-dpng',['../Presentation/' oname '_x-zVis.png'])

%% x-y plane turbulent intensity
figure;
pcolor(x(is(1):i),y(is(2):ie(2)-1)+3.5d6,Ty05(is(1):i,is(2):ie(2)-1)')
shading interp
colormap jet
caxis([0 1.0])
colorbar
xlabel('Easting [m]')
ylabel('Northing [m]')
title('Mean of 90th percentile of TI in x-y plane')
print('-r600','-dpng',['../Presentation/' oname '_x-yTI.png'])

% %% x-z plane Turb Vis
% figure;
% pcolor(x(is(1):i),z(is(3):k(end)),Tz05(is(1):i,is(3):k(end))')
% shading interp
% colormap jet
% caxis([0 1.0])
% colorbar
% xlabel('Easting [m]')
% ylabel('z [m]')
% title('Mean of 90th percentile of \nu_t in x-z plane')
% print('-r600','-dpng',['../Presentation/' oname '_x-zTI.png'])

%% Look at axial velocity profiles 
Axial = zeros(size(Uz05)); Axialz = zeros(size(Uz05)); 
Along = zeros(size(Uy05)); 
Alongy1 = zeros(size(Uy05)); Alongy2 = zeros(size(Uy05)); 
for ii = 1:length(Um)
    Axial(ii,:) = Uz05(ii,:)/Um(ii);
    Axialz(ii,:)  = z(1:k(end))/zhalf(ii);
    %
    Along(ii,:) = Uy05(ii,:)/Um(ii);
    Alongy1(ii,:) = (y0-y(1:ie(2)-1))/yhalf1(ii); 
    Alongy2(ii,:) = (y(1:ie(2)-1)-y0)/yhalf2(ii); 
end

if strcmp(direc,'Satakev8.0_Norm_Real_ke/') || strcmp(direc,'Satakev8.0_3D_Real_BW/')
%     figure;
%     I = find(xh >= 18 & xh <= 50);
% %     plot(Axial(I,:)',Axialz(I,:)','.')
%     plot(Along(I,:)',Alongy1(I,:)','.')  
%     hold on
%     plot(Along(I,:)',Alongy2(I,:)','.') 
%     xlabel('U/U_{m}')
%     ylabel('y/y_{1/2}')
%     ylim([-0.4 1.6])
%     ax = gca;
%     ax.YTick = -0.4:0.4:1.6;
%     title('Axial velocity profile in characteristic decay region')
%     print('-r600','-dpng',['../Presentation/' oname '_charVel.png'])  
    
    figure;
    %I = knnsearch(xh',[65 70 75 80]');
    I = find(xh >= 18 & xh <= 60);
    plot(Along(I,:)',Alongy1(I,:)','.')  
    hold on
    plot(Along(I,:)',Alongy2(I,:)','.') 
    ylim([-0.4 1.6])
    ax = gca;
    ax.YTick = -0.4:0.4:1.6;

    grid on
%     plot(Axial(I,:)',Axialz(I,:)','.')  
%     hold on
%     zz = 0:0.01:2; B = 1; A =5; n = 0.7;
%     uu = A*(zz+0.3).^(1/n).*(1-erf(B*(zz+0.3)));
%     plot(uu,zz)
    
elseif strcmp(direc,'Satakev8.0_ke/')
%     figure;
%     I = find(xh >= 15 & xh <= 35);
%     plot(Axial(I,:)',Axialz(I,:)','.')
%     xlabel('U/U_{m}')
%     ylabel('z/z_{1/2}')
%     title('Axial velocity profile in characteristic decay region')
%     print('-r600','-dpng',['Presentation/' oname '_charVel.png'])
% 
%     figure;
%     I = find(xh >= 45 & xh <= 60);
%     plot(Axial(I,:)',Axialz(I,:)','.')  
%     xlabel('U/U_{m}')
%     ylabel('z/z_{1/2}')
%     title('Axial velocity profile in radial decay region')
%     print('-r600','-dpng',['Presentation/' oname '_radialVel.png'])

    figure;
    %I = knnsearch(xh',[65 70 75 80]');
    I = find(xh >= 15 & xh <= 40);
    plot(Along(I,:)',Alongy1(I,:)','.')  
    hold on
    plot(Along(I,:)',Alongy2(I,:)','.') 
    ylim([-0.4 1.6])
    ax = gca;
    ax.YTick = -0.4:0.4:1.6;

    grid on

end
 xlabel('U/U_{m}')
 ylabel('y/y_{1/2}')
 title('Transverse velocity profiles')
% print('-r600','-dpng',['../Presentation/' oname '_transverseVel.png'])
