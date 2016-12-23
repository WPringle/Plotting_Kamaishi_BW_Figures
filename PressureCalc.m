% The script calculates the forces and moments on the submerged caisson due
% to the tsunami
clearvars;
clc;
close all;

g=9.807;   %gravity
rho=1.025; %specific gravity of salt water
Weight = (2.3*2.7 + 2.1*8.3)*g; %kN/m^2
Fc = 0.6;  % friction coefficient
%% 2DH part
direc = 'Satakev8.0_2DH_BW_FD/'; %'Satakev8.0_3D_Real_BW/';
fname = 'Satakev8.0_2DH_BW_FD'; %'Satakev8.0_3D_Real_BW_';
load([direc 'MAT_files/' fname 'xyzF.mat'])

is = find ( x >= 408130 & x <= 408175);
js = find ( y >= 845050 & y <= 845280);
load([direc 'MAT_files/' fname 'VelEta.mat'])

[time2D,tI] = sort(time); Eta = Eta(tI,:); time2D = time2D';

Pdiff = zeros(length(time2D),length(js));
Hmom = zeros(length(time2D),length(js));
Lift = zeros(length(time2D),length(js));
for t = 1:length(time2D)
    for j = 1:length(js)
        % find the i values where nf = -1 and get either side
        nn = squeeze(mn(js(j),is));
        I = find(Zk(nn) > -22);
        im = is(I(1)-1);
        ip = is(I(end)+1);
        nnm = mn(js(j),im);
        nnp = mn(js(j),ip);
        nnt = mn(js(j),is(I));
        Zkt = mean(Zk(nnt));
        % Get the pressure forces (kN/m) either side
        % integral(p.dz,zt,zb) 
        Pm = Eta(t,nnm)*(Zkt - Zk(nnm)) - 0.5*(Zkt^2 - Zk(nnm)^2);
        Pp = Eta(t,nnp)*(Zkt - Zk(nnp)) - 0.5*(Zkt^2 - Zk(nnp)^2);
        % Get the pressure diff
        Pdiff(t,j) =rho*g*(Pp - Pm);
        % Get Lift force
        Lift(t,j) = rho*g*(Zkt - 0.5*(Zk(nnm)+Zk(nnp)));
        % Get the moment force from caisson heel
        % integral(p*(z-zb).dz,zt,zb) 
        Mm = 0.5*Eta(t,nnm)*(Zkt^2 - Zk(nnm)^2) - ...
             Eta(t,nnm)*Zk(nnm)*(Zkt - Zk(nnm)) - ...
             1/3*(Zkt^3 - Zk(nnm)^3) + 0.5*Zk(nnm)*(Zkt^2 - Zk(nnm)^2);
        Mp = 0.5*Eta(t,nnp)*(Zkt^2 - Zk(nnp)^2) - ...
             Eta(t,nnp)*Zk(nnp)*(Zkt - Zk(nnp)) - ...
             1/3*(Zkt^3 - Zk(nnp)^3) + 0.5*Zk(nnp)*(Zkt^2 - Zk(nnp)^2);
        Hmom(t,j) = rho*g*(Mp - Mm);
    end
end
% pressure force
Phmin = min(Pdiff,[],2); %kN/m
Phmax = max(Pdiff,[],2); %kN/m
Phmean = mean(Pdiff,2);  %kN/m

% moment force
Mhmin = min(Hmom,[],2);
Mhmax = max(Hmom,[],2);
Mhmean = mean(Hmom,2);

% (weight - lift) * width :resistance force
Rwhmin = (Weight - max(Lift,[],2))*13; %kN/m
Rwhmax = (Weight - min(Lift,[],2))*13; %kN/m

% (weight - lift) * pivto arm : resistance moment
Rmhmin = (Weight - min(Lift,[],2))*13*6.5; %kN
Rmhmax = (Weight - max(Lift,[],2))*13*6.5; %kN

%% 3D part
direc = 'Satakev8.0_3D_Real_BW/';
load([direc 'MAT_files/Satake_BW_Real_ke.xyznF.mat'])

zs = find ( z >= -32 & z < -20.42);
is = find ( x >= 408130 & x <= 408175);
js = find ( y >= 845050 & y <= 845280);
load([direc 'MAT_files/Satake_BW_Real_ke.P.mat'])

[time,tI] = sort(time); P3 = P3(tI,:);

Pdiff = zeros(length(time),length(js));
Pmom = zeros(length(time),length(js));
Lift = zeros(length(time),length(js));
for t = 1:length(time)
    for j = 1:length(js)
        % find the i values where nf = -1 and get either side
        nfi = squeeze(nf(js(j),zs(1),is));
        I = find(nfi == -1);
        im = is(I(1)-1);
        ip = is(I(end)+1);
        nnm = squeeze(mn(js(j),zs,im));
        nnp = squeeze(mn(js(j),zs,ip));
        % dz values
        dz = diff(z([zs, zs(end)+1]));
        % Get the pressure forces (kN/m) either side
        Pm = sum(P3(t,nnm).*dz);
        Pp = sum(P3(t,nnp).*dz);
        % Get the pressure diff
        Pdiff(t,j) = (Pp - Pm)*rho;
        % Get the pressure forces (kN/m) on top
        nnt = squeeze(mn(js(j),zs(end)+1,is(I)));
        Pt  = sum(P3(t,nnt).*diff(x([is(I), is(I(end))+1])))./...
              sum(diff(x([is(I), is(I(end))+1])));
        % Get the pressures on bottom, estimate from average of either side
        Pb  = 0.5*(P3(t,nnm(1))+P3(t,nnp(1)));
        % Get the lift force
        Lift(t,j) = (Pb - Pt)*rho;
        % Get the moment force from caisson heel
        % integral(p*(z-zb).dz,zt,zb) 
        % mid-z values
        zhalf = 0.5*(z(zs(1):zs(end)) + z(zs(1)+1:zs(end)+1));
        Mm = sum(P3(t,nnm).*dz.*(zhalf-z(zs(1))));
        Mp = sum(P3(t,nnp).*dz.*(zhalf-z(zs(1))));
        Pmom(t,j) = rho*(Mp - Mm);
    end
end
% Get the drag force
Pdmin = min(Pdiff,[],2); %kN/m
Pdmax = max(Pdiff,[],2); %kN/m
Pdmean = mean(Pdiff,2); %kN/m

% Get the overturning moment
Memin = min(Pmom,[],2);
Memax = max(Pmom,[],2);
Memean = mean(Pmom,2);

% Total weight force is weight - lift (need to multiply by width to get kN/m)
Wemin = (Weight - max(Lift,[],2))*13; %kN/m
Wemax = (Weight - min(Lift,[],2))*13; %kN/m
Wemean = (Weight - mean(Lift,2))*13; %kN/m

% (weight - lift) * pivot arm : resistance moment
Rmemin = (Weight - min(Lift,[],2))*13*6.5; %kN
Rmemax = (Weight - max(Lift,[],2))*13*6.5; %kN
Rmemean = (Weight - mean(Lift,2))*13*6.5; %kN

%% Plot the forces
% Plot the drag force
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.09 0.010], [0.13 0.022]);
figure;
subplot(2,1,1);
hold on
[~,~] = jbfill(time2D/60,abs(Phmax),abs(Phmin),'b','none',0,0.5);
[~,~] = jbfill(time2D/60,Rwhmax*Fc,Rwhmin*Fc,'m','m',0,0.5);
[~,~] = jbfill(time/60,abs(Pdmax),abs(Pdmin),'k','none',0,0.5);
[~,~] = jbfill(time/60,Wemax*Fc,Wemin*Fc,'r','none',0,0.5);
plot(time2D/60, abs(Phmean),'b');
plot(time/60, abs(Pdmean),'k');
plot(time/60, Fc*Wemean,'r');
xlim([0 120])
ylim([0 2200])
text(5,2100,'(a)')
set(gca,'xticklabel','')
ylabel('|F_D| [kN/m]')
legend('2DH (Hydrostatic) Forcing','Hydrostatic Resistance',...
       '2CLOWNS-3D Forcing','2CLOWNS-3D Resistance',...
       'location','best')
legend boxoff
set(gca,'fontsize',7)
%print('-r600','-dpng','../Paper/Real_ke_Drag.png')
   
% Plot the overturning moment
subplot(2,1,2);
hold on
[~,~] = jbfill(time2D/60,abs(Mhmax)/1000,abs(Mhmin)/1000,'b','none',0,0.5);
[~,~] = jbfill(time2D/60,Rmhmax/1000,Rmhmin/1000,'m','m',0,0.5);
[~,~] = jbfill(time/60,abs(Memax)/1000,abs(Memin)/1000,'k','none',0,0.5);
[~,~] = jbfill(time/60,Rmemax/1000,Rmemin/1000,'r','none',0,0.5);
plot(time2D/60, abs(Mhmean)/1000,'b');
plot(time/60, abs(Memean)/1000,'k');
plot(time/60, Rmemean/1000,'r');
text(5,12,'(b)')
xlim([0 120])
xlabel('time since earthquake [min]')
ylabel('|M_D| [MN.m/m]')
% legend('2DH (Hydrostatic) Overturning Moment',...
%        'Hydrostatic Resisting Moment',...
%        '2CLOWNS-3D Overturning Moment',...
%        '2CLOWNS-3D Resisting Moment','location','East')
legend boxoff
set(gca,'fontsize',7)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 9 10]) %...
    %'PaperPositionMode','manual');

print('-r600','-dpdf','../Paper/Submerged_caisson_force.pdf')

% temax = Memax./Wemax; thmax = Mhmax./Wemax; 
% temin = Memin./Wemin; thmin = Mhmin./Wemin;
% Pemax = zeros(length(time),1);Pemin = zeros(length(time),1);
% Phmax = zeros(length(time),1);Phmin = zeros(length(time),1);
% for t = 1:length(time)
%     if temax(t) <= 13/3
%         Pemax(t) = (2*Wemax(t))/(3*temax(t));
%     else
%         Pemax(t) = (2*Wemax(t)/13)*(2-3*temax(t)/13);
%     end
%     if temin(t) <= 13/3
%         Pemin(t) = (2*Wemin(t))/(3*temin(t));
%     else
%         Pemin(t) = (2*Wemin(t)/13)*(2-3*temin(t)/13);
%     end
%     if thmax(t) <= 13/3
%         Phmax(t) = (2*Wemax(t))/(3*thmax(t));
%     else
%         Phmax(t) = (2*Wemax(t)/13)*(2-3*thmax(t)/13);
%     end
%     if thmin(t) <= 13/3
%         Phmin(t) = (2*Wemin(t))/(3*thmin(t));
%     else
%         Phmin(t) = (2*Wemin(t)/13)*(2-3*thmin(t)/13);
%     end
% end\
% subplot1(3, 1,'Gap', [.05 .05],'Max',[1.02 1.03],'Min',[0.1 0.09],...
%          'XTickL', 'Margin', 'YTickL', 'Margin');
% subplot1(1);
% hold on
% [ph,msg] = jbfill(time/60,Fc*Wemax./max(abs(Pdiff),[],2),...
%                   Fc*Wemin./min(abs(Pdiff),[],2),...
%              'b','b',0,0.5);
%figure;

%[~,~] = jbfill(time/60,-Whmax*Fc,-Whmin*Fc,'m','m',0,0.5);
%[~,~] = jbfill(time/60,-Wemax*Fc,-Wemin*Fc,'r','r',0,0.5);




%set(gca,'XTickLabel','');
%ylabel('$F_D$ [kN/m]','interpreter','latex')
%ylabel('$\Delta \eta$ [m]','interpreter','latex')
% [LX,LY] = gca_latex(8,0,1);
% axes(ha(2));
% subplot1(2);
% hold on
% [ph,msg] = jbfill(time/60,max(Pmom,[],2)/1e3,min(Pmom,[],2)/1e3,...
%            'b','b',1,0.5);
% hold on
% [ph,msg] = jbfill(time/60,max(Hmom,[],2)/1e3,min(Hmom,[],2)/1e3,...
%            'k','k',1,0.5);
% hold on
% [ph,msg] = jbfill(time/60,(Weight-max(Lift,[],2))*6.5/1e3,...
%            (Weight-min(Lift,[],2))*6.5/1e3,...
%            'r','r',1,0.5);
% % xlim([time(1)/3600 time(length(time))/3600])
% xlim([17 time(length(time))/60])
% % set(gca,'XTickLabel','');
% ylabel('$M_D$ [MNm/m]','interpreter','latex')
% % [LX,LY] = gca_latex(8,0,1);
% % axes(ha(3));
% subplot1(3);
% [ph,msg] = jbfill(time/60,Pemax,Pemin,...
%            'b','b',1,0.5);
% hold on
% [ph,msg] = jbfill(time/60,Phmax,Phmin,...
%            'k','k',1,0.5);
% hold on
% plot([time(1)/60 time(length(time))/60],[600 600],'r')
% plot([time(1)/60 time(length(time))/60],[800 800],'--r')
% %xlim([time(1)/3600 time(length(time))/3600])
% xlim([17 time(length(time))/60])
% ylim([0 1000])
% ylabel('$P_e$ [kPa]','interpreter','latex')
% xlabel('time since earthquake [min]')
% % [LX,LY] = gca_latex(8,1,1);
% for i = 1:length(time)
%     std1(i) = std(Hdiff(i,:));
% end
% k = 0;
% for i = 1:length(time2D)
%     if time2D(i) < time(1); continue; end
%     if time2D(i) > time(length(time)); break; end
%     k = k +1;
%     timenew(k) = time2D(i);
%     std2(k) = std(Hdiff2D(i,:));
% end