%%%%%%%%%%%%% Get the ordered survey data and compare with calculation data
clearvars;
close all;

%% Get the data
% Survey measurements
%load('../Inputs/Survey_data_order.mat'); % Contains x, y and z data
%N = length(height_order);      % number of survey locations

M = csvread('../../Inputs/ttjt_survey_Ryoishi_Kamaishi_Toni_tidecorrected.csv');

N = length(M);

[B_x, B_y] = latlon_to_UTM(M(:,2), M(:,1), 0, 143);
B_x = B_x + 500000;
%B_y = B_y - 3.5d6;

% Calculation directory
dir_num{1} = '../Satakev8.0_2DH_BW_FD\DAT\'; 
dir_num{2} = '../Satakev8.0_ke\DAT\'; 
%dir_num{2} = 'Satakev8.0_3D_Real_BW\DAT\'; 
legend_name = {'2DH NSWE','2CLOWNS-3D'}; 
%
%
style = {'-k','--k',':k'};
color = {'b','r'};
line = {'-','--',':'};
marker = {'ok','+k','*k'};
L = length(dir_num);
hlines = 2; % header skipping

%% Open all the files and get the max data from the calculation
eta_max = zeros(L,N); depth_max = zeros(L,N);
for p = 1:L
    for m = 1:N
        m_str = num2str(m);
        num = dlmread([dir_num{p} 'LOC' m_str '.dat'],'',hlines,0);
        eta_max(p,m) = GetMaxValue( num, 0 , 1);
        depth_max(p,m) = GetMaxValue( num, 0 , 2);
        zk(p,m) = num(1,2);
        
        clear num
    end
end

depth_lim = 0.0;
%Delete starting and ending areas (of extremely large inundations)
% I = [1 2 3 4 8 82 83]';
% eta_max(:,I) = [];
% height_order(I) = [];
% Delete locations where survey is unreliable
%I = find(isnan(height_order));

% Find point in Kamaishi
I = find(B_y > 4.3425d6 & B_y < 4.348d6); 

eta_kam = eta_max(:,I);
depth_kam = depth_max(:,I);
zk_kam = zk(:,I);
m_kam = M(I,3);
B_x = B_x(I); B_y = B_y(I);


[~,K] = find(isnan(eta_kam) | depth_kam < depth_lim ...
             | [m_kam'; m_kam'] - zk_kam < depth_lim);
eta_kam(:,K) = [];
depth_kam(:,K) = [];
m_kam(K) = [];
zk_kam(K) = [];
B_x(K) = []; B_y(K) = []; 

 % Delete inside Kamaishi
eta_max(:,I) = [];
depth_max(:,I) = [];
zk(:,I) = [];
M(I,:) = [];

[~,K] = find(isnan(eta_max) | depth_max < depth_lim ...
             | [M(:,3)';M(:,3)'] - zk < depth_lim);
eta_max(:,K) = [];
depth_max(:,K) = [];
M(K,:) = [];
N = length(M);

%Plot the results with least squares line
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.07 0.06], [0.06 0.018]);
figure;
subplot(1,L,1);
for p = 1:L
    %plot(height_order,eta_max(p,:),'bd','MarkerFaceColor','b')
    if p == 1
        plot(M(:,3),eta_max(p,:),'bd',...
             'MarkerFaceColor',color{p},'MarkerSize',3)
        hold on
    else
        plot(M(:,3),eta_max(p,:),'rs',...
            'MarkerFaceColor',color{p},'MarkerSize',3)
    end
    %lsline

    ylabel('Simulation [m]')
    if p == L
        xlabel('Observation [m]')
        plot(0:0.1:35,0:0.1:35,'k--')
    else
        set(gca,'xticklabel','')
    end
    xlim([0 35])
    LogKi = log(M(:,3)./eta_max(p,:)');
    LogK  = mean(LogKi);
    K = exp(LogK);
    kappa = exp(sqrt(mean(LogKi.^2 - LogK^2)));
    text(2,33-(p-1)*4,['\it{K} = ' num2str(K,'%.2f') ...
        ', \kappa = ' num2str(kappa,'%.2f')],'color',color{p})
    %text(2,33-(p-1)*6,['Geometric Mean: \it{K} = ' num2str(K,'%.2f')],...
    %    'color',color{p})
    %text(2,30-(p-1)*6,['Geometric Standard deviation: \kappa = ' ...
    %   num2str(kappa,'%.2f')],'color',color{p})
    title('Toni, Ryoishi, Otsuchi Runup/Inundation')
    legend('2DH','2CLOWNS-3D','location','SouthEast')
    %legend boxoff
    set(gca,'fontsize',7)
    %xlim([min(height_order) 12])
    %ylim([min(height_order) 12])
end
% 
subplot(1,L,2);
for p = 1:L
    %plot(height_order,eta_max(p,:),'bd','MarkerFaceColor','b')
    if p == 1
        plot(m_kam,eta_kam(p,:),'bd',...
            'MarkerFaceColor',color{p},'MarkerSize',3)
        hold on
    else
        plot(m_kam,eta_kam(p,:),'rs',...
            'MarkerFaceColor',color{p},'MarkerSize',3)
    end
    %lsline
    %ylabel('Simulation [m]')
    set(gca,'yticklabel','')
    xlim([0 35])
    if p == L
        xlabel('Observation [m]')
        plot(0:0.1:35,0:0.1:35,'k--')
    else
        set(gca,'xticklabel','')
    end
    LogKi = log(m_kam./eta_kam(p,:)');
    LogK  = mean(LogKi);
    K = exp(LogK);
    kappa = exp(sqrt(mean(LogKi.^2 - LogK^2)));
    %text(2,33-(p-1)*6,['Geometric Mean: \it{K} = ' num2str(K,'%.2f')],...
    %     'color',color{p})
    %text(2,30-(p-1)*6,['Geometric Standard deviation: \kappa = ' ...
    %     num2str(kappa,'%.2f')],'color',color{p})
    text(2,33-(p-1)*4,['\it{K} = ' num2str(K,'%.2f') ...
        ', \kappa = ' num2str(kappa,'%.2f')],'color',color{p})
    title('Kamaishi Runup/Inundation')
    %legend('2DH','2CLOWNS-3D')
    %xlim([min(height_order) 12])
    %ylim([min(height_order) 12])
end
% set(gca,'fontsize',7)
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 7],'PaperPositionMode','manual');
% print('-r600','-depsc','../Paper/Max_Inun_Loc.eps');


% fid = fopen('../Inputs/Survey.xy','w');
% B_x = B_x * 1d-6; B_y = B_y * 1d-6; 
% % Print out interpolated results
% for k = 1:length(m_kam)
%     fprintf(fid,'%12.6f  %12.6f %12.6f \n',B_x(k),B_y(k),m_kam(k)) ;
% end
% fclose(fid);

% % Plot the inundation area
% 
% legend_name = {'2DH NSWE','2CLOWNS-3D','Measured Maximum (GSI, 2011)'};
%            {'2DH NSW Submerged','2CLOWNS-3D Submerged'...
%               '2DH NSW Overtopping','2CLOWNS-3D Overtopping'};
% axes(ha(2));
% for p = 1:L
%       num = dlmread([dir_num{p} 'AREA5.dat'],'',hlines,0);
%       hold on
%       plot(num(:,1)/60,(num(:,2)-num(1,2))/1d6,...
%            'LineStyle',line{p},'Color',color{p})
% end
% plot([0,90],[7,7],'-.k')
% Setting limits...
% xlim([min(time/3600),max(time/3600)]);
% xlim([0,90]);
% ylim([-250 400]);
% xlabel('time since earthquake [min]') %,'FontSize',8,'interpreter','latex');
% ylabel('Inundation Area [km^2]')
% ylabel('$Q$ [m$^2$s$^{-1}$]','FontSize',8,'interpreter','latex')%[m$^2$s$^{-1}$]
% legend(legend_name,... 
%        'location','SouthEast') %,'FontSize',8,'interpreter','latex');
% legend boxoff
% set(gca,'FontSize',8)
% [LX,LY] = gca_latex(8,1,1);