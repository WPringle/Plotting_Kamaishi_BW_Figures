% Script to analyse the jet through the submerged section of the breakwater
% as computed by 2CLOWNS-3D simulation with time, for making plots for paper. 
% By: William Pringle Oct 14 2016, updated Dec 23 2016

clearvars; clc; close all;

%percentile setting
perc_set = 80;

%% Set and load the data
for turb = 1:2
    if turb == 1
        direc = '../Satakev8.0_ke/';
        fname = 'Tohoku_SatakeT_HYB_ke_k1mm_wall.';
        oname = 'ke_wall';
    elseif turb == 2
        direc = '../Satakev8.0_3D_Real_BW/';
        fname = 'Satake_BW_Real_ke.';
        oname = 'Real_ke';
    end

    fPath = cd;
    List = dir(fullfile(fPath, [direc fname '*']));
    List([List.isdir]) = [];
    L = length(List)/4;
    % order list by datenum
    Listc = struct2cell(List);
    Ordered = Listc'; %sortrows(Listc',5);

    load([direc 'MAT_files/' fname 'xyzn.mat']);  

    %% Choose wanted points 
    if turb == 2
        L = floor(L);
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
    nu_t_all = zeros(L,ie(1)); k_all = zeros(L,ie(1)); 
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
        nu_t_n = zeros(ie(1),1); k_n = zeros(ie(1),1);
        uz05 = zeros(ie(1),k(end));  kz05 = zeros(ie(1),k(end)); 
        vz05 = zeros(ie(1),k(end));  tz05 = zeros(ie(1),k(end));
        uy05 = zeros(ie(1),ie(2)-1); ky05 = zeros(ie(1),ie(2)-1); 
        vy05 = zeros(ie(1),ie(2)-1); 
        for ii = is(1):ie(1)
            nn = mn(j,k,ii);
            ux(ii) = max(-u(nn,1));
            c_mu_n = nu_t(nn).*eps(nn)./rk(nn).^2; 
            nu_t_n(ii) = mean(mean(nu_t(nn)));
            %k_n(ii) = mean(mean(rk(nn)));
            k_n(ii) = mean(sqrt(2/3*reshape(rk(nn),...
                               size(rk(nn),1)*size(rk(nn),2),1))./...
                               sqrt(u(nn,1).^2 + u(nn,2).^2 + u(nn,3).^2));
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
        nu_t_all(t,:) = nu_t_n; k_all(t,:) = k_n;
        uz_all(t,:,:) = uz05; kz_all(t,:,:) = kz05; vz_all(t,:,:) = vz05; 
        uy_all(t,:,:) = uy05; ky_all(t,:,:) = ky05; vy_all(t,:,:) = vy05; 
    end

    [time,tI] = sort(time); ux_want = ux_want(tI,:);
    
    %% Plot max nu_t and min c_mu versus location downstream
    nu_t_max = max(nu_t_all); 
    k_max = max(k_all);
    if turb == 1
        %[ax, h1, h2] = plotyy(xh(1:i),nu_t_max(1:i),xh(1:i),k_max(1:i)); %,'k--')
        % save
        xh_1 = xh(1:i);
        nut_1 = nu_t_max(1:i);
        k_1 = k_max(1:i);
    elseif turb == 2
        figure(2);
        subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.02], ...
                                              [0.1 0.01], [0.12 0.05]);
        subplot(2,1,1)
        
        [ax, h1, h2] = plotyy(xh_1,nut_1,xh_1,k_1); %,'k--')
        set(ax(1), 'YColor', 'blue')
        set(ax(2), 'YColor', 'red')
        set(h1, 'Color', 'blue')
        set(h2, 'Color', 'red')
        hold(ax(1), 'on')
        hold(ax(2), 'on')
        h3 = plot(xh(1:i),nu_t_max(1:i), '--b', 'Parent', ax(1));
        h4 = plot(xh(1:i),k_max(1:i), '--r', 'Parent', ax(2)); %,'k-')
        %xlim([0 78])
        %ylim([0.03 0.09])
        %xlabel('\itx/h')
        set(gca,'XTickLabel','');
        axes(ax(1)); ylabel('\it\nu_t\rm [m^2s^{-1}]'); xlim([0 78])
        axes(ax(2)); ylabel('\itI'); xlim([0 78])
        %ylabel('\it\nu_t\rm [m^2s^{-1}]')
        
        legend([h1 h2 h3 h4], ...
            '\nu_t (standard)','I (standard)',...
            '\nu_t (realizable)','I (realizable)',...
             'Location','NorthWest');
        
%         legend('Standard \itk-\epsilon','Realizable \itk-\epsilon',...
%                'Location','NorthWest');
        %set(hl,'Interpreter','latex')
        set(ax(1),'fontsize',7) %,'Yticklabel',0.03:0.01:0.09)
        set(ax(2),'fontsize',7)
        
        % Plot min c_mu versus location downstream
        subplot(2,1,2)
        c_mu_min = min(c_mu_all); 
        plot(xh(1:i),c_mu_min(1:i),'k-')
        xlabel('\itx_j/h_j')
        xlim([0 78])
        ylim([0.03 0.09])
        ylabel('\itC_{\mu}')
        set(gca,'fontsize',7,'Yticklabel',0.03:0.01:0.09)

       %% Plot points versus time 
        figure(3);
        plot(time/60,ux_want,'x-')
        hold on
        for jj = 1:length(xh_want)
            idx = find(ux_want(:,jj) > prctile(ux_want(:,jj),80));
            plot(time(idx)/60,ux_want(idx,jj),'ko')
        end
        %plot([21.5 40],[ prctile(ux_want,80); prctile(ux_want,80)],'k--')
        legend(['\itx_j/h_j\rm = ' num2str(xh_want(1))],...
               ['\itx_j/h_j\rm = ' num2str(xh_want(2))],...
               ['\itx_j/h_j\rm = ' num2str(xh_want(3))],...
               ['\itx_j/h_j\rm = ' num2str(xh_want(4))],...
               'location','best') %['x/h = ' num2str(xh_want(5))]
        xlim([21.67 40])
        ylim([0 14])
        xlabel('time since earthquake rupture [min]')
        ylabel('\itU_{m}\rm [m/s]')
        set(gca,'fontsize',7)
    end
    
    %% Time average ux_all for decay analysis
    figure(1); 
    subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.03], ...
                                          [0.1 0.02], [0.07 0.01]);
    %
    U90 = prctile(ux_all,perc_set,1);
    Um = zeros(1,i); interval = zeros(length(time),i);
    for ii = 1:i
        intn = find(ux_all(:,ii) > U90(ii));
        if ~isempty(intn)
            Um(ii)  = mean(ux_all(intn,ii));
            interval(1:length(intn),ii) = intn;
        end
    end
    U0 = Um(i);
    
    % Plot the decay of centreline velocity
    subplot(2,2,turb)
    loglog(xh(1:i),Um/U0,'.')
    % 
    hold on
    II = find(Um/U0 >= 0.4 & Um/U0 <= 0.6);
    slope_p = zeros(size(II));
    for nn = 1:length(II)
        I = find(xh(1:length(Um)) >= 1 & Um/U0 > Um(II(nn))/U0);
        % Finding lines using SLM engine to get arbitrary knots
        slm = slmengine(log(xh(I)),log(Um(I)/U0),...
          'degree',1,'interiorknots','free','knots',4);
        I = find(xh >= exp(slm.knots(3)) & xh <= exp(slm.knots(4)));
        p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
        slope_p(nn) = p1(1);
    end
    slope_vary(turb) = mean(slope_p);
    slope_std(turb) = std(slope_p);

    I = find(xh(1:length(Um)) >= 1 & Um/U0 > 0.5);
    % Finding lines using SLM engine to get arbitrary knots
    slm = slmengine(log(xh(I)),log(Um(I)/U0),...
          'degree',1,'interiorknots','free','knots',4);
    for jj = 1:3
        I = find(xh >= exp(slm.knots(jj)) & xh <= exp(slm.knots(jj+1)));
        p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
        loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
        if jj == 1
            text(exp(mean(log(xh(I)))),0.8,...
             ['\it\proptox_j\rm^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)
        elseif jj == 2
            text(exp(mean(log(xh(I)))),exp(max(log(Um(I)/U0))),...
             ['\it\proptox_j\rm^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)      
        else
            text(exp(mean(log(xh(xh >= exp(slm.knots(2)) & ...
                 xh <= exp(slm.knots(3))))))-3,exp(min(log(Um(I)/U0))),...
             ['\it\proptox_j\rm^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)   
        end
    end
    grid on 
    %xlabel('\itx_j/h_j')
    if turb == 1
        title('Standard \itk-\epsilon')
        ylabel('\itU_{m}/U_0')
    else
        set(gca,'YTickLabel','');
        title('Realizable \itk-\epsilon')
    end
    xlim([1 80])
    %xlim([4 80])
    ylim([0.1 1.0])
    ax = gca;
    ax.XTick = [1 2 3 4 6 8 10 20 40 60];
    ax.YTick = [0.1 0.2 0.4 0.6 0.8 1.0];
    set(gca,'fontsize',7)
    
    %% Look at z and y growth rates 
    Uz90 = squeeze(prctile(uz_all,perc_set,1)); 
    Kz90 = squeeze(prctile(kz_all,perc_set,1));
    Vz90 = squeeze(prctile(vz_all,perc_set,1)); 
    Tz90 = squeeze(prctile(tz_all,perc_set,1));
    Uz05 = zeros(i,k(end)); Kz05 = zeros(i,k(end)); 
    Vz05 = zeros(i,k(end)); Tz05 = zeros(i,k(end));
    zhalf = zeros(i,1);
    Uy90 = squeeze(prctile(uy_all,perc_set,1)); 
    Ky90 = squeeze(prctile(ky_all,perc_set,1));
    Vy90 = squeeze(prctile(vy_all,perc_set,1)); 
    Ty90 = squeeze(prctile(ty_all,perc_set,1));
    Uy05 = zeros(i,ie(2)-1); Ky05 = zeros(i,ie(2)-1); 
    Vy05 = zeros(i,ie(2)-1); Ty05 = zeros(i,ie(2)-1);
    yhalf1 = zeros(i,1);
    yhalf2 = zeros(i,1);
    for ii = i:-1:1
        % z growth rate
        for kk = 1:k(end)
%             Uz05(ii,kk)  = mean(uz_all(uz_all(:,ii,kk) > Uz90(ii,kk),ii,kk));
%             Kz05(ii,kk)  = mean(kz_all(kz_all(:,ii,kk) > Kz90(ii,kk),ii,kk));
%             Vz05(ii,kk)  = mean(vz_all(vz_all(:,ii,kk) > Vz90(ii,kk),ii,kk));
%             Tz05(ii,kk)  = mean(tz_all(vz_all(:,ii,kk) > Tz90(ii,kk),ii,kk));
            intn = interval(:,ii); intn(intn == 0) = [];          
            Uz05(ii,kk)  = mean(uz_all(intn,ii,kk));
            Kz05(ii,kk)  = mean(kz_all(intn,ii,kk));
            Vz05(ii,kk)  = mean(vz_all(intn,ii,kk));
            Tz05(ii,kk)  = mean(tz_all(intn,ii,kk));
        end
        I = find(Uz05(ii,1:k(end)-2) <= 0.5*Um(ii)); 
        if ~isempty(I)
            c1 = 0.5*Um(ii) - Uz05(ii,I(end));
            c2 = Uz05(ii,I(end)+1) - 0.5*Um(ii);
            c3 = Uz05(ii,I(end)+1) - Uz05(ii,I(end));
            zhalf(ii) = z(I(end)) * c2/c3 + z(I(end)+1) * c1/c3;
        end
        % y growth rate
        for jj = is(2):ie(2)-1
%             Uy05(ii,jj)  = mean(uy_all(uy_all(:,ii,jj) > Uy90(ii,jj),ii,jj));
%             Ky05(ii,jj)  = mean(ky_all(ky_all(:,ii,jj) > Ky90(ii,jj),ii,jj));
%             Vy05(ii,jj)  = mean(vy_all(vy_all(:,ii,jj) > Vy90(ii,jj),ii,jj));
%             Ty05(ii,jj)  = mean(ty_all(ty_all(:,ii,jj) > Ty90(ii,jj),ii,jj));
            
            intn = interval(:,ii); intn(intn == 0) = [];          
            Uy05(ii,jj)  = mean(uy_all(intn,ii,jj));
            Ky05(ii,jj)  = mean(ky_all(intn,ii,jj));
            Vy05(ii,jj)  = mean(vy_all(intn,ii,jj));
            Ty05(ii,jj)  = mean(ty_all(intn,ii,jj));
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
    
    %Plotting the growth rates
    subplot(2,2,2+turb)
    loglog(xh(1:i),zhalf/h,'.')
    hold on
    loglog(xh(1:i),yhalf1/h,'.')
    loglog(xh(1:i),yhalf2/h,'.')
%     % Plot best fit
%     II = find(zhalf/h > 2.6); I = find(xh(II(end):end) >= 1); 
%     I = II(end) + I -1;
%     p1 = polyfit(log(xh(I)),log(zhalf(I)'/h),1);
%     loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
%     text(2.5,1.8,['\it\proptox_j\rm^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)

    I = find(xh > exp(slm.knots(2)) & xh < exp(slm.knots(end)));
    p1 = polyfit(log(xh(I)),log(yhalf2(I)'/h),1);
    loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
    text(min(xh(I)),mean(yhalf2(I)/h),...
         ['\it\proptox_j\rm^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)
%     if turb == 2
%         II = find(zhalf/h > 2.6); I = find(xh(II(end):end) >= 2); 
%         I = II(end) + I -1;
%         p1 = polyfit(log(xh(I)),log(zhalf(I)'/h),1);
%         loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
%         text(2.7,2.6 ,['\it\proptox_j\rm^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)
%     elseif turb == 1
% 
%     end
    legend('\itz_{1/2}/h_j','\ity^S_{1/2}/h_j','\ity^N_{1/2}/h_j',...
           'Location','North')

    grid on
    xlim([1 80])
    ylim([1 10])
    xlabel('\itx_j/h_j')
    if turb == 1
        ylabel('\itz_{1/2}/h_j,y_{1/2}/h_j')
    else
        set(gca,'YTickLabel','');
    end
    ax = gca;
    ax.XTick = [1 2 3 4 6 8 10 20 40 60];
    ax.YTick = [1 2 4 6 8 10];
    set(gca,'fontsize',7)
    
%     % Plot estimate nu_t
%     figure(2);
%     subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.02], ...
%                                           [0.1 0.01], [0.12 0.01]);
%     subplot(2,1,1)
%     hold on
%     plot(xh(1:i),0.016*Um(1:i).*yhalf2')
end
% figure(1);
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 10],...
%     'PaperPositionMode','manual');
% print('-r600','-depsc','../../Paper/Decay_growth.eps'); 
% 
% figure(2);
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 9 9],...
%     'PaperPositionMode','manual');
% print('-r600','-depsc','../../Paper/nu_t_max_cu_mu_min.eps'); 
% 
% figure(3);
% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 10 7],...
%     'PaperPositionMode','manual');
% print('-r600','-depsc','../../Paper/transient_u_max.eps'); 