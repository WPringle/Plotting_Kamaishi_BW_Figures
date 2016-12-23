% Script to analyse the jet through the submerged section of the breakwater
% as computed by 2CLOWNS-3D simulation with time, for making plots for paper. 
% By: William Pringle Oct 14 2016, updated Dec 23 2016

clearvars; clc; close all;

%Setting up the subplot margins

%% Set and load the data
for turb = 1:2
    if turb == 1
        direc = 'Satakev8.0_ke/';
        fname = 'Tohoku_SatakeT_HYB_ke_k1mm_wall.';
        oname = 'ke_wall';
    elseif turb == 2
        direc = 'Satakev8.0_3D_Real_BW/';
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
    ux_all = zeros(L,ie(1)); c_mu_all = zeros(L,ie(1)); nu_t_all = zeros(L,ie(1)); 
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
        ux = zeros(ie(1),1); c_mu = zeros(ie(1),1); nu_t_n = zeros(ie(1),1);
        uz05 = zeros(ie(1),k(end));  kz05 = zeros(ie(1),k(end)); 
        vz05 = zeros(ie(1),k(end));  tz05 = zeros(ie(1),k(end));
        uy05 = zeros(ie(1),ie(2)-1); ky05 = zeros(ie(1),ie(2)-1); 
        vy05 = zeros(ie(1),ie(2)-1); 
        for ii = is(1):ie(1)
            nn = mn(j,k,ii);
            ux(ii) = max(-u(nn,1));
            c_mu_n = nu_t(nn).*eps(nn)./rk(nn).^2; 
            nu_t_n(ii) = mean(mean(nu_t(nn)));
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
        ux_all(t,:)  = ux; c_mu_all(t,:) = c_mu; nu_t_all(t,:) = nu_t_n;
        uz_all(t,:,:) = uz05; kz_all(t,:,:) = kz05; vz_all(t,:,:) = vz05; 
        uy_all(t,:,:) = uy05; ky_all(t,:,:) = ky05; vy_all(t,:,:) = vy05; 
    end

    [time,tI] = sort(time); ux_want = ux_want(tI,:);
    
    %% Plot max nu_t and min c_mu versus location downstream
    figure(2);
    subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.02], ...
                                          [0.09 0.01], [0.12 0.01]);
    subplot(2,1,1)
    nu_t_max = max(nu_t_all); 
    if turb == 1
        plot(xh(1:i),nu_t_max(1:i),'k--')
        hold on
    elseif turb == 2
        plot(xh(1:i),nu_t_max(1:i),'k-')
        xlim([0 78])
        %ylim([0.03 0.09])
        %xlabel('\itx/h')
        set(gca,'XTickLabel','');
        ylabel('\it\nu_t')
        legend('Standard \itk-\epsilon','Realizable \itk-\epsilon',...
               'Location','NorthWest');
        %set(hl,'Interpreter','latex')
        set(gca,'fontsize',7) %,'Yticklabel',0.03:0.01:0.09)
        
        % Plot min c_mu versus location downstream
        subplot(2,1,2)
        c_mu_min = min(c_mu_all); 
        plot(xh(1:i),c_mu_min(1:i),'-')
        xlabel('\itx/h')
        xlim([0 78])
        ylim([0.03 0.09])
        ylabel('\itC_{\mu}')
        set(gca,'fontsize',7,'Yticklabel',0.03:0.01:0.09)
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 9 9],...
            'PaperPositionMode','manual');
        print('-r600','-depsc','../Paper/nu_t_max_cu_mu_min.eps'); 
    end
    %% Time average ux_all for decay analysis
    figure(1); 
    subplot = @(m,n,p) subtightplot (m, n, p, [0.06 0.03], ...
                                          [0.08 0.02], [0.07 0.01]);
    %
    U90 = prctile(ux_all,90,1);
    Um = zeros(1,i);
    for ii = 1:i
        Um(ii)  = mean(ux_all(ux_all(:,ii) > U90(ii),ii));
    end
    U0 = Um(i);
    
    % Plot the decay of centreline velocity
    subplot(2,2,turb)
    loglog(xh(1:i),Um/U0,'.')
    % 
    hold on
    if strcmp(direc,'Satakev8.0_Norm_Real_ke/') || strcmp(direc,'Satakev8.0_3D_Real_BW/')
        %I = find(xh >= 1 & xh <= 15);
        I = find(xh >= 1 & xh <= 20);
        p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
        loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
        text(3,0.83 ,['[\itx/h\rm]^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)

        %I = find(xh >= 15 & xh <= 60);
        I = find(xh >= 18 & xh <= 50);
        p2 = polyfit(log(xh(I)),log(Um(I)/U0),1);
        loglog(xh(I),exp(p2(2))*xh(I).^p2(1),'-','LineWidth',1.25)
        text(19,0.75 ,['[\itx/h\rm]^{' num2str(p2(1),'%.2f') '}'],'fontsize',7)

        %I = find(xh >= 60 & xh <= 85);
        I = find(xh >= 48 & xh <= 66);
        p3 = polyfit(log(xh(I)),log(Um(I)/U0),1);
        loglog(xh(I),exp(p3(2))*xh(I).^p3(1),'-','LineWidth',1.25)
        text(26,0.60 ,['[\itx/h\rm]^{' num2str(p3(1),'%.2f') '}'],'fontsize',7)  
    elseif strcmp(direc,'Satakev8.0_ke/')
        I = find(xh >= 1 & xh <= 17);
        p1 = polyfit(log(xh(I)),log(Um(I)/U0),1);
        loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
        text(3,0.80 ,['[\itx/h\rm]^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)

        I = find(xh >= 15 & xh <= 35);
        p2 = polyfit(log(xh(I)),log(Um(I)/U0),1);
        loglog(xh(I),exp(p2(2))*xh(I).^p2(1),'-','LineWidth',1.25)
        text(15,0.77 ,['[\itx/h\rm]^{' num2str(p2(1),'%.2f') '}'],'fontsize',7)

        I = find(xh >= 32 & xh <= 40);
        p3 = polyfit(log(xh(I)),log(Um(I)/U0),1);
        loglog(xh(I),exp(p3(2))*xh(I).^p3(1),'-','LineWidth',1.25)
        text(40,0.75 ,['[\itx/h\rm]^{' num2str(p3(1),'%.2f') '}'],'fontsize',7)    
    end

    grid on 
    xlabel('\itx/h')
    if turb == 1
        ylabel('\itU_{m}/U_0')
    else
        set(gca,'YTickLabel','');
    end
    xlim([1 80])
    %xlim([4 80])
    ylim([0.1 1.0])
    ax = gca;
    ax.XTick = [1 2 3 4 6 8 10 20 40 60];
    ax.YTick = [0.1 0.2 0.4 0.6 0.8 1.0];
    set(gca,'fontsize',7)
    
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
    
    %Plotting the growth rates
    subplot(2,2,2+turb)
    loglog(xh(1:i),zhalf/h,'.')
    hold on
    loglog(xh(1:i),yhalf1/h,'.')
    loglog(xh(1:i),yhalf2/h,'.')
    if strcmp(direc,'Satakev8.0_Norm_Real_ke/') || strcmp(direc,'Satakev8.0_3D_Real_BW/')
        I = find(xh >= 2 & xh <= 6.5); I(isnan(zhalf(I))) = [];
        p1 = polyfit(log(xh(I)),log(zhalf(I)'/h),1);
        loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
        text(2.7,2.6 ,['[\itx/h\rm]^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)

%         I = find(xh >= 2 & xh <= 6.5); I(isnan(zhalf(I))) = [];
%         p1 = polyfit(log(xh(I)),log(zhalf(I)'/h),1);
%         loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
%         text(2.7,1 ,['[x/h]^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)
    elseif strcmp(direc,'Satakev8.0_ke/')
        I = find(xh >= 1 & xh <= 3.5);
        p1 = polyfit(log(xh(I)),log(zhalf(I)'/h),1);
        loglog(xh(I),exp(p1(2))*xh(I).^p1(1),'-','LineWidth',1.25)
        text(2.5,1.8,['[\itx/h]^{' num2str(p1(1),'%.2f') '}'],'fontsize',7)
    end
    legend('\itz_{1/2}/h','\ity^S_{1/2}/h','\ity^N_{1/2}/h','Location','best')

    grid on
    xlim([1 80])
    ylim([1 10])
    xlabel('\itx/h')
    if turb == 1
        ylabel('\itz_{1/2}/h,y_{1/2}/h')
    else
        set(gca,'YTickLabel','');
    end
    ax = gca;
    ax.XTick = [1 2 3 4 6 8 10 20 40 60];
    ax.YTick = [1 2 4 6 8 10];
    set(gca,'fontsize',7)
end
figure(1);
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 10],...
    'PaperPositionMode','manual');
print('-r600','-depsc','../Paper/Decay_growth.eps'); 
