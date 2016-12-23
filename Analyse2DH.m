% Script to analyse the 2DH momentum fluxes over domain
% By: William Pringle Oct 31 2016

clearvars; 
clc; close all;
rho = 1003;
g = 9.807;
% Place for timeseries
x_want = 408150; %[407800 408500]; %[408000, 409080];
y_want = 845010:10:845300; %844700:10:845620; %843970:10:846190;
%% Set and load the data
for n =1:2
    if n == 1
        direc = 'Satakev8.0_3D_Real_BW/';
        fname = 'Satakev8.0_3D_Real_BW_';
        
        fname3 = 'Satake_BW_Real_ke.';
%        oname = 'BW_3D';
%         direc = 'Satakev8.0_3D_Nobw/';
%         fname = 'Satakev8.0_3D_NoBW_';
%         oname = '';  
        load([direc 'MAT_files/' fname3 'xyznF.mat'])
        x3 = x; y3 = y; mn3 = mn; z3 = z; is3 = is; ie3 = ie;
        load([direc 'MAT_files/' fname3 'VelF.mat'])
        [t3,I] = sort(time); Vel3 = Vel3(I,:,:); F3 = F3(I,:);
        kss = is(3) ; kee = ie(3)-1;
        load([direc 'MAT_files/' fname3 'P.mat'])
        P3 = P3(I,:);
    else
        direc = 'Satakev8.0_2DH_BW/';
        fname = 'Satakev8.0_2DH_BW_';
        %direc = 'Satakev8.0_2DH_Nobw/';
        %fname = 'Satakev8.0_2DH_';
        %oname = '';
    end
    load([direc 'MAT_files/' fname 'xyzF.mat']);  
    is2 = find(x == x3(is3(1))); ie2 = find(x == x3(ie3(1)));
    js2 = find(y == y3(is3(2))); je2 = find(y == y3(ie3(2)));
    load([direc 'MAT_files/' fname 'VelEta.mat'])  
    [time,I] = sort(time); Vel = Vel(I,:,:); Eta = Eta(I,:);
    ts3 = find(time == t3(2));
    %% Choose plane and wanted cross-section
    % Make desired target into grid format for x-y plane
    xx = 0.5*(x(is(1):ie(1)-1) + x(is(1)+1:ie(1)));
    yy = 0.5*(y(is(2):ie(2)-1) + y(is(2)+1:ie(2))) + 3.5d6;
    xx = xx/1d6;
    yy = yy/1d6;
    
%     %if n == 2
%     	Zk = z(4)-(z(4) - z(3))*Zk;
%     %end

    js = is(2)  ; je = ie(2)-1;
    iss = is(1) ; iee = ie(1)-1;

    QX = zeros(length(time),length(y),length(x));
    DP = zeros(length(time),length(y),length(x));
    GX = zeros(length(time),length(y),length(x));
    QU = zeros(length(time),length(yy),length(xx));
    EP = zeros(length(time),length(y),length(x));
    EK = zeros(length(time),length(y),length(x));
    PR = zeros(length(time),length(y),length(x));
    % For getting correct indices
    ix = 1; jy = 1; kz = 1;
    %
    for j = js:je
        for i = iss:iee
            if nf(j,i) == - 1; continue; end
            nn = mn(j,i); 
            if nn ~=0
                % Make QX
                nnmi = mn(j,i-ix); 
                if nnmi ~=0
                    %if n == 1 
                    DPx = max(0,0.5*(Eta(:,nnmi) + Eta(:,nn) ...
                                     - Zk(nn) - Zk(nnmi)));
                    %elseif n == 2
                    %    DPx = max(0,0.5*(Eta(:,nnmi) + Eta(:,nn)) ...
                    %                 - max(Zk(nn),Zk(nnmi)));
                    %end
                    Grad = (Eta(:,nn) - Eta(:,nnmi))/10;
                else
                    DPx = max(0,Eta(:,nn) - Zk(nn));
                    Grad = zeros(length(time),1);
                end
                GX(:,j,i) = Grad;
                QX(:,j,i) = Vel(:,nn,1).*DPx;
                DP(:,j,i) = max(0,Eta(:,nn) - Zk(nn));
                
                % Make the array
                %if strcmp(output,'u')
                nnpi = mn(j,i+ix); 
                nnpj = mn(j+jy,i);
                if nnpi ~=0
                    u1 = 0.5 * (Vel(:,nn,1) + Vel(:,nnpi,1)); 
                else
                    u1 = 0.5 * Vel(:,nn,1); 
                end
                if nnpj ~=0
                    u2 = 0.5 * (Vel(:,nn,2) + Vel(:,nnpj,2));
                else
                    u2 = 0.5 * Vel(:,nn,2); 
                end 
                QU(:,j-2,i-2) = (u1.^2 + u2.^2).*DP(:,j,i);
                if n == 1 && i >= is2 && i < ie2 && j >= js2 && j < je2
                    i3 = i - is2 + 3; j3 = j - js2 + 3;
                    for k = kss:kee
                        nn3 = mn3(j3,k,i3);
                        if nn3 == 0 || nn3 > length(fb); continue; end
                        nnpi3 = mn3(j3,k,i3+ix); 
                        nnpj3 = mn3(j3+jy,k,i3);
                        nnpk3 = mn3(j3,k+kz,i3);
                        if nnpi3 ~=0
                            u1 = 0.5 * (Vel3(2:end,nn3,1) + Vel3(2:end,nnpi3,1)); 
                        else
                            u1 = 0.5 * Vel3(2:end,nn3,1);
                        end
                        if nnpj3 ~=0
                            u2 = 0.5 * (Vel3(2:end,nn3,2) + Vel3(2:end,nnpj3,2)); 
                        else
                            u2 = 0.5 * Vel3(2:end,nn3,2);
                        end
                        if nnpk3 ~=0
                            u3 = 0.5 * (Vel3(2:end,nn3,3) + Vel3(2:end,nnpk3,3)); 
                        else
                            u3 = 0.5 * Vel3(2:end,nn3,3);
                        end  %+ u3.^2
                        EK(ts3:end,j,i) = (u1.^2 + u2.^2).*u1.... 
                     .*F3(2:end,nn3)*fb(nn3)*(z3(k+1)-z3(k)) + EK(ts3:end,j,i);
                 
                        EP(ts3:end,j,i) = 0.5*(z3(k+1)+z3(k))*u1.*F3(2:end,nn3)*fb(nn3)...
                              *(z3(k+1)-z3(k)) + EP(ts3:end,j,i);
                        PR(ts3:end,j,i) = P3(2:end,nn3).*u1.*F3(2:end,nn3)*fb(nn3)...
                              *(z3(k+1)-z3(k)) + PR(ts3:end,j,i);
                    end
                    EK(:,j,i) = 0.5*rho*EK(:,j,i);
                    EP(:,j,i) = rho*g.*EP(:,j,i);
                    PR(:,j,i) = rho*PR(:,j,i);
                else
                    EP(:,j,i) = 0.5*rho*g.*u1.*(Eta(:,nn).^2-Zk(nn)^2);
                    EK(:,j,i) = 0.5*rho*(u1.^2 + u2.^2).*DP(:,j,i).*u1;
                    PR(:,j,i) = rho*g*u1.*(0.5*Eta(:,nn).^2 + ...
                                Eta(:,nn).*Zk(nn) + 0.5*Zk(nn)^2);
                end
            end
        end
    end
    QUmax = squeeze(max(QU));
    QUmax(QUmax == 0) = NaN;
    %
    figure;
    pcolor(xx,yy,QUmax)
    shading interp
    colormap(jet)
    xlim([0.402 xx(end)])
    caxis([0 2000])
    %title('Velocity')
    colorbar
    %print('-r600','-dpng',['Presentation/' oname '_' plane 'Vel.png'])

    %
    I = knnsearch(x',x_want');
    J = knnsearch(y',y_want');
    QX_want = QX(:,J,I);
    DP_want = DP(:,J,I);
    GX_want = GX(:,J,I);
    EP_want = EP(:,J,I);
    EK_want = EK(:,J,I);
    PR_want = PR(:,J,I);
    
    QXs = squeeze(sum(QX_want,2));
    DPs = squeeze(mean(DP_want,2));
    GXs = squeeze(mean(GX_want,2));
    EPs = squeeze(mean(EP_want,2));
    EKs = squeeze(mean(EK_want,2));
    PRs = squeeze(mean(PR_want,2));
    kkk = 2;
    for iii = 1:length(x_want)
        kkk = kkk + 1;
        figure(kkk);
        plot(time,QXs(:,iii));
        hold on
        if n == 2; legend('Real','2DH'); end
        kkk = kkk + 1;
        figure(kkk);
        plot(time,DPs(:,iii));
        hold on
        if n == 2; legend('Real','2DH'); end
        kkk = kkk + 1;
        figure(kkk);
        plot(time,GXs(:,iii));
        hold on
        if n == 2; legend('Real','2DH'); end
        kkk = kkk + 1;
        figure(kkk);
        plot(time,EPs(:,iii),time,EKs(:,iii),time,PRs(:,iii),...
             time,EPs(:,iii)+EKs(:,iii)+PRs(:,iii));
        hold on
        if n == 2; 
            legend('3D EP','3D EK','3D PR','3D Tot',...
                   '2DH EP','2DH EK','2DH PR','2DH Tot'); 
        end
    end
end