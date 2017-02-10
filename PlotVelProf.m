%%%%%% Plot the profile of velocities and pressure %%%%%%%%%%%%%%%%%%%%%%
clearvars;
close all;
clc;

g = 9.807;

% directory and filename
direc = '../Satakev8.0_3D_Real_BW/';
fname = 'Satake_BW_Real_ke.';
% y cross-section point
yp = 845140;

% the times to get plots for
time = [1899.0, 2427.1];

% Setting up the subplot margins
% subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], ...
%                                           [0.08 0.03], [0.045 0.01]);

subplot = @(m,n,p) subtightplot (m, n, p, [0.03 0.01], ...
                                          [0.05 0.02], [0.045 0.01]);
figure;
ty = 0;
for plottype = {'whole','closeup'}
% plottype
%plottype = 'closeup';
    if strcmp(plottype,'closeup')
        %Positions to measure
        xp = [408000 408050 408100 408150 408200 408250 408300];
        % multiplier for ruler 
        pm = [4,10,0.5];
        %
        diff = 40;
        mult = 0.4;
    elseif strcmp(plottype,'whole')
        %Positions to measure
        xp = [406450 406875 407300 407725 408150 408575 409000];
        % multiplier for ruler 
        pm = [25,50,5];
        %
        diff = 250;
        mult = 0.6;
    end

    %% Load data and produce plots
    % Loading the xyz and fb data
    load([direc 'MAT_files/' fname 'xyznF.mat']);  

    zz = 0.5*(z(is(3):ie(3)-1) + z(is(3)+1:ie(3)));
    % Get index of require y value
    I = knnsearch(x',xp');
    J = find(y == yp);
    % Get ground level
    zg = zeros(ie(1)+1,1);
    ii = 0;
    for i = is(1):ie(1)
        ii = ii + 1;
        for k = is(3):ie(3)-1
           if mn(J,k,i) > 0 
               if a(mn(J,k,i),1) > 0
                  zg(i) = z(k+1) - a(mn(J,k,i),1)*(z(k+1)-z(k));
                  break;
               end
           end
        end
    end
    [~,mI] = max(zg(zg < 0));
    x = [x(1:mI+1),x(mI+2),x(mI+2:mI+3),x(mI+3),x(mI+4:end)];
    zg = [zg(1:mI+1);zg(mI+1);zg(mI+2:mI+3);zg(mI+4);zg(mI+4:end)]; 
    %
    Ip = knnsearch(x',xp');
    Ip(xp == 408150) = Ip(xp == 408150) + 1;
    % Loop over the times creating new figures
    x = x*1e-6;
    xp = xp*1e-6;
    pm = pm*1e-6;
    tt = 2;
    %figure;
    for t = time
        tt = tt - 1;
        % Load current time
        load([direc 'MAT_files/' fname 't' num2str(t,'%.1f') '.mat'])

        for i = 1:length(xp)
            for k = is(3):ie(3)-1
               if mn(J,k,I(i)) > 0 
                   if F(mn(J,k,I(i))) > 0
                      eta(i) = z(k) + F(mn(J,k,I(i)))*fb(mn(J,k,I(i)))...
                                                       *(z(k+1)-z(k));
                   end
               end   
            end
            for k = is(3):ie(3)-1
               if mn(J,k,I(i)) > 0 
                   Up(i,k) = u(mn(J,k,I(i)),1);
                   Wp(i,k) = u(mn(J,k,I(i)),3);
                   Pp(i,k) = p(mn(J,k,I(i))) - g*(eta(i)-(z(k)+z(k+1))*0.5);
               end
            end
        end
        Pp(Up == 0) = NaN;
        Up(Up == 0) = NaN;
        Wp(Wp == 0) = NaN;

        % Loop over plot types (velocity and pressure)
        for p = 1:3
            %subplot(3,2,2*p-tt)
            subplot(6,2,(2*p-tt) + ty)
            % Draw the ground
            area(x,zg,-80,'FaceColor',[101,77,33]/255,'FaceAlpha',0.7)
            hold on
            % Drawing the middle positions
            if strcmp(plottype,'whole')
                xs = x(is(1));
                xlim([xs x(ie(1))])
            else
                xs = xp(1)-50*1d-6;
                xlim([xs xp(end)+50*1d-6])
            end
            for i = 1:length(xp)
                if p == 1
                    % Draw the middle positions
                    plot([xp(i) xp(i)],[zg(Ip(i)) 20],'--k')
                    % Write text of current time
                    if i == 1 && strcmp(plottype,'whole')
                        title(['\itt \rm= ' num2str(t/60,3) ' min'])
                    end
                    %text(xp(end)+diff*0.04*1e-6,8,['\itt \rm= ' num2str(t/60,3) ' min'],'fontsize',8)
                else
                    % Draw the middle positions
                    plot([xp(i)+5*1e-6 xp(i)+5*1e-6],[zg(Ip(i)) 20],'--k')
                end
            end
            ylim([-80 20])
            % get handle to current axes
            %set(gca,'Ticklength',[0,5])
            %
            for i = 1:length(xp)
                if i == 1
                    % Make ruler
                    plot([xs+(xp(1)-xs)*1/8 ...
                          xs+(xp(1)-xs)*1/8 + diff*1e-6],[0,0],'-k')
                    plot([xs+(xp(1)-xs)*1/8 ...
                          xs+(xp(1)-xs)*1/8],[0,4],'-k')
                    plot([xs+(xp(1)-xs)*1/8 + diff*1e-6 ...
                          xs+(xp(1)-xs)*1/8 + diff*1e-6],[0,4],'-k')
                    plot([xs+(xp(1)-xs)*1/8 + diff/2*1e-6 ...
                          xs+(xp(1)-xs)*1/8 + diff/2*1e-6],[0,4],'-k')
                    if tt == 1
                        ylabel('z [m]')
                    else
                        set(gca,'YTickLabel','');
                    end
                    text(xs+(xp(1)-xs)*1/8-0.06*diff*1e-6,10,num2str(0),'fontsize',7)
                    text(xs+(xp(1)-xs)*1/8-0.14*diff*1e-6 + diff*1e-6,10,num2str(diff*1e-6/pm(p)),'fontsize',7)
                    text(xs+(xp(1)-xs)*1/8-0.14*diff*1e-6 + diff/2*1e-6,10,num2str(diff/2*1e-6/pm(p)),'fontsize',7)
                end
                if p == 1
                    % Draw the velocity profiles
                    plot(Up(i,is(3):ie(3)-1)*pm(p)+xp(i),zz,'-k')  
                    %set(gca,'XTickLabel','');
                    if i == 1
                        text(xs+(xp(1)-xs)*1/8-diff*mult*1e-6 + diff/2*1e-6,...
                             -7,'\itu \rm[m/s]','fontsize',7)
                    end
                elseif p == 2
                    % Draw the velocity profiles
                    plot(Wp(i,is(3):ie(3)-1)*pm(p)+xp(i)+5*1e-6,zz,'-k')
                    %set(gca,'XTickLabel','');
                    if i == 1
                    text(xs+(xp(1)-xs)*1/8-diff*mult*1e-6 + diff/2*1e-6,...
                         -7,'\itw \rm[m/s]','fontsize',7)
                    end
                elseif p == 3
                    % Draw the pressure profiles
                    plot(Pp(i,is(3):ie(3)-1)*pm(p)+xp(i)+5*1e-6,zz,'-k')
                    if i == 1
                    text(xs+(xp(1)-xs)*1/8-diff*mult*1e-6 + diff/2*1e-6,...
                         -7,'\itp_d \rm[kPa]','fontsize',7)
                    end
                    if strcmp(plottype,'closeup')
                        xlabel('Easting [10^6 m]')
                    end
                end
                set(gca, 'XTick', xp)
            end
            set(gca,'fontsize',7)
            % get handle to current axes
            a = gca;
            % set box property to off and remove background color
            set(a,'box','off','color','none')
            % create new, empty axes with box but without ticks
            b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
            % set original axes as active
            axes(a)
            % link axes in case of zooming
            linkaxes([a b])
        end

    end
    ty = ty + 6;

end
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 19 17],...
    'PaperPositionMode','manual');
print('-r600','-depsc',['../../Paper/Vel_Prof_all.eps']);
