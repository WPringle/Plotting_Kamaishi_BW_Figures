%%%%%% Plot the map of free-surface heights, velocities etc %%%%%%%%%%%%%%
clear all;
close all;
k = 0;
%Get relevant simulation results
dir_num{1} = 'Tohoku_HYB_10m_Fujima_CaseC_z2.2_BWk1mm_no_lim2D.t7200.9MAX_5.dat';
%use the function GetMaxData2015 to extract the vectors and variable
L = length(dir_num);

for Var = 3:3 %:4 %2:3
%Var = 1; %ETA
%Var = 2; %DEP
%Var = 3; %VEL
%Var = 4; %VOR
if Var == 1
    word = 'ETA';
    %c_min = -2;
    %c_max = 2;
    c_min = 0;
    c_max = 20;
elseif Var == 2
%    word = 'VEL';
    c_min = 0;
    c_max = 20;
%     word = 'DEP';
%     c_min = -20;
%     c_max = 100; 
elseif Var <= 3 
%    word = 'DIF';
%    c_min = -1;
%    c_max = 1;
     word = 'VEL';
     c_min = 0;
     c_max = 20;
elseif Var == 4
    word = 'VOR';
    c_min = -0.2;
    c_max = 0.2;
elseif Var == 5
    word = 'MOM';
    c_min = 0;
    c_max = 5000;
end

if Var == 1
    %Get relevant survey data
    load('Survey_data_order.mat'); %Contains x, y and z data
    %Delete starting and ending areas (of extremely large inundations)
    I = [1 2 3 4 82 83]'; %8 
    height_order(I) = []; Xp_order(I) = [];  Yp_order(I) = [];
    %Delete locations where survey is unreliable
    I = find(isnan(height_order));
    height_order(I) = []; Xp_order(I) = [];  Yp_order(I) = [];
    %load No_inun_data.mat
    %height_order(K) = []; Xp_order(K) = [];  Yp_order(K) = [];
    N = length(height_order);
end

for i = 1:L
    k = k+1;
    %subplot(2,1,k)
    %axes(ha(k));
    if Var == 5
        [ xx, yy, Value ] = GetMaxData2015( dir_num{i}, 3 );
        [ ~, ~, Value1 ] = GetMaxData2015( dir_num{i}, 2 );   
        h = pcolor(xx,yy+3.5d6,Value.^2.*Value1); 
    elseif Var == 3
        [ xx, yy, Value ] = GetMaxData2015( dir_num{i}, Var );
        h = pcolor(xx,yy+3.5d6,Value); 
    elseif Var == 1
        [ xx, yy, Value ] = GetMaxData2015( dir_num{i}, Var );
        %[ ~, ~, Value1 ] = GetMaxData2015( dir_num{i-1}, Var );
        %h = pcolor(xx,yy+3.5d6,Value1-Value); 
        h = pcolor(xx,yy+3.5d6,Value); 
        %[ ~, ~, U2{i} ] = GetMaxData2015( dir_num{i}, 5 );
        %[ ~, ~, V2{i} ] = GetMaxData2015( dir_num{i}, 6 );
%         [ xx, yy, zz, U{i} ] = GetData3D2015( dir_num{i}, 3 );
%         [ ~, ~, ~, V{i} ] = GetData3D2015( dir_num{i}, 4 );
%         zs = find ( zz == -27);
%         U2{i} = flipud(imrotate(squeeze(U{i}(:,:,zs)),90));
%         V2{i} = flipud(imrotate(squeeze(V{i}(:,:,zs)),90));
        %h = pcolor(xx,yy+3.5d6,sqrt(U2{i}.^2+V2{i}.^2)); 
%     elseif Var == 3
%         MaxV = sqrt(U.^2 + V.^2);
%         [MaxV, I] = max(MaxV,[],1);
%         MaxV = squeeze(MaxV); 
%         MaxU = zeros(size(MaxV)); MaxW = zeros(size(MaxV));
%         iL = length(xx); jL = length(yy); 
%         for i = 1:iL
%             for j = 1:jL
%                 MaxU(i,j) = U(I(1,i,j),i,j);
%                 MaxW(i,j) = V(I(1,i,j),i,j);
%             end
%         end
%         [qx,qy] = meshgrid(xx1,yy1+3.5d6); s=size(qx);
%         x = reshape(qx,s(1)*s(2),1);
%         y = reshape(qy,s(1)*s(2),1);
%         z = reshape(eta,s(1)*s(2),1);
%     %get the interpolant
%         F = TriScatteredInterp(x,y,z,'linear'); % for chikei;
%         [qx,qy] = meshgrid(xx,yy);
%     % %find the mesh values
%         qz = F(qx,qy);
%         qz1 = flipud(imrotate(MaxV,90));
%         %h = pcolor(xx,yy,flipud(imrotate(MaxV,90)));
%         h = pcolor(xx,yy,(qz1-qz)./qz);
    elseif Var == 4
         S = size(U2{i}); Vor = NaN(S); 
         for j = 1:S(2)-1
             for z = 1:S(1)-1
                 Vor(z,j) = (U2{i}(z+1,j)-U2{i}(z,j))/(yy(z+1)-yy(z)) + ...
                              (V2{i}(z,j+1)-V2{i}(z,j))/(xx(j+1)-xx(j));
             end
         end
%         MaxVor = squeeze(max(Vor,[],1));
%         MinVor = squeeze(min(Vor,[],1));
%         VorM = NaN(size(MaxVor));
%         for i = 1:iL-1
%             for j = 1:jL-1
%                 if abs(MaxVor(i,j)) >= abs(MinVor(i,j))
%                     VorM(i,j) = MaxVor(i,j);
%                 else
%                     VorM(i,j) = MinVor(i,j);
%                 end
%             end
%         end
        h = pcolor(xx,yy+3.5d6,Vor);
        %h = pcolor(xx,yy,flipud(imrotate(VorM,90)));
    end

%Umax = squeeze(max(abs(U),[],1)); Vmax = squeeze(max(abs(V),[],1)); 
%eta = sqrt(Umax.^2+Vmax.^2);
%h = pcolor(xx,yy,flipud(imrotate(MaxV,90)));
set(h, 'EdgeColor', 'none');
caxis([c_min c_max])
hold on
% S = size(U2{i});
% [x,y] = meshgrid(xx,yy);
% gapy = 5;
% gapx = 20;
% I = find(U2{i} < -10); U2{i}(I) = -10;
% q = quiver(x(1:gapy:S(1),1:gapx:S(2)),y(1:gapy:S(1),1:gapx:S(2))+3.5d6,...
%     U2{i}(1:gapy:S(1),1:gapx:S(2)),V2{i}(1:gapy:S(1),1:gapx:S(2)),2.5);
%      %flipud(imrotate(U2(1:3:S(1),1:3:S(2)),90)),...
%      %flipud(imrotate(V2(1:3:S(1),1:3:S(2)),90)),1.5);
%  set(q,'color','black')
%plot the survey results
if Var == 1
scatter( Xp_order, Yp_order, 14, height_order, ...
        'o', 'filled', 'MarkerEdgeColor','k'  );
end
axis equal
%axis([406000 409500 4.3435d6 4.3465d6])
axis([402000 max(xx) min(yy)+3.5d6 max(yy)+3.5d6])
%rectangle('Position', [407350 844660+3.5d6 408600-407350 846000-844660]...
%          ,'LineStyle','--')

%     if k < 2
%        set(gca,'XTickLabel','');
%        if i == 1
%          ylabel('Northing [m]','Fontsize',6,'interpreter','latex');
%          [LX,LY] = gca_latex(6,0,1);
%        else
%           set(gca,'YTickLabel',''); 
%        end  
%     else
if i == L
    xlabel('Easting [m]') %,'FontSize',6,'interpreter','latex')
%     rectangle('Position', [406500 844800+3.5d6 409000-406500 845500-844800]...
%             ,'LineStyle','--','EdgeColor','black')
    rectangle('Position',[407100 844700+3.5d6 2100 930],'LineStyle','--')
    colormap parula
    t = colorbar('SouthOutside');%^2H
    set(get(t,'xlabel'),'string','$\eta_{max}$ [m]','Interpreter','Latex')%Max. Tsunami Height
else
    set(gca,'XTickLabel',''); 
end
%        if i == 2
          ylabel('Northing [m]') %,'Fontsize',6,'interpreter','latex');
%          [LX,LY] = gca_latex(6,1,1);
%        else
%           set(gca,'YTickLabel',''); 
%           [LX,LY] = gca_latex(6,1,0);
%        end  
%         t = colorbar; %('SouthOutside');% Inserts a colorbar. a handle is created
%         set(gca,'FontSize',6)
%         colorbar_latex({'0','5','10','15','20'},6,'Y');
        %colorbar_latex({'0','0.05','0.10','0.15','0.20','0.25'},6,'Y');
        %set(get(t,'ylabel'),'string','$|\omega_{max}|$ [s$^{-1}$]',...
        %    'Fontsize',6,'Interpreter','Latex')
        %            set(get(t,'ylabel'),'string','$|u_{s,max}|$ [ms$^{-1}$]',...
        %        'Fontsize',6,'Interpreter','Latex')
    %     if Var == 1
    %         set(get(t,'ylabel'),'string','$\eta_{max}$ [m]',...
    %             'Fontsize',6,'Interpreter','Latex')
    %     else
    %         set(get(t,'xlabel'),'string','$|U_{max}|$ [ms$^{-1}$]',...
    %             'Fontsize',6,'Interpreter','Latex')
    %     end
%    end   
    %if i == L
      
%        if k == 1
%            colorbar_latex({'0','2','4','6','8','10'},6,'X');
%             set(get(t,'xlabel'),'string','$|U|$ [ms$^{-1}$]',...
%                 'Fontsize',6,'Interpreter','Latex')
%        else
%            colorbar_latex({'-0.2','0.1','0','0.1','0.2'},6,'Y');
%             set(get(t,'ylabel'),'string','$\omega$ [s$^{-1}$]',...
%                 'Fontsize',6,'Interpreter','Latex')
%        end
    %end
end
end
