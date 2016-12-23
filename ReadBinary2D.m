clearvars; clc;
%% ReadBinary output
% direc = 'Satakev8.0_Real_ke/';
% fname = 'SatakeT_HYB_Real_ke2D.t';
% oname = 'Satakev8.0_Real_';

% direc = 'Satakev8.0_3D_Real_BW/';
% fname = 'Satake_BW_Real_ke2D.t';
% oname = 'Satakev8.0_3D_Real_BW_'; 

% direc = 'Satakev8.0_3D_Nobw/';
% fname = 'Satake_NOBW_Real_ke2D.t';
% oname = 'Satakev8.0_3D_NoBW_'; 

% direc = 'Satakev8.0_2DH_Nobw/';
% fname = 'Satakev8.0_2DH_nobw2D.t';
% oname = 'Satakev8.0_2DH_';
% 
direc = 'Satakev8.0_2DH_BW_FD/';
fname = 'Satakev8.0_2DH_BW_FD2D.t';
oname = 'Satakev8.0_2DH_BW_FD';

% direc = 'Satakev8.0_2DH/';
% fname = 'Tohoku_2DH_Satake_T12D.t';
% oname = 'Satakev8.0_2DH_'; %'TUv1.2_source'; 

fname_i = '../Inputs/Tohoku_2D_';
layerc = {'5'}; %{'1','2','3','4','5'};
output = 'collect'; %bathy'; %'bathy'; %
for l =  1:length(layerc)
    layer = layerc{l};
%     if l == 5
%         print_topo = 1;
%     else
%         print_topo = 0;
%     end
    if strcmp(output,'bathy') 
        L = 1;
        maindata = [fname_i layer ];
        fdata = [fname_i layer '.fdata'];
    else
        fPath = cd;
        List = dir(fullfile(fPath, [direc fname '*' layer ]));
        List([List.isdir]) = [];
        L = length(List);
        % order list by datenum
        Listc = struct2cell(List);
        Ordered = sortrows(Listc',1);
        maindata = [direc fname '0.000_' layer];
        fdata = [direc fname '302.6_' layer '.fdata'];
    end
    if strcmp(output,'out') || strcmp(output,'bathy')
        %% Read xyz data
        [ inne, is, ie, x, y, z, nf, nfb , mn, in ] = ...
                                    read_xyzn2D( [fname_i layer '.xyzn'] );

        %% Read the fdata
        [ Hx, Man, Zk ] = read_fdata_2D( fdata , inne, z(4) ); 
        Zk  = z(4)-(z(4) - z(3))*Zk;
        Hx  = z(4)-(z(4) - z(3))*Hx;
        
        save([direc 'MAT_files/' oname 'xyzF.mat'],'inne','is','ie','x',...
                              'y','z','nf','nfb','mn','in','Hx','Man','Zk')
                                     
        % Make desired target into grid format
        xx = 0.5*(x(is(1):ie(1)-1) + x(is(1)+1:ie(1)));
        yy = 0.5*(y(is(2):ie(2)-1) + y(is(2)+1:ie(2))) + 3.5d6;
        xx = xx/1d6;
        yy = yy/1d6;
    else
        % Load data
        load([direc 'MAT_files/' oname 'xyzF.mat'])

        Vel = zeros(L,inne,2);
        Eta = zeros(L,inne);
    end

%     if strcmp(output,'out')
%         s     = length(xx)*length(yy);
%         zeta  = -9999*ones(s,L); 
%         qx    = zeros(s,1); 
%         qy    = zeros(s,1); 
%         bathy = NaN(s,1); 
%         time  = zeros(L,1);
%     elseif strcmp(output,'bathy')
%         bathy = NaN(length(yy),length(xx)); 
%     end
    %% Now loop over all files
    for t = 1:L
        t
        if strcmp(output,'out')
            % Get t
            File = Ordered{t}; FL = length(File);
            tstrng = File(FL-7:FL-2);
            if strcmp(tstrng(1),'t')
                tstrng(1) = [];
            end
            time(t) = str2double(tstrng);
            
            % Need to reset nf and nfb
            load([direc 'MAT_files/' oname 'xyzF.mat']);  
            
            % Read the main data
            [ eta, u, F, nf, nfb ] = read_main_data2D(  ...
                                        [direc File], inne, in, nf, nfb ); 

            save([direc 'MAT_files/' File '.mat'],'nf','nfb','F','u','eta');
        elseif strcmp(output,'collect')
            % Get t
            File = Ordered{t}; FL = length(File);
            tstrng = File(FL-7:FL-2);
            if strcmp(tstrng(1),'t')
                tstrng(1) = [];
            end
            time(t) = str2double(tstrng);
            
            % Load data
            load([direc 'MAT_files/' File '.mat']);
            
            %depth = max(0,eta-Zk);
            Vel(t,:,:) = u;
            Eta(t,:)   = eta;

%             % Make what we want
%             for j = is(2):ie(2)-1
%                 for i = is(1):ie(1)-1
%                     nn = mn(j,i);
%                     if nn ~=0
%                         if F(nn) == 0 ; continue; end
%                         nnpi = mn(j,i+1); 
%                         nnpj = mn(j+1,i);
%                         if nnpi ~=0
%                             u1 = 0.5 * (u(nn,1) + u(nnpi,1)); 
%                         else
%                             u1 = 0.5 * u(nn,1); 
%                         end
%                         if nnpj ~=0
%                             u2 = 0.5 * (u(nn,2) + u(nnpj,2));
%                         else
%                             u2 = 0.5 * u(nn,2); 
%                         end
%                         Depth = max(0,eta(nn) - Zk(nn));
%                         QX(t,j,i) = Depth * u1;
%                         QY(t,j,i) = Depth * u2;
%                         QU(t,j,i) = Depth * (u1^2+u2^2);
%                     end
%                 end
%             end

        end
%         if strcmp(output,'out')
%             Eta = NaN(length(yy),length(xx));
%             jj = 0;
%             for j = is(2):ie(2)-1
%                 jj = jj + 1; ii = 0;
%                 for i = is(1):ie(1)-1
%                     ii = ii + 1;
%                     nn = mn(j,i);
%                     nnip = mn(j,i+1);
%                     nnjp = mn(j+1,i);
%                     if nn ~=0
%                         % Make the array
%                         if strcmp(output,'out')
%                             if F(nn) > 0
%                                 %Eta(jj,ii) = eta(nn);
%                                 if nnip ~=0
%                                     umid = 0.5*(u(nn,1) + u(nnip,1));
%                                 else
%                                     umid = 0.5*u(nn,1);
%                                 end
%                                 if nnjp ~=0
%                                     vmid = 0.5*(u(nn,2) + u(nnjp,2));
%                                 else
%                                     vmid = 0.5*u(nn,2);
%                                 end
%                                 Eta(jj,ii) = sqrt(umid^2+vmid^2);
%                             end
%                         elseif strcmp(output,'bathy')
%                             Znow = (z(4) - z(3)) * ( 1 - Zk(nn)) + z(3);
%                             if print_topo == 0 && Znow > 0; continue; end
%                             bathy(jj,ii) = Znow;
%                         end
%                     end
%                 end
%             end
%             grdwrite2(xx,yy,Eta,[direc 'NETCDF_files/VEL/' File '.nc'])
%         end
%          if strcmp(output,'bathy')
%              grdwrite2(xx,yy,bathy,[oname layer '_' output '.nc'])
%          elseif strcmp(output,'out')
%              grdwrite2(xx,yy,eta,[direc oname layer '_' tstrng '.nc'])
%     %         I = find(isnan(bathy)); 
%     %         bathy(I) = []; qx(I) = []; qy(I) = []; zeta(I,:) = [];
%     %         tri = delaunay(qx,qy);
%     %         adcirc_zeta_write(time,qx,qy,tri,bathy,zeta,[oname layer '_' output '.nc'])
%     %     elseif strcmp(output,'u')
%     %         adcirc_vel_write(qx,qy,time,u,v,[oname layer '_' output '.nc'])
%          end
    end
end
if strcmp(output,'collect')
    save([direc 'MAT_files/' oname 'VelEta.mat'],'time','Vel','Eta','-v7.3');
end