clearvars; clc; close all;
%% ReadBinary output
direc = 'Satakev8.0_3D_Real_BW/';
fname_o = 'Satake_BW_Real_ke.';
fname_i = '../Inputs/Satakev8.0_bw_3D_z2.2';
oname = 'Satakev8.0_3D_Real_BW_'; 
%time = 't1837.3';

output = 'collect';

fPath = cd;
List = dir(fullfile(fPath, [direc fname_o '*']));
List([List.isdir]) = [];
L = length(List)/4;
% order list by datenum
Listc = struct2cell(List);
Ordered = Listc'; %sortrows(Listc',5);
% 
if strcmp(output,'out')
    % Read xyz data
    [ inne, inne2, is, ie, x, y, z, nf, nfb, mn, in ] = ...
                                          read_xyzn3D( [fname_i '.xyzn'] );
    
    % Read fdata
    [ fb, a ] = read_fdata_3D( [fname_i '.fdata'] , inne, inne2);                            
    %                                      
    save([direc 'MAT_files/' fname_o 'xyznF.mat'],'inne', 'inne2', ...
                   'is', 'ie','x', 'y', 'z','mn','in','nf','nfb','fb','a');  
else
    % Load data
    load([direc 'MAT_files/' fname_o 'xyznF.mat'])
    xx = 0.5*(x(is(1):ie(1)-1)+x(is(1)+1:ie(1)));
    yy = 0.5*(y(is(2):ie(2)-1)+y(is(2)+1:ie(2)));
    xx = xx*1e-6;
    yy = yy*1e-6 + 3.5;
    kk = knnsearch(z',-10);
end
time  = zeros(L,1);
count = 0;
for t = 1:L
    t
    
    Ordern = Ordered(t*4-3:t*4,:);
    Ordern = sortrows(Ordern,-3);
    % Get time
    File = Ordern{1}; FL = length(File);
    tstrng = File(FL-5:FL);
    if strcmp(tstrng(1),'t')
        tstrng(1) = [];
    end
    timen = str2double(tstrng);
    %if timen < 1850 || timen > 2000 ; continue; end
    
    count = count + 1;
    
    time(count) = timen;
    
    if strcmp(output,'out')
        % Need to reset nf and nfb
        load([direc 'MAT_files/' fname_o 'xyznF.mat']);  

        %% Read the main data       
        [ p, u, F, nf, nfb ] = read_main_data3D( [direc File], ...
                                                  inne2, in, nf, nfb ); 
        % 
        %% Read the rdata
        [ rk, eps, nu_t, pro ] = read_rdata3D( [direc File '.rdata'], ...
                                               inne, in, nf );  

        save([direc 'MAT_files/' File '.mat'],'nf','nfb','F','u','p',...
                                             'rk','eps','nu_t','pro');
    elseif strcmp(output,'collect')       
       %
       load([direc 'MAT_files/' File '.mat'])
     
%        if ~exist('P3','var')
% %            Vel3 = zeros(L,size(u,1),size(u,2)); 
% %            F3   = zeros(L,length(F)); 
%            P3   = zeros(L,length(p)); 
% %            turb_vis = zeros(L,length(nu_t)); 
% %            turb_pro = zeros(L,length(pro)); 
% %            turb_eps = zeros(L,length(eps)); 
% %            turb_k = zeros(L,length(rk)); 
%        end
        
        Eta = NaN(length(yy),length(xx));
        jj = 0;
        for j = is(2):ie(2)-1
            jj = jj + 1; ii = 0;
            for i = is(1):ie(1)-1
                ii = ii + 1;
                nn = mn(j,kk,i);
                nnip = mn(j,kk,i+1);
                nnjp = mn(j+1,kk,i);
                if nn ~=0
                    % Make the array
                    if F(nn) > 0
                         Eta(jj,ii) = rk(nn);
%                         Eta(jj,ii) = nu_t(nn);
%                             %Eta(jj,ii) = eta(nn);
%                             if nnip ~=0
%                                 umid = 0.5*(u(nn,1) + u(nnip,1));
%                             else
%                                 umid = 0.5*u(nn,1);
%                             end
%                             if nnjp ~=0
%                                 vmid = 0.5*(u(nn,2) + u(nnjp,2));
%                             else
%                                 vmid = 0.5*u(nn,2);
%                             end
%                             Eta(jj,ii) = sqrt(umid^2+vmid^2);
                    end
                end
            end
        end
        grdwrite2(xx,yy,Eta,[direc 'NETCDF_files/TKE/' File '.nc'])
%        F3(count,:) = F;
%        Vel3(count,:,:) = u;
%       P3(count,:) = p;
%    turb_vis(count,:) = nu_t;
%    turb_pro(count,:) = pro;
%    turb_eps(count,:) = eps;
%    turb_k(count,:) = rk;
    end
end
% vel = vel(1:count,:,:);
% turb_vis = turb_vis(1:count,:);
% turb_pro = turb_pro(1:count,:);
% turb_eps = turb_eps(1:count,:);
% turb_k = turb_k(1:count,:);
% time = time(1:count);
%     
% load([direc 'MAT_files/' fname_o 'xyzn.mat']);  
% u    = squeeze(mean(vel,1));
% nu_t = squeeze(mean(turb_vis,1))';
% turb_pro = squeeze(mean(turb_pro,1))';
% turb_k = squeeze(mean(turb_k,1))';
% turb_eps = squeeze(mean(turb_eps,1))';
% 
% %save([direc fname_o 'pro.mat'],'time','turb_pro');
% 
%if strcmp(output,'collect')   
%    save([direc 'MAT_files/' fname_o 'P.mat'],'time','P3','-v7.3');
%end