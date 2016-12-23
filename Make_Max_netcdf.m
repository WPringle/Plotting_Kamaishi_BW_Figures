% Make the grd netcdf for GMT
clearvars; clc; close all;

dir_num = 'Satakev8.0_2DH_BW_FD\Satakev8.0_2DH_BW_FD2D.t7200.9MAX_5.dat';
%dir_num = 'Satakev8.0_3D_Nobw\Satake_NOBW_Real_ke2D.t7200.9MAX_5.dat';
%dir_num = 'Satakev8.0_3D_Real_BW\Satake_BW_Real_ke2D.t7200.9MAX_5.dat';
%dir_num = 'Satakev8.0_ke\Tohoku_SatakeT_HYB_ke_k1mm_wall2D.t7200.9MAX_5.dat';
%dir_num = 'Satakev8.0_Norm_Real_ke\SatakeT_Norm_HYB_realizable_ke2D.t7200.9MAX_5.dat';
%dir_num = 'Satakev8.0_NoRAN\SatakeT_Norm_HYB_NoRAN2D.t7200.9MAX_5.dat';

[ xx, yy, ETA_max ] = GetMaxData2015( dir_num, 3 );
xx = xx * 1d-6; yy = yy * 1d-6 + 3.5;

% Plot
figure;
pcolor(xx,yy,ETA_max); 
shading interp
caxis([0 20])

% M = csvread('../Inputs/ttjt_survey_Ryoishi_Kamaishi_Toni_tidecorrected.csv');
% 
% [B_x, B_y] = latlon_to_UTM(M(:,2), M(:,1), 0, 143);
% B_x = B_x + 500000;
% %B_y = B_y - 3.5d6;

%hold on
%scatter(B_x*1d-6,B_y*1d-6,[],M(:,3),'filled','MarkerEdgeColor','k')
% Write out to .nc
grdwrite2(xx,yy,ETA_max,'Satakev8.0_2DH_BW_FD_MaxVel.nc')
