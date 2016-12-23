% Compare fluxes
clearvars; clc; close all

direc{1} = 'Satakev8.0_2DH_BW_wall/DAT/';
direc{2} = 'Satakev8.0_2DH_BW_FD/DAT/';
direc{3} = 'Satakev8.0_3D_Real_BW/DAT/';
color = {'r--','b-.','k-'}; nn = 0;
figure;
hold on
for d = direc
    t = dlmread([d{1} 'XSECM_MIDX.dat'],'',2,0); 
    time = t(2:end,1);
    x    = t(1,2:end);
    flux = t(2:end,2:end);
    nn = nn + 1;
%     % Plot average per unit width
%     [~,~] = jbfill(time/60,max(flux,[],2),min(flux,[],2),...
%                    color{nn},color{nn},0,0.5);
    plot(time/60,-mean(flux,2),color{nn}); 
end
xlim([0 120])
ylim([-200 300])
xlabel('time since earthquake [min]')
ylabel('Q_x [m^2s^{-1}]')
legend('2DH, F_D = 0','2DH, F_D = 0.5','2CLOWNS-3D',...
       'location','NorthEast')
legend boxoff
set(gca,'box','on')
set(gca,'fontsize',7)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 9 7],'PaperPositionMode','manual');
print('-r600','-depsc','../Paper/Flux_Compare.eps'); %'-fillpage',

%legend(direc);