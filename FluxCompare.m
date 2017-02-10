% Compare fluxes
clearvars; clc; close all

direc{1} = '../Satakev8.0_2DH_BW_wall/DAT/';
direc{2} = '../Satakev8.0_2DH_BW_FD/DAT/';
direc{3} = '../Satakev8.0_3D_Real_BW/DAT/';
color = {'r--','b-.','k-'}; nn = 0;
figure;
hold on
for d = direc
    t = dlmread([d{1} 'XSECM_MIDX.dat'],'',2,0); 
    time = t(2:end,1);
    x    = t(1,2:end);
    flux = t(2:end,2:end);
    nn = nn + 1;
    
   if strcmp(d{1}(15),'3') 
       ym = x;
       load([d{1} '../MAT_files/Satakev8.0_3D_Real_BW_xyzF.mat'])
       j = knnsearch(y',ym'-5);
       nnm = mn(j,1018);
       
       fPath = cd;
       fname = '../Satake_BW_Real_ke2D.';
       List = dir(fullfile(fPath, [d{1} fname '*']));
       List([List.isdir]) = [];
       L = length(List);
       % order list by datenum
       Listc = struct2cell(List);
       Ordered = Listc'; %sortrows(Listc',5);
       
       for t = 1:L

            File = Ordered{t}; FL = length(File);
            tstrng = File(FL-7:FL-2);
            if strcmp(tstrng(1),'t')
                tstrng(1) = [];
            end
            timen = str2double(tstrng);
            if timen > time(end)
                time(end+1) = timen;

                load([d{1} '../MAT_files/' File '.mat'])
                %load([d{1} '../MAT_files/Satake_BW_Real_ke2D.t7200.9_5.mat'])
                flux_t = 0.5*(eta(nnm) - Zk(nnm) + eta(nnm-1) - Zk(nnm-1) ).*u(nnm,1);
                flux(end+1,:) = flux_t;
            end
        end
        plot(time/60,-mean(flux,2),color{nn}); 
   else
%     % Plot average per unit width
%     [~,~] = jbfill(time/60,max(flux,[],2),min(flux,[],2),...
%                    color{nn},color{nn},0,0.5);
       plot(time/60,-mean(flux,2),color{nn}); 
   end
end
xlim([0 120])
ylim([-200 300])
xlabel('time since earthquake [min]')
ylabel('\itQ_x\rm [m^2s^{-1}]')
legend('2DH, \itF_D\rm = 0','2DH, \itF_D\rm = 0.5','2CLOWNS-3D',...
       'location','NorthEast')
legend boxoff
set(gca,'box','on')
set(gca,'fontsize',7)
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 9 7],'PaperPositionMode','manual');
print('-r600','-depsc','../../Paper/Flux_Compare.eps'); %'-fillpage',

%legend(direc);