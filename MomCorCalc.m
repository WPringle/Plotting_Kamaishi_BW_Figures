% The script calculates the correction coefficient on momentum at the
% breakwater section
clearvars;
clc;
close all;

%% Open the 3D files
direc = 'Satakev8.0_3D_Real_BW/';
load([direc 'MAT_files/Satake_BW_Real_ke.xyznF.mat'])

is = find ( x == 408150 );
js = find ( y >= 845050 & y <= 845280);
zs = find ( z > -22.5); zs(end) = [];
load([direc 'MAT_files/Satake_BW_Real_ke.t1908.3.mat'])

time = 1908;
dz = 2.2;

%% Calc Bxx
Bxx = zeros(length(time),length(js));
for t = 1:length(time)
    for j = 1:length(js)
        % find the i values where nf = -1 and get either side
        mn_w = mn(js(j),zs,is); mn_wm = mn(js(j),zs,is-1); 
        I = find( F(mn_w) <= 0 ); J = find( F(mn_wm) <= 0 );
        mn_w(max(length(I),length(J))) = []; 
        mn_wm(max(length(I),length(J))) = [];
        % Calc depth
        D = 0.5 * dz * sum(F(mn_w) + F(mn_wm)); 
        % Calc Mom
        M = 0.5 * dz * sum(u(mn_w,1) .* ( F(mn_w) + F(mn_wm) ) );
        % Calc Depth-averaged U
        U = M/D;
        % Calc integral(u'^2.dz,eta,bot)
        Udiff = 0.5 * dz * sum( (u(mn_w,1) - U).^2 .* ( F(mn_w) + F(mn_wm) ) );
        % Calc Bxx
        Bxx(t,j)   = 1 + D/M^2*Udiff;
    end
end
% Get the drag force
Bxxmin = min(Bxx,[],2); %kN/m
Bxxmax = max(Bxx,[],2); %kN/m

% Plot Bxx
figure(1);
[~,~] = jbfill(time/60,Bxxmax,Bxxmin,'b','b',0,0.5);
xlim([0 120])
xlabel('time since earthquake [min]')
ylabel('Bxx')