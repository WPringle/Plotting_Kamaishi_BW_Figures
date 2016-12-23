%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      This M-File will read a DAT file that is written using            %
%        the "data_output" subroutine in the FORTRAN program             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%BY: WILLIAM PRINGLE 8 AUG 2014
%
clearvars; close all;
%Set Directory of Figure output
%Directory of DATA
dir_num = 'TUv1.2/'; 
%dir_num1 = 'Satakev8.0_s/';
dir_num2 = 'Satakev8.0_ke/DAT/';
dir_meas =  '../Inputs/NOWPHAS_Tsunami_data/';
line_color = {'-k','-r','-b','-g','-y'};
legend_name = {'Central Miyagi','South Iwate','North Miyagi','Cental Iwate'};
L = 4;
hlines = 2; %specify the number of lines of headers
hlinesmeas = 3;
%Directory of MEASUREMENT or second DATA
%Input variables here%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
h=0;              % initial water depth - we use this to nondimensionlise
                  % the variables (set zero for dimensionless values)
t=0;              % time for output (cross-section data)
                  % (maybe a vector to create movies). Set to zero when
                  % plotting wave guage data. Set to dimensionless if
                  % using dimensionless values and dimensional otherwise
t1=0;             % Initial time value (adjust between simulation 
                  % and experiment etc). This is always dimensional 
var = 1;          % Variable time for wave guage data:
                  % var = 1 : Free-surface
                  %       2 : Depth
                  %       3 : U velocity
                  %       4 : V velocity
                  %       5 : Specific momentum flux
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.05], [0.11 0.01], [0.08 0.018]);
figure % create new figure
for i = 1:L
    num = dlmread([dir_num '80' num2str(i) '.dat'],'',hlines,0);
    %num1 = dlmread([dir_num1 '80' num2str(i) '.dat'],'',hlines,0);
    num2 = dlmread([dir_num2 '80' num2str(i) '.dat'],'',hlines,0);
    meas = dlmread([dir_meas '2011TET80' num2str(i) 'G.txt'],',',hlinesmeas,1);
    start = dlmread([dir_meas '2011TET80' num2str(i) 'G.txt'],'',[1 2 1 2]);
    [ l, ~ ] = size(meas);
    %Get the desired arrays from the full matrix data
    [ time, eta] = GetArraysWG( num, var, t1, h);
    %[ time1, eta1] = GetArraysWG( num1, var, t1, h);
    [ time2, eta2] = GetArraysWG( num2, var, t1, h);
    [ xm, etam ] = GetArraysWG( meas, var+2, t1, h );
    
    I = find(xm == 0); 
    if start == 11
        etam(I(2):l) = []; xm(I(2):l) = [];
    elseif start == 10
        etam(I(3):l) = []; xm(I(3):l) = [];
        etam(1:I(2)-1) = []; xm(1:I(2)-1) = [];
    end
    I = find(etam == 9999.99); etam(I) = []; xm(I) = [];
    etam = etam * 0.01;
    ym = num2str(xm);
    timem = zeros(length(xm),1);
    for n=1:length(ym)
        for m=5:-2:1
            if ~isnan(str2double(ym(n,m:m+1)))
                if m == 5; p = 1;end
                if m == 3; p = 60;end
                if m == 1; p = 3600;end
                timem(n) = str2double(ym(n,m:m+1)) * p + timem(n);
            end
        end
    end           %14:46:23  
    timem = timem - 53183; %subtract time of earthquake
    time  = time + 155; % Extra time for instantaneous assumption
    %time1  = time1 + 180; % Extra time for instantaneous assumption
    %% Plot the free surface time series
    subplot(L,1,i);
    plot(timem/60,etam,'-k')
    hold on
    %plot(time/60,eta,'-g')
    %plot(time1/60,eta1,'-b')
    plot(time2/60,eta2,'-.r')
    
    ylim([-10 10]);
    xlim([0 120]);
    ylabel('\eta [m]')
    if i ~= L
       set(gca,'XTickLabel','');
    else
       legend('Observed','Satake v8.0 Kinematic Simulated') %'Satake v8.0 Static',
       xlabel('time since earthquake rupture [min]')
    end
    text(80,6,legend_name{i})
    
end

