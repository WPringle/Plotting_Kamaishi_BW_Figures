clearvars; clc;
%% ReadBinary output
% direc = 'Satakev8.0_Real_ke/';
% fname = 'SatakeT_HYB_Real_ke2D.t';
% oname = 'Satakev8.0_Real_';

% direc = 'Satakev8.0_3D_Nobw/';
% fname = 'Satake_NOBW_Real_ke2D.t';
% oname = 'Satakev8.0_3D_NoBW_'; 

direc = 'Satakev8.0_3D_Real_BW/';
fname = 'Satake_BW_Real_ke.t';
oname = 'Satakev8.0_3D_Real_BW_'; 

% direc = 'Satakev8.0_2DH_Nobw/';
% fname = 'Satakev8.0_2DH_nobw2D.t';
% oname = 'Satakev8.0_2DH_';

% direc = 'Satakev8.0_2DH_BW/';
% fname = 'Satakev8.0_2DH_BW2D.t';
% oname = 'Satakev8.0_2DH_BW_';

% direc = 'Satakev8.0_2DH/';
% fname = 'Tohoku_2DH_Satake_T12D.t';
% oname = 'Satakev8.0_2DH_'; %'TUv1.2_source'; 

fname_i = '../Inputs/Tohoku_2D_';
layerc = {''}; %{'1','2','3','4','5'};
type = 'series';
for l =  1:length(layerc)
    layer = layerc{l};
    fPath = cd;
    List = dir(fullfile(fPath, [direc fname '*' layer ]));
    List([List.isdir]) = [];
    L = length(List);
    % order list by datenum
    Listc = struct2cell(List);
    if strcmp(layer,'')
        L = L/4;
        Ordered = Listc';
    else
        Ordered = sortrows(Listc',1);
    end
    for t = 1:L
        % Get t
        if strcmp(layer,'')
            Ordern = Ordered(t*4-3:t*4,:);
            Ordern = sortrows(Ordern,-3);
            % Get time
            File = Ordern{1}; FL = length(File);
            tstrng = File(FL-5:FL);
        else
            File = Ordered{t}; FL = length(File);
            tstrng = File(FL-7:FL-2);
        end
        if strcmp(tstrng(1),'t')
            tstrng(1) = [];
        end
        if strcmp(type,'series')
            fid = fopen([direc 'Series_files/' File '.nc.series'],'w');
            if strcmp(layer,'1')
                fprintf(fid,'%s',['1.4 5.1 12 0 0 LT ' tstrng '[s]']);
            elseif strcmp(layer,'5')
                fprintf(fid,'%s',['0.4041 4.3487 12 0 0 LT ' tstrng ' [s]']);
            elseif strcmp(layer,'')
                fprintf(fid,'%s',['0.40684 4.345725 12 0 0 LT ' tstrng ' [s]']);
            end
            fclose(fid);
        elseif strcmp(type,'bar')
            fid = fopen([direc 'Bar_files/' File '.nc.bar'],'w');
            width = 2*str2double(tstrng)/7200;
            if strcmp(layer,'1')
                start = 1.3686 + 0.5*width/7.6165;
                fprintf(fid,'%s',[num2str(start) ' 4.3485 ' num2str(width) ' 0.5']);
            elseif strcmp(layer,'5')
                start = 0.4021 + 0.5*width/1115;
                fprintf(fid,'%s',[num2str(start) ' 4.3485 ' num2str(width) ' 0.5']);
            elseif strcmp(layer,'')
                start = 0.406075 + 0.5*width/3097;
                fprintf(fid,'%s',[num2str(start,8) ' 4.34569 ' ...
                                  num2str(width,8) ' 0.5']);
            end
            fclose(fid);
        end
    end
end