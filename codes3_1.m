clear
clc
% tip: 卫星是柱浓度，模式分层浓度
load('coastlines.mat');  % 加载海岸线数据
plot(coastlon,coastlat)
hold on

info_model = ncinfo('mod\GEOSChem.SpeciesConc.20140101_0000z.nc4');
info_sat = ncinfo('sat\20140101-C3S-L2_GHG-GHG_PRODUCTS-TANSO-GOSAT-OCFP-DAILY-v7.2.nc');

%% ------------------------------------1月份比较-----------------------------------
files_sat = dir('sat\*.nc');
files_model = dir('mod\*.nc4');
% 模式数据
ch4_model = ncread(['mod\',files_model(1).name],'SpeciesConc_CH4'); %  {'mol mol-1 dry'}
lon_model = ncread(['mod\',files_model(1).name],'lon');
lat_model = ncread(['mod\',files_model(1).name],'lat');
lev_model = ncread(['mod\',files_model(1).name],'lev'); % 模式数据气压层
time_model = ncread(['mod\',files_model(1).name],'time');
temp_model = [];
temp_sat = [];
for i = 1 : 31   % 1月31个文件
    disp(['正在处理1月',files_sat(i).name(7:8),'日数据'])
    lon_sat = ncread(['sat\',files_sat(i).name],'longitude');
    lat_sat = ncread(['sat\',files_sat(i).name],'latitude');
    lev_sat = ncread(['sat\',files_sat(i).name],'pressure_levels'); % sat气压层
    pw_sat = ncread(['sat\',files_sat(i).name],'pressure_weight');
    ch4_pa_sat = ncread(['sat\',files_sat(i).name],'ch4_profile_apriori');
    xch4_ak_sat = ncread(['sat\',files_sat(i).name],'xch4_averaging_kernel');
    xch4_sat = ncread(['sat\',files_sat(i).name],'xch4');  % 单位：ppb
    scatter(lon_sat,lat_sat)
    [m_lon,m_lat] = meshgrid(lon_model,lat_model);
    % 从高度插值后的模式数据选取出卫星观测位置的数据
    for n = 1 : length(lon_sat)
        for yy = 1 : size(m_lon,2)
            if lon_sat(n,1) < lon_model(yy,1)
                position_mark(n,1) = yy;  % 记录经度所在格点位置
                break
            else
                position_mark(n,1) = 144;
            end
        end
    end
    for n = 1 : length(lon_sat)
        for xx = 1 : size(m_lon,1)
            if lat_sat(n,1) < lat_model(xx,1)
                position_mark(n,2) = xx;  % 记录纬度所在格点位置
                break
            end
        end
    end
    % 转换成柱浓度
    for lon_x = 1 : length(lon_model)
        for lat_y = 1 : length(lat_model)
            ch4_interp_model(lon_x,lat_y,:,i) = ...
                interp1(lev_model*1000,squeeze(ch4_model(lon_x,lat_y,:,i)),lev_sat(:,1),'method'); % 将模式各层数据插值到同卫星同层
        end
    end
    for j = 1 : length(xch4_sat) % 每天观测位置数
        for k = 1 : 20
            ch4_slip_model(k,j) = ...
                (ch4_pa_sat(k,j) + (ch4_interp_model(position_mark(j,1),position_mark(j,2),k,i)*10^9 - ch4_pa_sat(k,j)) * xch4_ak_sat(k,j)) * pw_sat(k,j); % 转为柱浓度
        end
    end
    eval(['ch4_column_model',files_sat(i).name(5:8),'= sum(ch4_slip_model,1);']); % 每天各观测点的模式数据柱浓度
    ch4_column_model = sum(ch4_slip_model,1);
    % 卫星--模式数据比较
    RMSe1(i,1) = sqrt(sum((ch4_column_model - xch4_sat').^2)/length(xch4_sat));  %RMSE评估两者数据的差异
    
    temp_model = [temp_model;ch4_column_model'];
    temp_sat = [temp_sat;xch4_sat];
    
    clear position_mark ch4_slip_model;
end
ave_model1 = mean(temp_model);
ave_sat1 = mean(temp_sat);
%% ----------------------------------------2月份比较-----------------------------------
files_sat = dir('sat2\*.nc');
files_model = dir('mod\*.nc4');
% mod
ch4_model = ncread(['mod\',files_model(2).name],'SpeciesConc_CH4'); %  {'mol mol-1 dry'}
lon_model = ncread(['mod\',files_model(2).name],'lon');
lat_model = ncread(['mod\',files_model(2).name],'lat');
lev_model = ncread(['mod\',files_model(2).name],'lev'); % 模式数据气压层
time_model = ncread(['mod\',files_model(2).name],'time');
temp_model = [];
temp_sat = [];
for i = 1 : 26   % 1月31个文件
    disp(['正在处理2月',files_sat(i).name(7:8),'日数据'])
    lon_sat = ncread(['sat2\',files_sat(i).name],'longitude');
    lat_sat = ncread(['sat2\',files_sat(i).name],'latitude');
    lev_sat = ncread(['sat2\',files_sat(i).name],'pressure_levels'); % 卫星数据气压层
    pw_sat = ncread(['sat2\',files_sat(i).name],'pressure_weight');
    ch4_pa_sat = ncread(['sat2\',files_sat(i).name],'ch4_profile_apriori');
    xch4_ak_sat = ncread(['sat2\',files_sat(i).name],'xch4_averaging_kernel');
    xch4_sat = ncread(['sat2\',files_sat(i).name],'xch4');  % 单位：ppb
    %     scatter(lon_sat,lat_sat)
    [m_lon,m_lat] = meshgrid(lon_model,lat_model);
    % 从高度插值后的模式数据选取出卫星观测位置的数据
    for n = 1 : length(lon_sat)
        for yy = 1 : size(m_lon,2)
            if lon_sat(n,1) < lon_model(yy,1)
                position_mark(n,1) = yy;  % 记录经度所在格点位置
                break
            else
                position_mark(n,1) = 144;
            end
        end
    end
    for n = 1 : length(lon_sat)
        for xx = 1 : size(m_lon,1)
            if lat_sat(n,1) < lat_model(xx,1)
                position_mark(n,2) = xx;  % 记录纬度所在格点位置
                break
            end
        end
    end
    % 转换成柱浓度
    for lon_x = 1 : length(lon_model)
        for lat_y = 1 : length(lat_model)
            ch4_interp_model(lon_x,lat_y,:,i) = ...
                interp1(lev_model*1000,squeeze(ch4_model(lon_x,lat_y,:,i)),lev_sat(:,1),'method'); % 将模式各层数据插值到同卫星同层
        end
    end
    for j = 1 : length(xch4_sat) % 每天观测位置数
        for k = 1 : 20
            ch4_slip_model(k,j) = ...
                (ch4_pa_sat(k,j) + (ch4_interp_model(position_mark(j,1),position_mark(j,2),k,i)*10^9 - ch4_pa_sat(k,j)) * xch4_ak_sat(k,j)) * pw_sat(k,j); % 转为柱浓度
        end
    end
    eval(['ch4_column_model',files_sat(i).name(5:8),'= sum(ch4_slip_model,1);']); % 每天各观测点的mod柱浓度
    ch4_column_model = sum(ch4_slip_model,1);
    % 卫星--mod比较
    RMSe2(i,1) = sqrt(sum((ch4_column_model - xch4_sat').^2)/length(xch4_sat));  %RMSE评估两者数据的差异
    temp_model = [temp_model;ch4_column_model'];
    temp_sat = [temp_sat;xch4_sat];
    clear position_mark ch4_slip_model;
end
ave_model2 = mean(temp_model,'omitnan');
ave_sat2 = mean(temp_sat);

%% ----------------------------------------3月份比较----------------------------------
files_sat = dir('sat3\*.nc');
files_model = dir('mod\*.nc4');
% 模式数据
ch4_model = ncread(['mod\',files_model(3).name],'SpeciesConc_CH4'); %  {'mol mol-1 dry'}
lon_model = ncread(['mod\',files_model(3).name],'lon');
lat_model = ncread(['mod\',files_model(3).name],'lat');
lev_model = ncread(['mod\',files_model(3).name],'lev'); % 模式数据气压层
time_model = ncread(['mod\',files_model(3).name],'time');
temp_model = [];
temp_sat = [];
for i = 1 : 31   % 1月31个文件
    disp(['正在处理3月',files_sat(i).name(7:8),'日数据'])
    lon_sat = ncread(['sat3\',files_sat(i).name],'longitude');
    lat_sat = ncread(['sat3\',files_sat(i).name],'latitude');
    lev_sat = ncread(['sat3\',files_sat(i).name],'pressure_levels'); % 卫星数据气压层
    pw_sat = ncread(['sat3\',files_sat(i).name],'pressure_weight');
    ch4_pa_sat = ncread(['sat3\',files_sat(i).name],'ch4_profile_apriori');
    xch4_ak_sat = ncread(['sat3\',files_sat(i).name],'xch4_averaging_kernel');
    xch4_sat = ncread(['sat3\',files_sat(i).name],'xch4');  % 单位：ppb
    %     scatter(lon_sat,lat_sat)
    [m_lon,m_lat] = meshgrid(lon_model,lat_model);
    % 从高度插值后的模式数据选取出卫星观测位置的数据
    for n = 1 : length(lon_sat)
        for yy = 1 : size(m_lon,2)
            if lon_sat(n,1) < lon_model(yy,1)
                position_mark(n,1) = yy;  % 记录经度所在格点位置
                break
            else
                position_mark(n,1) = 144;
            end
        end
    end
    for n = 1 : length(lon_sat)
        for xx = 1 : size(m_lon,1)
            if lat_sat(n,1) < lat_model(xx,1)
                position_mark(n,2) = xx;  % 记录纬度所在格点位置
                break
            end
        end
    end
    % 转换成柱浓度
    for lon_x = 1 : length(lon_model)
        for lat_y = 1 : length(lat_model)
            ch4_interp_model(lon_x,lat_y,:,i) = ...
                interp1(lev_model*1000,squeeze(ch4_model(lon_x,lat_y,:,i)),lev_sat(:,1),'method'); % 将模式各层数据插值到同卫星同层
        end
    end
    for j = 1 : length(xch4_sat) % 每天观测位置数
        for k = 1 : 20
            ch4_slip_model(k,j) = ...
                (ch4_pa_sat(k,j) + (ch4_interp_model(position_mark(j,1),position_mark(j,2),k,i)*10^9 - ch4_pa_sat(k,j)) * xch4_ak_sat(k,j)) * pw_sat(k,j); % 转为柱浓度
        end
    end
    eval(['ch4_column_model',files_sat(i).name(5:8),'= sum(ch4_slip_model,1);']); % 每天各观测点的模式数据柱浓度
    ch4_column_model = sum(ch4_slip_model,1);
    % 卫星--模式数据比较
    RMSe3(i,1) = sqrt(sum((ch4_column_model - xch4_sat').^2)/length(xch4_sat));  %RMSE评估两者数据的差异
    temp_model = [temp_model;ch4_column_model'];
    temp_sat = [temp_sat;xch4_sat];
    clear position_mark ch4_slip_model;
end
ave_model3 = mean(temp_model);
ave_sat3 = mean(temp_sat);
%% ----------------------------------------4月份比较----------------------------------
files_sat = dir('sat4\*.nc');
files_model = dir('mod\*.nc4');
% 模式数据
ch4_model = ncread(['mod\',files_model(3).name],'SpeciesConc_CH4'); %  {'mol mol-1 dry'}
lon_model = ncread(['mod\',files_model(3).name],'lon');
lat_model = ncread(['mod\',files_model(3).name],'lat');
lev_model = ncread(['mod\',files_model(3).name],'lev'); % 模式数据气压层
time_model = ncread(['mod\',files_model(3).name],'time');
temp_model = [];
temp_sat = [];
for i = 1 : 30   % 1月31个文件
    disp(['正在处理4月',files_sat(i).name(7:8),'日数据'])
    lon_sat = ncread(['sat4\',files_sat(i).name],'longitude');
    lat_sat = ncread(['sat4\',files_sat(i).name],'latitude');
    lev_sat = ncread(['sat4\',files_sat(i).name],'pressure_levels'); % sat气压层
    pw_sat = ncread(['sat4\',files_sat(i).name],'pressure_weight');
    ch4_pa_sat = ncread(['sat4\',files_sat(i).name],'ch4_profile_apriori');
    xch4_ak_sat = ncread(['sat4\',files_sat(i).name],'xch4_averaging_kernel');
    xch4_sat = ncread(['sat4\',files_sat(i).name],'xch4');  % 单位：ppb
    %     scatter(lon_sat,lat_sat)
    [m_lon,m_lat] = meshgrid(lon_model,lat_model);
    % 从高度插值后的模式数据选取出卫星观测位置的数据
    for n = 1 : length(lon_sat)
        for yy = 1 : size(m_lon,2)
            if lon_sat(n,1) < lon_model(yy,1)
                position_mark(n,1) = yy;  % 记录经度所在格点位置
                break
            else
                position_mark(n,1) = 144;
            end
        end
    end
    for n = 1 : length(lon_sat)
        for xx = 1 : size(m_lon,1)
            if lat_sat(n,1) < lat_model(xx,1)
                position_mark(n,2) = xx;  % 记录纬度所在格点位置
                break
            end
        end
    end
    % 转换成柱浓度
    for lon_x = 1 : length(lon_model)
        for lat_y = 1 : length(lat_model)
            ch4_interp_model(lon_x,lat_y,:,i) = ...
                interp1(lev_model*1000,squeeze(ch4_model(lon_x,lat_y,:,i)),lev_sat(:,1),'method'); % 将模式各层数据插值到同卫星同层
        end
    end
    for j = 1 : length(xch4_sat) % 每天观测位置数
        for k = 1 : 20
            ch4_slip_model(k,j) = ...
                (ch4_pa_sat(k,j) + (ch4_interp_model(position_mark(j,1),position_mark(j,2),k,i)*10^9 - ch4_pa_sat(k,j)) * xch4_ak_sat(k,j)) * pw_sat(k,j); % 转为柱浓度
        end
    end
    eval(['ch4_column_model',files_sat(i).name(5:8),'= sum(ch4_slip_model,1);']); % 每天各观测点的模式数据柱浓度
    ch4_column_model = sum(ch4_slip_model,1);
    % 卫星--模式数据比较
    RMSe4(i,1) = sqrt(sum((ch4_column_model - xch4_sat').^2)/length(xch4_sat));  %RMSE评估两者数据的差异
    temp_model = [temp_model;ch4_column_model'];
    temp_sat = [temp_sat;xch4_sat];
    clear position_mark ch4_slip_model;
end
ave_model4 = mean(temp_model);
ave_sat4 = mean(temp_sat);
%% ----------------------------------------5月份比较----------------------------------
files_sat = dir('sat5\*.nc');
files_model = dir('mod\*.nc4');
% 模式数据
ch4_model = ncread(['mod\',files_model(5).name],'SpeciesConc_CH4'); %  {'mol mol-1 dry'}
lon_model = ncread(['mod\',files_model(5).name],'lon');
lat_model = ncread(['mod\',files_model(5).name],'lat');
lev_model = ncread(['mod\',files_model(5).name],'lev'); % 模式数据气压层
time_model = ncread(['mod\',files_model(5).name],'time');
temp_model = [];
temp_sat = [];
for i = 1 : 25   % 1月31个文件
    
    disp(['正在处理5月',files_sat(i).name(7:8),'日数据'])
    lon_sat = ncread(['sat5\',files_sat(i).name],'longitude');
    lat_sat = ncread(['sat5\',files_sat(i).name],'latitude');
    lev_sat = ncread(['sat5\',files_sat(i).name],'pressure_levels'); % sat气压层
    pw_sat = ncread(['sat5\',files_sat(i).name],'pressure_weight');
    ch4_pa_sat = ncread(['sat5\',files_sat(i).name],'ch4_profile_apriori');
    xch4_ak_sat = ncread(['sat5\',files_sat(i).name],'xch4_averaging_kernel');
    xch4_sat = ncread(['sat5\',files_sat(i).name],'xch4');  % 单位：ppb
    %     scatter(lon_sat,lat_sat)
    [m_lon,m_lat] = meshgrid(lon_model,lat_model);
    % 从高度插值后的模式数据选取出卫星观测位置的数据
    for n = 1 : length(lon_sat)
        for yy = 1 : size(m_lon,2)
            if lon_sat(n,1) < lon_model(yy,1)
                position_mark(n,1) = yy;  % 记录经度所在格点位置
                break
            else
                position_mark(n,1) = 144;
            end
        end
    end
    for n = 1 : length(lon_sat)
        for xx = 1 : size(m_lon,1)
            if lat_sat(n,1) < lat_model(xx,1)
                position_mark(n,2) = xx;  % 记录纬度所在格点位置
                break
            end
        end
    end
    % 转换成柱浓度
    for lon_x = 1 : length(lon_model)
        for lat_y = 1 : length(lat_model)
            ch4_interp_model(lon_x,lat_y,:,i) = ...
                interp1(lev_model*1000,squeeze(ch4_model(lon_x,lat_y,:,i)),lev_sat(:,1),'method'); % 将模式各层数据插值到同卫星同层
        end
    end
    for j = 1 : length(xch4_sat) % 每天观测位置数
        for k = 1 : 20
            ch4_slip_model(k,j) = ...
                (ch4_pa_sat(k,j) + (ch4_interp_model(position_mark(j,1),position_mark(j,2),k,i)*10^9 - ch4_pa_sat(k,j)) * xch4_ak_sat(k,j)) * pw_sat(k,j); % 转为柱浓度
        end
    end
    eval(['ch4_column_model',files_sat(i).name(5:8),'= sum(ch4_slip_model,1);']); % 每天各观测点的模式数据柱浓度
    ch4_column_model = sum(ch4_slip_model,1);
    % 卫星--模式数据比较
    RMSe5(i,1) = sqrt(sum((ch4_column_model - xch4_sat').^2)/length(xch4_sat));  %RMSE评估两者数据的差异
    temp_model = [temp_model;ch4_column_model'];
    temp_sat = [temp_sat;xch4_sat];
    clear position_mark ch4_slip_model;
end
ave_model5 = mean(temp_model);
ave_sat5 = mean(temp_sat);
%% ----------------------------------------6月份比较----------------------------------
files_sat = dir('sat6\*.nc');
files_model = dir('mod\*.nc4');
% mod
ch4_model = ncread(['mod\',files_model(6).name],'SpeciesConc_CH4'); %  {'mol mol-1 dry'}
lon_model = ncread(['mod\',files_model(6).name],'lon');
lat_model = ncread(['mod\',files_model(6).name],'lat');
lev_model = ncread(['mod\',files_model(6).name],'lev'); % 模式数据气压层
time_model = ncread(['mod\',files_model(6).name],'time');
temp_model = [];
temp_sat = [];
for i = 1 : 30   % 1月31个文件
    disp(['正在处理6月',files_sat(i).name(7:8),'日数据'])
    lon_sat = ncread(['sat6\',files_sat(i).name],'longitude');
    lat_sat = ncread(['sat6\',files_sat(i).name],'latitude');
    lev_sat = ncread(['sat6\',files_sat(i).name],'pressure_levels'); % sat气压层
    pw_sat = ncread(['sat6\',files_sat(i).name],'pressure_weight');
    ch4_pa_sat = ncread(['sat6\',files_sat(i).name],'ch4_profile_apriori');
    xch4_ak_sat = ncread(['sat6\',files_sat(i).name],'xch4_averaging_kernel');
    xch4_sat = ncread(['sat6\',files_sat(i).name],'xch4');  % 单位：ppb
    %     scatter(lon_sat,lat_sat)
    [m_lon,m_lat] = meshgrid(lon_model,lat_model);
    % 从高度插值后的模式数据选取出卫星观测位置的数据
    for n = 1 : length(lon_sat)
        for yy = 1 : size(m_lon,2)
            if lon_sat(n,1) < lon_model(yy,1)
                position_mark(n,1) = yy;  % 记录经度所在格点位置
                break
            else
                position_mark(n,1) = 144;
            end
        end
    end
    for n = 1 : length(lon_sat)
        for xx = 1 : size(m_lon,1)
            if lat_sat(n,1) < lat_model(xx,1)
                position_mark(n,2) = xx;  % 记录纬度所在格点位置
                break
            end
        end
    end
    % 转换成柱浓度
    for lon_x = 1 : length(lon_model)
        for lat_y = 1 : length(lat_model)
            ch4_interp_model(lon_x,lat_y,:,i) = ...
                interp1(lev_model*1000,squeeze(ch4_model(lon_x,lat_y,:,i)),lev_sat(:,1),'method'); % 将模式各层数据插值到同卫星同层
        end
    end
    for j = 1 : length(xch4_sat) % 每天观测位置数
        for k = 1 : 20
            ch4_slip_model(k,j) = ...
                (ch4_pa_sat(k,j) + (ch4_interp_model(position_mark(j,1),position_mark(j,2),k,i)*10^9 - ch4_pa_sat(k,j)) * xch4_ak_sat(k,j)) * pw_sat(k,j); % 转为柱浓度
        end
    end
    eval(['ch4_column_model',files_sat(i).name(5:8),'= sum(ch4_slip_model,1);']); % 每天各观测点的模式数据柱浓度
    ch4_column_model = sum(ch4_slip_model,1);
    % 卫星--模式数据比较
    RMSe6(i,1) = sqrt(sum((ch4_column_model - xch4_sat').^2)/length(xch4_sat));  %RMSE评估两者数据的差异
    temp_model = [temp_model;ch4_column_model'];
    temp_sat = [temp_sat;xch4_sat];
    clear position_mark ch4_slip_model;
end
ave_model6 = mean(temp_model);
ave_sat6 = mean(temp_sat);
rmse = [mean(RMSe1);mean(RMSe2);mean(RMSe3);mean(RMSe4);mean(RMSe5);mean(RMSe6)];
disp('计算完毕！！！')
figure(2)
ave_model = [ave_model1,ave_model2,ave_model3,ave_model4,ave_model5,ave_model6];
ave_sat = [ave_sat1,ave_sat2,ave_sat3,ave_sat4,ave_sat5,ave_sat6];
plot(1:6,ave_model,'r-');
hold on
plot(1:6,ave_sat,'b-');

