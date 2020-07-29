clear
clc

path = '20200702\';
files = dir(path);

ae_info = ncinfo([path,files(3).name,'\ae20120522_20181031.public.nc']);
ae_xch4 = ncread([path,files(3).name,'\ae20120522_20181031.public.nc'],'xch4_ppm');
ae_time(:,1) = ncread([path,files(3).name,'\ae20120522_20181031.public.nc'],'year');
ae_time(:,2) = ncread([path,files(3).name,'\ae20120522_20181031.public.nc'],'day');
ae_lon = ncread([path,files(3).name,'\ae20120522_20181031.public.nc'],'long_deg');
ae_lat = ncread([path,files(3).name,'\ae20120522_20181031.public.nc'],'lat_deg');

[ae_lonmax,ae_latmax,ae_lonmin,ae_latmin] = getxy(ae_lon,ae_lat);

Bremen_info = ncinfo([path,files(4).name,'\br20100122_20190425.public.nc']);
br_time(:,1) = ncread([path,files(4).name,'\br20100122_20190425.public.nc'],'year');
br_time(:,2) = ncread([path,files(4).name,'\br20100122_20190425.public.nc'],'day');
br_lon = ncread([path,files(4).name,'\br20100122_20190425.public.nc'],'long_deg');
br_lat = ncread([path,files(4).name,'\br20100122_20190425.public.nc'],'lat_deg');
br_xch4 = ncread([path,files(4).name,'\br20100122_20190425.public.nc'],'xch4_ppm');

[br_lonmax,br_latmax,br_lonmin,br_latmin] = getxy(br_lon,br_lat);

darwin_info = ncinfo([path,files(8).name,'\db20050828_20190328.public.nc']);
db_time(:,1) = ncread([path,files(8).name,'\db20050828_20190328.public.nc'],'year');
db_time(:,2) = ncread([path,files(8).name,'\db20050828_20190328.public.nc'],'day');
db_lon = ncread([path,files(8).name,'\db20050828_20190328.public.nc'],'long_deg');
db_lat = ncread([path,files(8).name,'\db20050828_20190328.public.nc'],'lat_deg');
db_xch4 = ncread([path,files(8).name,'\db20050828_20190328.public.nc'],'xch4_ppm');

[db_lonmax,db_latmax,db_lonmin,db_latmin] = getxy(db_lon,db_lat);

Izana_info = ncinfo([path,files(5).name,'\iz20070518_20200128.public.nc']);
iz_time(:,1) = ncread([path,files(5).name,'\iz20070518_20200128.public.nc'],'year');
iz_time(:,2) = ncread([path,files(5).name,'\iz20070518_20200128.public.nc'],'day');
iz_lon = ncread([path,files(5).name,'\iz20070518_20200128.public.nc'],'long_deg');
iz_lat = ncread([path,files(5).name,'\iz20070518_20200128.public.nc'],'lat_deg');
iz_xch4 = ncread([path,files(5).name,'\iz20070518_20200128.public.nc'],'xch4_ppm');

[iz_lonmax,iz_latmax,iz_lonmin,iz_latmin] = getxy(iz_lon,iz_lat);

Lamont_info = ncinfo([path,files(6).name,'\oc20080706_20200101.public.nc']);
Lamont_time(:,1) = ncread([path,files(6).name,'\oc20080706_20200101.public.nc'],'year');
Lamont_time(:,2) = ncread([path,files(6).name,'\oc20080706_20200101.public.nc'],'day');
Lamont_lon = ncread([path,files(6).name,'\oc20080706_20200101.public.nc'],'long_deg');
Lamont_lat = ncread([path,files(6).name,'\oc20080706_20200101.public.nc'],'lat_deg');
Lamont_xch4 = ncread([path,files(6).name,'\oc20080706_20200101.public.nc'],'xch4_ppm');

[Lamont_lonmax,Lamont_latmax,Lamont_lonmin,Lamont_latmin] = getxy(Lamont_lon,Lamont_lat);

Tsukuba_info = ncinfo([path,files(7).name,'\tk20110804_20190423.public.nc']);
tk_time(:,1) = ncread([path,files(7).name,'\tk20110804_20190423.public.nc'],'year');
tk_time(:,2) = ncread([path,files(7).name,'\tk20110804_20190423.public.nc'],'day');
tk_lon = ncread([path,files(7).name,'\tk20110804_20190423.public.nc'],'long_deg');
tk_lat = ncread([path,files(7).name,'\tk20110804_20190423.public.nc'],'lat_deg');
tk_xch4 = ncread([path,files(7).name,'\tk20110804_20190423.public.nc'],'xch4_ppm');

[tk_lonmax,tk_latmax,tk_lonmin,tk_latmin] = getxy(tk_lon,tk_lat);

wollongong_info = ncinfo([path,files(9).name,'\wg20080626_20190430.public.nc']);
wg_time(:,1) = ncread([path,files(9).name,'\wg20080626_20190430.public.nc'],'year');
wg_time(:,2) = ncread([path,files(9).name,'\wg20080626_20190430.public.nc'],'day');
wg_lon = ncread([path,files(9).name,'\wg20080626_20190430.public.nc'],'long_deg');
wg_lat = ncread([path,files(9).name,'\wg20080626_20190430.public.nc'],'lat_deg');
wg_xch4 = ncread([path,files(9).name,'\wg20080626_20190430.public.nc'],'xch4_ppm');

[wg_lonmax,wg_latmax,wg_lonmin,wg_latmin] = getxy(wg_lon,wg_lat);


%% 显示站点位置
load('coastlines.mat'); % 加载海岸线数据
plot(coastlon,coastlat)
hold on
lon = [ae_lon(1,1);br_lon(1,1);db_lon(1,1);iz_lon(1,1);Lamont_lon(1,1);tk_lon(1,1);wg_lon(1,1)];
lat = [ae_lat(1,1);br_lat(1,1);db_lat(1,1);iz_lat(1,1);Lamont_lat(1,1);tk_lat(1,1);wg_lat(1,1)];
scatter(lon,lat,'r') % 画出站点位置

%% 按月存取数据
[ae1,ae2,ae3,ae4,ae5,ae6] = obtime2mon(ae_time,ae_lon,ae_lat,ae_xch4);
[br1,br2,br3,br4,br5,br6] = obtime2mon(br_time,br_lon,br_lat,br_xch4);
[db1,db2,db3,db4,db5,db6] = obtime2mon(db_time,db_lon,db_lat,db_xch4);
[iz1,iz2,iz3,iz4,iz5,iz6] = obtime2mon(iz_time,iz_lon,iz_lat,iz_xch4);
[Lamont1,Lamont2,Lamont3,Lamont4,Lamont5,Lamont6] = obtime2mon(Lamont_time,Lamont_lon,Lamont_lat,Lamont_xch4);
[tk1,tk2,tk3,tk4,tk5,tk6] = obtime2mon(tk_time,tk_lon,tk_lat,tk_xch4);
[wg1,wg2,wg3,wg4,wg5,wg6] = obtime2mon(wg_time,wg_lon,wg_lat,wg_xch4);

%% 计算3月情况
files_sat = dir('sat3\*.nc');
files_model = dir('mod\*.nc4');
% 模式数据
ch4_model = ncread(['mod\',files_model(1).name],'SpeciesConc_CH4'); %  {'mol mol-1 dry'}
lon_model = ncread(['mod\',files_model(1).name],'lon');
lat_model = ncread(['mod\',files_model(1).name],'lat');
lev_model = ncread(['mod\',files_model(1).name],'lev'); % 模式数据气压层
time_model = ncread(['mod\',files_model(1).name],'time');

% 提取各站点一天多个值提取为一天一个值
[ae_station] = unique_stations(ae3); % 获取站点某日柱浓度数据
[br_station] = unique_stations(br3); % 获取站点某日柱浓度数据
[db_station] = unique_stations(db3);
[iz_station] = unique_stations(iz3);
[Lamont_station] = unique_stations(Lamont3);
[tk_station] = unique_stations(tk3);
[wg_station] = unique_stations(wg3);

% 结合卫星、模式数据开始计算
for i = 1 : 31   % 3月31个文件
    disp(['正在处理3月',files_sat(i).name(7:8),'日数据'])
    lon_sat = ncread(['sat3\',files_sat(i).name],'longitude');
    lat_sat = ncread(['sat3\',files_sat(i).name],'latitude');
    lev_sat = ncread(['sat3\',files_sat(i).name],'pressure_levels'); % sat气压层
    pw_sat = ncread(['sat3\',files_sat(i).name],'pressure_weight');
    ch4_pa_sat = ncread(['sat3\',files_sat(i).name],'ch4_profile_apriori');
    xch4_ak_sat = ncread(['sat3\',files_sat(i).name],'xch4_averaging_kernel');
    xch4_sat = ncread(['sat3\',files_sat(i).name],'xch4');  % 单位：ppb
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
    
    %----------卫星、模式、站点相关性-------------
    % 第一个站（ae）三种数据归纳整理
%     [modsatsta] = output3data(i,lon_sat,lat_sat,ae_lonmax,ae_latmax,ae_lonmin,ae_latmin,ae_station,ch4_column_model,xch4_sat);
    [p] = sectxy(lon_sat,lat_sat,ae_lonmax,ae_latmax,ae_lonmin,ae_latmin) % 获取每日落在ae站点4*4位置处的卫星观测经纬度所在行
    if  p > 0 
        for sta = 1 : size(ae_station,1)
            if i == ae_station(sta,2)-59
                ae_modsatsta(i,:) = [ch4_column_model(p),xch4_sat(p),ae_station(sta,5)]; % 按日期整理卫星和模式柱浓度
            end
        end
    end
     % 第二个站（br）三种数据归纳整理
    [p] = sectxy(lon_sat,lat_sat,br_lonmax,br_latmax,br_lonmin,br_latmin) % 获取每日落在br站点4*4位置处的卫星观测经纬度所在行
    if  p > 0 
        for sta = 1 : size(br_station,1)
            if i == br_station(sta,2)-59
                br_modsatsta(i,:) = [ch4_column_model(p),xch4_sat(p),br_station(sta,5)]; % 按日期整理卫星和模式柱浓度
            end
        end
    end
    % 第3个站（db）三种数据归纳整理
    [p] = sectxy(lon_sat,lat_sat,db_lonmax,db_latmax,db_lonmin,db_latmin) % 获取每日落在br站点4*4位置处的卫星观测经纬度所在行
    if  p > 0 
        for sta = 1 : size(db_station,1)
            if i == db_station(sta,2)-59
                db_modsatsta(i,:) = [ch4_column_model(p),xch4_sat(p),db_station(sta,5)]; % 按日期整理卫星和模式柱浓度
            end
        end
    end
    % 第4个站（iz）三种数据归纳整理
    [p] = sectxy(lon_sat,lat_sat,iz_lonmax,iz_latmax,iz_lonmin,iz_latmin) % 获取每日落在br站点4*4位置处的卫星观测经纬度所在行
    if  p > 0 
        for sta = 1 : size(iz_station,1)
            if i == iz_station(sta,2)-59
                iz_modsatsta(i,:) = [ch4_column_model(p),xch4_sat(p),iz_station(sta,5)]; % 按日期整理卫星和模式柱浓度
            end
        end
    end
    % 第5个站（Lamont）三种数据归纳整理
    [p] = sectxy(lon_sat,lat_sat,Lamont_lonmax,Lamont_latmax,Lamont_lonmin,Lamont_latmin) % 获取每日落在br站点4*4位置处的卫星观测经纬度所在行
    if  p > 0 
        for sta = 1 : size(Lamont_station,1)
            if i == Lamont_station(sta,2)-59
                Lamont_modsatsta(i,:) = [ch4_column_model(p),xch4_sat(p),Lamont_station(sta,5)]; % 按日期整理卫星和模式柱浓度
            end
        end
    end
    % 第6个站（tk）三种数据归纳整理
    [p] = sectxy(lon_sat,lat_sat,tk_lonmax,tk_latmax,tk_lonmin,tk_latmin) % 获取每日落在br站点4*4位置处的卫星观测经纬度所在行
    if  p > 0 
        for sta = 1 : size(tk_station,1)
            if i == tk_station(sta,2)-59
                tk_modsatsta(i,:) = [ch4_column_model(p),xch4_sat(p),tk_station(sta,5)]; % 按日期整理卫星和模式柱浓度
            end
        end
    end
    % 第7个站（wg）三种数据归纳整理
    [p] = sectxy(lon_sat,lat_sat,wg_lonmax,wg_latmax,wg_lonmin,wg_latmin) % 获取每日落在br站点4*4位置处的卫星观测经纬度所在行
    if  p > 0 
        for sta = 1 : size(wg_station,1)
            if i == wg_station(sta,2)-59
                wg_modsatsta(i,:) = [ch4_column_model(p),xch4_sat(p),wg_station(sta,5)]; % 按日期整理卫星和模式柱浓度
            end
        end
    end
    
    clear position_mark ch4_slip_model;
end
%% 求相关系数

% ae_modsatsta(ae_modsatsta(:,1)==0,:) = [];
% ae_cof_mod = corrcoef(ae_modsatsta(:,1),ae_modsatsta(:,3)*1000);
% ae_cof_sat = corrcoef(ae_modsatsta(:,2),ae_modsatsta(:,3)*1000);

br_modsatsta(br_modsatsta(:,1)==0,:) = [];
br_cof_mod = corrcoef(br_modsatsta(:,1),br_modsatsta(:,3)*1000);
br_cof_sat = corrcoef(br_modsatsta(:,2),br_modsatsta(:,3)*1000);

db_modsatsta(db_modsatsta(:,1)==0,:) = [];
db_cof_mod = corrcoef(db_modsatsta(:,1),db_modsatsta(:,3)*1000);
db_cof_sat = corrcoef(db_modsatsta(:,2),db_modsatsta(:,3)*1000);

iz_modsatsta(iz_modsatsta(:,1)==0,:) = [];
iz_cof_mod = corrcoef(iz_modsatsta(:,1),iz_modsatsta(:,3)*1000);
iz_cof_sat = corrcoef(iz_modsatsta(:,2),iz_modsatsta(:,3)*1000);

Lamont_modsatsta(Lamont_modsatsta(:,1)==0,:) = [];
Lamont_cof_mod = corrcoef(Lamont_modsatsta(:,1),Lamont_modsatsta(:,3)*1000);
Lamont_cof_sat = corrcoef(Lamont_modsatsta(:,2),Lamont_modsatsta(:,3)*1000);

tk_modsatsta(tk_modsatsta(:,1)==0,:) = [];
tk_cof_mod = corrcoef(tk_modsatsta(:,1),tk_modsatsta(:,3)*1000);
tk_cof_sat = corrcoef(tk_modsatsta(:,2),tk_modsatsta(:,3)*1000);

wg_modsatsta(wg_modsatsta(:,1)==0,:) = [];
wg_cof_mod = corrcoef(wg_modsatsta(:,1),wg_modsatsta(:,3)*1000);
wg_cof_sat = corrcoef(wg_modsatsta(:,2),wg_modsatsta(:,3)*1000);




















