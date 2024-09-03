%% Conversion of dataset from Multi Observation Global Ocean to cell array
clearvars;
clc;
t = double(ncread("J4702_13_2deg.nc",'time')) ; % time = 1 week, seconds since 1970
lat = ncread("J4702_13_2deg.nc",'latitude'); % in degrees north
lon = ncread("J4702_13_2deg.nc",'longitude'); % in degrees east
depth = double(ncread("J4702_13_2deg.nc",'depth')); % in m
mixlay = ncread("J4702_13_2deg.nc",'mlotst'); % in m
temp = ncread("J4702_13_2deg.nc",'to'); % in °C
ut = ncread("J4702_13_2deg.nc",'ugo'); %in m/s
vt = ncread("J4702_13_2deg.nc",'vgo'); %in m/s

%% Example runs
% Current prediction
cur = curpred(lon,lat,depth,t,ut,vt,-20.84,47.67,100,1100,14,1995);
writetable(cur,"J4702_13_curpred.xlsx");

% Temperature at origin
T_total = temppred(lon,lat,depth,mixlay,t,ut,vt,temp,-20.84,47.67,0,100,100,1100,22,14);
% Average temperature
T_totalavg = temppred_avg(lon,lat,depth,mixlay,t,ut,vt,temp,-20.84,47.67,0,100,100,1100,22,14,1);
% Temperature at sediment trap location
T_total2 = temppred2(lon,lat,depth,mixlay,t,temp,-20.84,47.67,0,100,100,1100,22,14);