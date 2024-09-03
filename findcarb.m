%% Function to derive carbonate parameters around given location
function d13C = findcarb(londata,latdata,depthmin,depthmax)
load GLODAPv2.2023_Merged_Master_File.mat

lon = G2longitude;
lat = G2latitude;
c13 = G2c13;
o18 = G2o18;
pH = G2phtsinsitutp;
Alk = G2talk;
year = G2year;
month = G2month;
day = G2day;
depth = G2depth;
%% Define year and months into matrix
yelim = year;
molim = month;
time = [yelim molim];
%% Defines depths
%depthmin = 0;
%depthmax = 150;
deplim = 0;
for i = 1:length(depth)
    if depth(i) <= depthmax && depth(i) >= depthmin
        deplim(i) = depth(i);
    else
        deplim(i) = NaN;
    end
end
deplim = deplim.';
%% Gives area around set location
%londata = 58.59;
%latdata = 19.39;
lon_lim_max = londata + 10;
lon_lim_min = londata - 10;
lat_lim_max = latdata + 10;
lat_lim_min = latdata - 10;
lonlim = 0;
for i = 1:length(lon)
    if lon(i) <= lon_lim_max && lon(i) >= lon_lim_min
        lonlim(i) = lon(i);
    else
        lonlim(i) = NaN;
    end
end
lonlim = lonlim.';
latlim = 0;
for i = 1:length(lat)
    if lat(i) <= lat_lim_max && lat(i) >= lat_lim_min
        latlim(i) = lat(i);
    else
        latlim(i) = NaN;
    end
end
latlim = latlim.';
%% Arranges into a matrix
data = [time lonlim latlim deplim c13 o18 pH Alk];
dat2 = data;
dat2(any(isnan(data(:,1:5)),2),:) = [];
%% Finds and defines values for given depths
few2 = dat2;
few2(any(isnan(dat2(:,6)),2),:) = [];
for i = 1:length(few2)
    if few2(i,5) <= 50
        fe1(i,:) = few2(i,:);
    elseif few2(i,5) > 50 && few2(i,5) <= 100
        fe2(i,:) = few2(i,:);
    elseif few2(i,5) > 100 && few2(i,5) <= 150
        fe3(i,:) = few2(i,:);
    else
        fe1(i,:) = NaN;
        fe2(i,:) = NaN;
        fe3(i,:) = NaN;
    end
end

for i = 1 : length(fe1)
    if fe1(i,1) == 0
        fe1(i,:) = NaN;
    else
        fe1(i,:) = fe1(i,:);
    end
end
dC0_50 = fe1;
dC0_50(any(isnan(fe1(:,1)),2),:) = [];
for ii = 1 : length(fe2)
    if fe2(ii,1) == 0
        fe2(ii,:) = NaN;
    else
        fe2(ii,:) = fe2(ii,:);
    end
end
dC50_100 = fe2;
dC50_100(any(isnan(fe2(:,1)),2),:) = [];
for iii = 1 : length(fe3)
    if fe3(iii,1) == 0
        fe3(iii,:) = NaN;
    else
        fe3(iii,:) = fe3(iii,:);
    end
end
dC100_150 = fe3;
dC100_150(any(isnan(fe3(:,1)),2),:) = [];

dat3 = [dC0_50;dC50_100;dC100_150];
%% Arrange data
dat4(1,:) = mean(dC0_50(:,5:9),'omitnan');
dat4(2,:) = std(dC0_50(:,5:9),'omitnan');
dat4(3,:) = mean(dC50_100(:,5:9),'omitnan');
dat4(4,:) = std(dC50_100(:,5:9),'omitnan');
dat4(5,:) = mean(dC100_150(:,5:9),'omitnan');
dat4(6,:) = std(dC100_150(:,5:9),'omitnan');
%% Export data
d13C = array2table(dat4,'VariableNames',{'Depth','d13CDIC','d18Osw','pH','Alkalinity'});
end