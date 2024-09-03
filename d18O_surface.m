%% d18O of the surface ocean txt file from NASA
%ncdisp("d18O_dataset.nc") %displays file

%% Put everything in Matlab format
clearvars;
clc;
lat = ncread("d18O_dataset.nc",'lat'); % in degrees north
lon = ncread("d18O_dataset.nc",'lon'); % in degrees east
depth = double(ncread("d18O_dataset.nc",'depth')); % in m
d18O = ncread("d18O_dataset.nc",'d18o'); % in m
for i = 1:length(lon)
    for ii = 1:length(lat)
        for iii = 1:length(depth)
            if d18O(i,ii,iii) <= -1.0000000e+29
                d18O(i,ii,iii) = NaN;
            else 
                d18O(i,ii,iii) = d18O(i,ii,iii);
            end
        end
    end
end
%% Selection of area
londata = -20.84; %to input
latdata = 47.67; %to input
lon_lim_max = londata + 5;
lon_lim_min = londata - 5;
lat_lim_max = latdata + 5;
lat_lim_min = latdata - 5;
d18Olonlim = lon_lim_min:1:lon_lim_max;
d18Olonlim1 = discretize(lon,d18Olonlim);
d18Olatlim = lat_lim_min:1:lat_lim_max;
d18Olatlim1 = transpose(discretize(lat,d18Olatlim));
d18Olonlim2 = repmat(d18Olonlim1,1,length(lat));
d18Olatlim2 = repmat(d18Olatlim1,length(lon),1);
d18Oselect = d18Olonlim2.*d18Olatlim2;
d18Oselect1 = repmat(d18Oselect,1,1,length(depth));
d18Oselect2 = ~isnan(d18Oselect1);

%% 0-100 m d18O
d18Odis = d18Oselect2.*d18O;
d18Odis(d18Odis == 0) = NaN;
d18Oupper = d18Odis(:,:,1:7);
meand18Ou = mean(mean(mean(d18Oupper,'omitnan'),'omitnan'),'omitnan');
stdd18Ou = mean(mean(std(d18Oupper,'omitnan'),'omitnan'),'omitnan');
%% Export for one location as matrix
writematrix(meand18Ou,"J4702_13_d18Omean.xlsx");
writematrix(stdd18Ou,"J4702_13_d18Ostd.xlsx");
