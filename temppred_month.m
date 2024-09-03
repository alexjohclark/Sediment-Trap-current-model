%% Function to find what months fit D47 temperatures
function te = temppred_month(lon,lat,depth,t,temp,long,lati,mini,maxi,nrsamp)
%% Temperature
% Temperature range for maximum production depths
temptot = zeros(length(depth),length(t),length(lon),length(lat));
temp = flip(flip(temp,1),2);
for i = 1:length(lon)
    for ii = 1:length(lat)
        temptot(:,:,i,ii) = reshape(temp(i,ii,:,:),length(depth),height(t)); %depth = rows, time is columns, last two are lon and lat
    end
end

startx = discretize(long,lon);
starty = discretize(lati,lat);

%% Depth range of D47 temperature at sedi trap
D47_temptable = readtable("Supplementary Table S7.csv",'VariableNamingRule','preserve');
D47_temp = table2array(D47_temptable(:,'D47 Temp'));
D47_temp_SD = table2array(D47_temptable(:,'D47 Temp SD'));
%nrsamp = 7;
depthD47 = zeros(length(depth),length(t),length(lon),length(lat));
depthD47range = zeros(2,length(t));
for j = 1:length(depth)
    for ji = 1 : length(t)
        if temptot(j,ji,startx,starty) <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && temptot(j,ji,startx,starty) >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
            depthD47(j,ji,:,:) = depth(j);
        elseif temptot(j,ji,startx,starty) <= D47_temp(nrsamp) && temptot(j,ji,startx,starty) >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
            depthD47(j,ji,:,:) = depth(j);
        elseif temptot(j,ji,startx,starty) <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && temptot(j,ji,startx,starty) >= D47_temp(nrsamp)
            depthD47(j,ji,:,:) = depth(j);
        else 
            depthD47(j,ji,:,:) = NaN;
        end
            depthD47range(1,ji) = min(depthD47(:,ji,startx,starty),[],'omitnan');
            depthD47range(2,ji) = max(depthD47(:,ji,startx,starty),[],'omitnan');
    end
end

for i = 1 : length(t)
    if depthD47range(1,i) > mini && depthD47range(2,i) < maxi
        mdepD47(i,1) = depthD47range(1,i);
        mdepD47(i,2) = depthD47range(2,i);
    elseif depthD47range(1,i) < mini && depthD47range(2,i) < maxi
        mdepD47(i,1) = depthD47range(1,i);
        mdepD47(i,2) = depthD47range(2,i);
    elseif depthD47range(1,i) < maxi && depthD47range(2,i) < maxi
        mdepD47(i,1) = depthD47range(1,i);
        mdepD47(i,2) = depthD47range(2,i);
    elseif depthD47range(1,i) > 200 && depthD47range(2,i) > 200
        mdepD47(i,1) = depthD47range(1,i);
        mdepD47(i,2) = depthD47range(2,i);
    else
        mdepD47(i,1) = NaN;
        mdepD47(i,2) = NaN;
    end
end
%% Export
mdep(1,:) = discretize(depthD47range(1,:),depth);
mdep(2,:) = discretize(depthD47range(2,:),depth);
for i = 1 : length(t)
mdepD47(i,3) = mean(temptot(~isnan(mdep(1:2,i)),i,startx,starty),'omitnan');
mdepD47(i,4) = std(temptot(~isnan(mdep(1:2,i)),i,startx,starty),'omitnan');
end

te = array2table(mdepD47,'VariableNames',{'Min depth (m)','Max depth (m)','Average Temperature (°C)','Standard deviation Temperature (°C)'});
end