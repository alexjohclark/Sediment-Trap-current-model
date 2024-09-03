%% Function deriving the sinking path of a particle and temperature from origin
function te = temppred(lon,lat,depth,mixlay,t,ut,vt,temp,long,lati,mini,maxi,S,Tdep,nrsamp,day)
% Velocity for all locations into a more workable form
utot = zeros(height(depth),height(t),height(lon),height(lat));
vtot = zeros(height(depth),height(t),height(lon),height(lat));
for i = 1:height(lon)
    for ii = 1:height(lat)
        utot(:,:,i,ii) = reshape(ut(i,ii,:,:),height(depth),height(t)); %depth = rows, time is columns, last two are lon and lat
        vtot(:,:,i,ii) = reshape(vt(i,ii,:,:),height(depth),height(t));
    end
end 
utot = flip(utot,1);
vtot = flip(vtot,1);

%% Time duration of time opening
Ss = S/86400; % is the sinking speed in m/s (from m/d)
dep = Ss*60*60; % depth after an hour for faster computations, m/h
ttrap = Tdep/dep; %time to trap in hours
tweek = t/60/60/24/7;
dur = day/7; %in weeks
Sw = S*7;
twtrap = Tdep/Sw; %time to trap in weeks
sttrap = tweek(end)-twtrap-dur;
if sttrap <= tweek(1)
    optrap = 1;
else
    optrap = discretize(sttrap,tweek);
end
entrap = tweek(end);
cltrap = discretize(entrap,tweek);

z = 0;
for iii = 1 : ttrap
    z(iii+1) = z(iii) + dep;
end
z_column = transpose(z);

ddd = discretize(z_column,depth,'IncludedEdge','left');

x = -60*60*utot(ddd,:,:,:); %speed at each depth interval for duration of trap in hours
y = -60*60*vtot(ddd,:,:,:);
%% Correction if NaN exists
x(isnan(x))=0;
y(isnan(y))=0;

%% Distance to travel per time inverval and location
xt = zeros(height(x),width(x),height(lon),height(lat));
yt = zeros(height(x),width(y),height(lon),height(lat));
for iv = 1:height(lon)
    for v = 1:height(lat)
        for j = optrap : cltrap
            for ji = 1 : length(x)
                xt(ji+1,j,iv,v) = sum(x(1:ji,j,iv,v));
                yt(ji+1,j,iv,v) = sum(y(1:ji,j,iv,v));
                xt(1,j,iv,v)=0;
                yt(1,j,iv,v)=0;
            end
        end
    end
end

%% Movement of a particle
% Reshape into lon-lat vs depth and time
startx = discretize(long,lon);
starty = discretize(lati,lat);

R = 6371000;%depth of the Earth in m
xcoor = R.*cosd(lat).*cosd(lon); %x-coordinates of lon-lat grid
ycoor = R.*cosd(lat).*sind(lon); %y-coordinates of lon-lat grid
for vi = 1:length(xcoor)/2 - 1
    for vii = 1:length(xt)
        for ti = optrap : cltrap
            if abs(xcoor(startx) + xt(vii,ti,startx,starty)) > abs(xcoor(startx+vi)) || abs(xcoor(startx) - xt(vii,ti,startx,starty)) < abs(xcoor(startx-vi))
                xm(vii,:,:,:) = xt(vii,ti,startx,starty);
            elseif abs(xcoor(startx) + xt(vii,ti,startx,starty)) <= abs(xcoor(startx+vi))
                xm(vii,:,:,:) = xt(vii,ti,startx+vi,starty);
            elseif abs(xcoor(startx) - xt(vii,ti,startx,starty)) >= abs(xcoor(startx-vi))
                xm(vii,:,:,:) = xt(vii,ti,startx-vi,starty);
            end
        end
    end
end
for viii = 1:length(ycoor)/2 - 1
    for ix = 1:length(yt)
        for ti = optrap : cltrap
            if abs(ycoor(starty) + yt(ix,ti,startx,starty)) > abs(ycoor(starty+viii)) || abs(ycoor(starty) - yt(ix,ti,startx,starty)) < abs(ycoor(starty-viii))
                ym(ix,:,:,:) = yt(ix,ti,startx,starty);
            elseif abs(ycoor(starty) + yt(ix,ti,startx,starty)) <= abs(ycoor(starty+viii))
                ym(ix,:,:,:) = yt(ix,ti,startx,starty+viii);
            elseif abs(ycoor(starty) - yt(ix,ti,startx,starty)) >= abs(ycoor(starty-viii))
                ym(ix,:,:,:) = yt(ix,ti,startx,starty-viii);
            end
        end
    end
end
%% Corrections to previous
for jii = 1:length(xt)-2
    for ji = jii + 1
        if abs(xm(ji)) < abs(xm(jii))
            xm(ji) = xm(jii) + (xm(ji+1)-xm(ji));
        elseif isnan(abs(xm(ji))) && isnan(abs(xm(ji+1)))
            xm(ji) = xm(jii);
        else
            xm(jii) = xm(jii);
        end
    end
end
for jiii = 1:length(yt)-2
        if abs(ym(jiii+1)) < abs(ym(jiii))
            ym(jiii+1) = ym(jiii) + (ym(jiii+2)-ym(jiii+1));
        elseif isnan(abs(ym(jiii+1))) && isnan(abs(ym(jiii+2)))
            ym(jiii+1) = ym(jiii);
        else
            ym(jiii) = ym(jiii);
        end
end
% Remove last row as it won't work with the above
xm(end) = [];
ym(end) = [];
% Total distance travelled
for i = 1:length(xm)-1
    len(i) = sqrt((xm(i+1)-xm(i))^2 + (ym(i+1)-ym(i))^2 + (z_column(i+1)-z_column(i))^2);
    lentot(i+1) = sum(len(1:i));
end
lengt = transpose(lentot);
ylat = (ym./110574.0)+lati;
xlon = xm.*(1/(111320*cosd(lati)))+long;

%% Temperature at location
% Temperature range
temptot = zeros(length(depth),length(t),length(lon),length(lat));
temp = flip(flip(temp,1),2);
for i = 1:length(lon)
    for ii = 1:length(lat)
        temptot(:,:,i,ii) = reshape(temp(i,ii,:,:),length(depth),height(t)); %depth = rows, time is columns, last two are lon and lat
    end
end

depthmin = discretize(mini,depth);
depthmax = discretize(maxi,depth);
locx = xlon(end);
locy = ylat(end);
%Correction if current path goes outside of model dataset - will be visible
%on plots to then obtain a wider dataset
if locx >= lon(end)
    actlocx = discretize(lon(end),lon);
elseif locx <= lon(1)
    actlocx = discretize(lon(1),lon);
else
    actlocx = discretize(locx,lon);
end
if locy >= lat(end)
    actlaty = discretize(lat(end),lat);
elseif locy <= lat(1)
    actlaty = discretize(lat(1),lat);
else
    actlaty = discretize(locy,lat);
end

%% Above mixed layer depth capture cone for open time of trap at provenance
mixlay(isnan(mixlay))=0;
for i = -1 : 1 : 1
    for ii = -1 : 1 : 1
        for iii = i + 2
            for iv = ii + 2
                for v = optrap : cltrap
                    mixdeptest(iii,iv,v) = discretize(mixlay(actlocx + i,actlaty + ii,v),depth);
                    mixdep(iii,iv,v) = mixlay(actlocx + i,actlaty + ii,v);
                    for vi = 1 : mixdeptest(iii,iv,v)
                    mixtemper(vi,v,iii,iv) = temptot(vi, v, actlocx + i,actlaty + ii);
                    mix_dep_tempt(iii,iv,v) = mean(nonzeros(mixtemper(:,v,iii,iv)),'omitnan');
                    mix_dep_tempts(iii,iv,v) = std(nonzeros(mixtemper(:,v,iii,iv)),'omitnan');
                    end
                end
            end
        end
    end
end
mix_dep_tempt(mix_dep_tempt==0)=NaN;
mix_dep_tempts(mix_dep_tempts==0)=NaN;
% Average and stdev for capture cone around provenance for entire time
mix_dep_temp_mean = mean(mix_dep_tempt,3,'omitnan');
mix_dep_temp_std = std(mix_dep_tempt,0,3,'omitnan');

mix_dep_temp = mean(nonzeros(mix_dep_temp_mean),'omitnan');
std_dep_temp = mean(nonzeros(mix_dep_temp_std),'omitnan');

%% Maximum production depth capture cone for open time of trap at provenance
for i = -1 : 1 : 1
    for ii = -1 : 1 : 1
        for iii = i + 2
            for iv = ii + 2
                for v = optrap : cltrap
                    for vi = depthmin : depthmax
                    maxtemper(vi,v,iii,iv) = temptot(vi, v, actlocx + i,actlaty + ii);
                    max_prod_tempt(iii,iv,v) = mean(nonzeros(maxtemper(:,v,iii,iv)),'omitnan');
                    max_prod_tempts(iii,iv,v) = std(nonzeros(maxtemper(:,v,iii,iv)),'omitnan');
                    end
                end
            end
        end
    end
end
max_prod_tempt(max_prod_tempt==0)=NaN;
max_prod_tempts(max_prod_tempts==0)=NaN;
% Average and stdev for capture cone around provenance for entire time
max_prod_temp_mean = mean(max_prod_tempt,3,'omitnan');
max_prod_temp_std = std(max_prod_tempt,0,3,'omitnan');

max_prod_temp = mean(nonzeros(max_prod_temp_mean),'omitnan');
std_prod_temp = mean(nonzeros(max_prod_temp_std),'omitnan');

%% Temperature for whole photic zone 0-200 m capture cone at provenance
for i = -1 : 1 : 1
    for ii = -1 : 1 : 1
        for iii = i + 2
            for iv = ii + 2
                for v = optrap : cltrap
                    for vi = 1 : 22
                    photemper(vi,v,iii,iv) = temptot(vi, v, actlocx + i,actlaty + ii);
                    phot_zo_tempt(iii,iv,v) = mean(nonzeros(photemper(:,v,iii,iv)),'omitnan');
                    phot_zo_tempts(iii,iv,v) = std(nonzeros(photemper(:,v,iii,iv)),'omitnan');
                    end
                end
            end
        end
    end
end
phot_zo_tempt(phot_zo_tempt==0)=NaN;
phot_zo_tempts(phot_zo_tempts==0)=NaN;
% Average and stdev for capture cone around provenance for entire time
phot_zo_temp_mean = mean(phot_zo_tempt,3,'omitnan');
phot_zo_temp_std = std(phot_zo_tempt,0,3,'omitnan');

max_prod_temp1 = mean(nonzeros(phot_zo_temp_mean),'omitnan');
std_prod_temp1 = mean(nonzeros(phot_zo_temp_std),'omitnan');

%% Sea Surface Temperature capture cone at provenance
for i = -1 : 1 : 1
    for ii = -1 : 1 : 1
        for iii = i + 2
            for iv = ii + 2
                for v = optrap : cltrap
                    for vi = 1 : 2
                    ssttemper(vi,v,iii,iv) = temptot(vi, v, actlocx + i,actlaty + ii);
                    sst_tempt(iii,iv,v) = mean(nonzeros(ssttemper(:,v,iii,iv)),'omitnan');
                    sst_tempts(iii,iv,v) = std(nonzeros(ssttemper(:,v,iii,iv)),'omitnan');
                    end
                end
            end
        end
    end
end
sst_tempt(sst_tempt==0)=NaN;
sst_tempts(sst_tempts==0)=NaN;
% Average and stdev for capture cone around provenance for entire time
sst_temp_mean = mean(sst_tempt,3,'omitnan');
sst_temp_std = std(sst_tempt,0,3,'omitnan');

meantempsst = mean(nonzeros(sst_temp_mean),'omitnan');
stdtempsst = mean(nonzeros(sst_temp_std),'omitnan');

%% Depth range of D47 temperature at provenance
D47_temptable = readtable("Supplementary Table S7.csv",'VariableNamingRule','preserve');
D47_temp = table2array(D47_temptable(:,'D47 Temp'));
D47_temp_SD = table2array(D47_temptable(:,'D47 Temp SD'));
%nrsamp = 18;
depthD47 = zeros(length(depth),ceil(dur),length(lon),length(lat));
for j = 1:length(depth)
    for ji = optrap:cltrap
        if temptot(j,ji,actlocx,actlaty) <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && temptot(j,ji,actlocx,actlaty) >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
            depthD47(j,ji,:,:) = depth(j);
        elseif temptot(j,ji,actlocx,actlaty) <= D47_temp(nrsamp) && temptot(j,ji,actlocx,actlaty) >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
            depthD47(j,ji,:,:) = depth(j);
        elseif temptot(j,ji,actlocx,actlaty) <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && temptot(j,ji,actlocx,actlaty) >= D47_temp(nrsamp)
            depthD47(j,ji,:,:) = depth(j);
        else 
            depthD47(j,ji,:,:) = NaN;
        end
            depthD47range(1,ji) = min(depthD47(:,ji,actlocx,actlaty),[],'omitnan');
            depthD47range(2,ji) = max(depthD47(:,ji,actlocx,actlaty),[],'omitnan');
    end
end
mdepD47 = mean(depthD47range,2,'omitnan');
%% Does D47 temperature fit with which temperature range at provenance
if mix_dep_temp + std_dep_temp <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && mix_dep_temp + std_dep_temp >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,1) = mix_dep_temp;
    D47_comp(2,1) = std_dep_temp;
    D47_comp(3,1) = 1;
    D47_comp(4,1) = locx;
    D47_comp(5,1) = locy;
    D47_comp(6,1) = mean(mean(mean(mixdep,'omitnan'),'omitnan'),'omitnan');
elseif mix_dep_temp - std_dep_temp <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && mix_dep_temp - std_dep_temp >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,1) = mix_dep_temp;
    D47_comp(2,1) = std_dep_temp;
    D47_comp(3,1) = 1;
    D47_comp(4,1) = locx;
    D47_comp(5,1) = locy;
    D47_comp(6,1) = mean(mean(mean(mixdep,'omitnan'),'omitnan'),'omitnan');
else
    D47_comp(1,1) = mix_dep_temp;
    D47_comp(2,1) = std_dep_temp;
    D47_comp(3,1) = 0;
    D47_comp(4,1) = locx;
    D47_comp(5,1) = locy;
    D47_comp(6,1) = mean(mean(mean(mixdep,'omitnan'),'omitnan'),'omitnan');
end
if max_prod_temp + std_prod_temp <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && max_prod_temp + std_prod_temp >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,2) = max_prod_temp;
    D47_comp(2,2) = std_prod_temp;
    D47_comp(3,2) = 1;
    D47_comp(4,2) = locx;
    D47_comp(5,2) = locy;
elseif max_prod_temp - std_prod_temp <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && max_prod_temp - std_prod_temp >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,2) = max_prod_temp;
    D47_comp(2,2) = std_prod_temp;
    D47_comp(3,2) = 1;
    D47_comp(4,2) = locx;
    D47_comp(5,2) = locy;
else
    D47_comp(1,2) = max_prod_temp;
    D47_comp(2,2) = std_prod_temp;
    D47_comp(3,2) = 0;
    D47_comp(4,2) = locx;
    D47_comp(5,2) = locy;
end
if max_prod_temp1 + std_prod_temp1 <= D47_temp(nrsamp) + D47_temp_SD(nrsamp)  && max_prod_temp1 + std_prod_temp1 >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,3) = max_prod_temp1;
    D47_comp(2,3) = std_prod_temp1;
    D47_comp(3,3) = 1;
    D47_comp(4,3) = locx;
    D47_comp(5,3) = locy;
elseif max_prod_temp1 - std_prod_temp1 <= D47_temp(nrsamp) + D47_temp_SD(nrsamp)  && max_prod_temp1 - std_prod_temp1 >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,3) = max_prod_temp1;
    D47_comp(2,3) = std_prod_temp1;
    D47_comp(3,3) = 1;
    D47_comp(4,3) = locx;
    D47_comp(5,3) = locy;
else
    D47_comp(1,3) = max_prod_temp1;
    D47_comp(2,3) = std_prod_temp1;
    D47_comp(3,3) = 0;
    D47_comp(4,3) = locx;
    D47_comp(5,3) = locy;
end
if meantempsst + stdtempsst<= D47_temp(nrsamp) + D47_temp_SD(nrsamp)  && meantempsst + stdtempsst >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,4) = meantempsst;
    D47_comp(2,4) = stdtempsst;
    D47_comp(3,4) = 1;
    D47_comp(4,4) = locx;
    D47_comp(5,4) = locy;
elseif meantempsst - stdtempsst <= D47_temp(nrsamp) + D47_temp_SD(nrsamp)  && meantempsst - stdtempsst >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,4) = meantempsst;
    D47_comp(2,4) = stdtempsst;
    D47_comp(3,4) = 1;
    D47_comp(4,4) = locx;
    D47_comp(5,4) = locy;
else
    D47_comp(1,4) = meantempsst;
    D47_comp(2,4) = stdtempsst;
    D47_comp(3,4) = 0;
    D47_comp(4,4) = locx;
    D47_comp(5,4) = locy;
end
D47_comp(1,5) = D47_temp(nrsamp);
D47_comp(2,5) = D47_temp_SD(nrsamp);
D47_comp(3,5) = 1;
D47_comp(4,5) = mdepD47(1);
D47_comp(5,5) = mdepD47(2);

te = array2table(D47_comp,'VariableNames',{'Above mixed layer Temp (°C)','Maximum Production Temp (°C)', 'Photic Zone Temp (°C)', 'Sea Surface Temp (°C)', 'D47 Temp (°C)'});
%% For assembling and exporting into a matrix and tables
for i = 1 : 3
    for ii = 1 : 3
        for iii = ii + 3
            for iv = i + 3
                for v = iii + 3
               D47_mix(i,ii) = mix_dep_tempt(i,ii,optrap);
               D47_max(i,ii) = max_prod_tempt(i,ii,optrap);
               D47_phot(i,ii) = phot_zo_tempt(i,ii,optrap);
               D47_sst(i,ii) = sst_tempt(i,ii,optrap);
               D47_mix(i,iii) = mix_dep_temp_mean(i,ii);
               D47_max(i,iii) = max_prod_temp_mean(i,ii);
               D47_phot(i,iii) = phot_zo_temp_mean(i,ii);
               D47_sst(i,iii) = sst_temp_mean(i,ii);
               D47_mix(i,v) = mix_dep_tempt(i,ii,cltrap);
               D47_max(i,v) = max_prod_tempt(i,ii,cltrap);
               D47_phot(i,v) = phot_zo_tempt(i,ii,cltrap);
               D47_sst(i,v) = sst_tempt(i,ii,cltrap);
               D47_mix(iv,ii) = mix_dep_tempts(i,ii,optrap);
               D47_max(iv,ii) = max_prod_tempts(i,ii,optrap);
               D47_phot(iv,ii) = phot_zo_tempts(i,ii,optrap);
               D47_sst(iv,ii) = sst_tempts(i,ii,optrap);
               D47_mix(iii,iv) = mix_dep_temp_std(i,ii);
               D47_max(iii,iv) = max_prod_temp_std(i,ii);
               D47_phot(iii,iv) = phot_zo_temp_std(i,ii);
               D47_sst(iii,iv) = sst_temp_std(i,ii);
               D47_mix(iv,v) = mix_dep_tempts(i,ii,cltrap);
               D47_max(iv,v) = max_prod_tempts(i,ii,cltrap);
               D47_phot(iv,v) = phot_zo_tempts(i,ii,cltrap);
               D47_sst(iv,v) = sst_tempts(i,ii,cltrap);
                end
            end
        end
    end
end

D47_mix1 = array2table(D47_mix);
D47_max1 = array2table(D47_max);
D47_phot1 = array2table(D47_phot);
D47_sst1 = array2table(D47_sst);
writetable(D47_mix1,'temps_mixlay.xlsx');
writetable(D47_max1,'temps_maxpro.xlsx');
writetable(D47_phot1,'temps_phozo.xlsx');
writetable(D47_sst1,'temps_SST.xlsx');

%% Reconfigure for figure plotting
temptot1 = mean(mean(temptot(depthmin:depthmax,optrap:cltrap,:,:),'omitnan'),'omitnan');
temptot1 = flip(flip(temptot1,1),2);
for i = 1:length(lon)
    for ii = 1:length(lat)
        temptot2(i,ii,:,:) = reshape(temptot1(:,:,i,ii),1,1); %depth = rows, time is columns, last two are lon and lat
    end
end
% Figure to double check
figure
s = pcolor(lon,lat,temptot2);
s.FaceColor = 'interp';
hold
plot(xlon,ylat,"Color",'red')
grid on
xlabel('Longitude (°)')
ylabel('Latitude (°)')
scatter(long,lati,"filled",'o')
text(xlon(1)+0.05,ylat(1)+0.05,'Sediment Trap')
c = colorbar;
c.Label.String = 'Average maximum production Temp (°C)';
end