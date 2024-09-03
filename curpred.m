%% Function of a sinking particle until depth for given time period and year
function f = curpred(lon,lat,depth,t,ut,vt,long,lati,S,Tdep,day,year)
% Velocity for all locations into a more workable form
utot = zeros(height(depth),height(t),height(lon),height(lat));
vtot = zeros(height(depth),height(t),height(lon),height(lat));
for i = 1:height(lon)
    for ii = 1:height(lat)
        utot(:,:,i,ii) = reshape(ut(i,ii,:,:),height(depth),height(t)); %depth = rows, time is columns, last two are lon and lat
        vtot(:,:,i,ii) = reshape(vt(i,ii,:,:),height(depth),height(t));
    end
end
% Flipped from sediment trap upwards
utot = flip(utot,1);
vtot = flip(vtot,1);

%% Movement of particle in xyz direction for all
%S = 100; Tdep = 4000; day = 17.38; long = 13.5628; lati = -41.1361; year = 1995;%test input parameters;

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

% Correction if NaN exists
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
R = 6371000;%depth of the Earth in m
xcoor = R.*cosd(lat).*cosd(lon); %x-coordinates of lon-lat grid
ycoor = R.*cosd(lat).*sind(lon); %y-coordinates of lon-lat grid
zcoor = R.*sind(lat);

startx = discretize(long,lon);
starty = discretize(lati,lat);
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
        if abs(xm(jii+1)) < abs(xm(jii))
            xm(jii+1) = xm(jii) + (xm(jii+2)-xm(jii+1));
        elseif isnan(abs(xm(jii+1))) && isnan(abs(xm(jii+2)))
            xm(jii+1) = xm(jii);
        else
            xm(jii) = xm(jii);
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
%% 'Time of capture'
opt = t(optrap)/(60*60*24)-((year-1970)*365);
clt = t(cltrap)/(60*60*24)-((year-1970)*365);
if opt >= 0 && opt <= 31
    opt1 = 1;
elseif opt > 31 && opt <= 59
    opt1 = 2;
    opt = opt - 31;
elseif opt > 59 && opt <= 90
    opt1 = 3;
    opt = opt - 59;
elseif opt > 90 && opt <= 120
    opt1 = 4;
    opt = opt - 90;
elseif opt > 120 && opt <= 151
    opt1 = 5;
    opt = opt - 120;
elseif opt > 151 && opt <= 181
    opt1 = 6;
    opt = opt - 151;
elseif opt > 181 && opt <= 212
    opt1 = 7;
    opt = opt - 181;
elseif opt > 212 && opt <= 243
    opt1 = 8;
    opt = opt - 212;
elseif opt > 243 && opt <= 273
    opt1 = 9;
    opt = opt - 243;
elseif opt > 273 && opt <= 304
    opt1 = 10;
    opt = opt - 273;
elseif opt > 304 && opt <= 334
    opt1 = 11;
    opt = opt - 304;
elseif opt > 334 && opt <= 365
    opt1 = 12;
    opt = opt - 334;
end
if clt >= 0 && clt <= 31
    clt1 = 1;
elseif clt > 31 && clt <= 59
    clt1 = 2;
    clt = clt - 31;
elseif clt > 59 && clt <= 90
    clt1 = 3;
    clt = clt - 59;
elseif clt > 90 && clt <= 120
    clt1 = 4;
    clt = clt - 90;
elseif clt > 120 && clt <= 151
    clt1 = 5;
    clt = clt - 120;
elseif clt > 151 && clt <= 181
    clt1 = 6;
    clt = clt - 151;
elseif clt > 181 && clt <= 212
    clt1 = 7;
    clt = clt - 181;
elseif clt > 212 && clt <= 243
    clt1 = 8;
    clt = clt - 212;
elseif clt > 243 && clt <= 273
    clt1 = 9;
    clt = clt - 243;
elseif clt > 273 && clt <= 304
    clt1 = 10;
    clt = clt - 273;
elseif clt > 304 && clt <= 334
    clt1 = 11;
    clt = clt - 304;
elseif clt > 334 && clt <= 365
    clt1 = 12;
    clt = clt - 334;
end
for ii = length(xm)
    findur(ii,1) = opt;
    findur(ii-1,1) = opt1;
    findur(ii,2) = clt;
    findur(ii-1,2) = clt1;
end

f = array2table([xm xlon ym ylat z_column lengt findur],'VariableNames',{'x distance (m)', 'x longitude (°)','y distance (m)', 'y latitude (°)','depth (m)','Total distance from origin (m)', 'Duration colums', 'Duration columns'});

%% Plot A1
%figure
%plot(-ut(:,:,1,1))
%xticklabels(lon)
%xlabel('Longitude (°)')
%ylabel('Eastward current velocity (m/s)')
%% Plot A2
%figure
%plot(-vt(:,:,1,1))
%xticklabels(lon)
%xlabel('Longitude (°)')
%ylabel('Northward current velocity (m/s)')
%% Plot A potential
%ut1 = flip(flip(ut,1),2);
%utest = -ut1(:,:,1,1);
%vt1 = flip(flip(vt,1),2);
%vtest = -vt1(:,:,1,1);
%hyptest = hypot(utest,vtest);
%dirtest = atan2d(vtest,utest);
%figure
%pcolor(lon,lat,hyptest)
%xlabel('Longitude (°)')
%ylabel('Latitude (°)')
%c = colorbar;
%c.Label.String = 'Current velocity magnitude (m/s)';
%figure
%pcolor(lon,lat,dirtest)
%xlabel('Longitude (°)')
%ylabel('Latitude (°)')
%c = colorbar;
%c.Label.String = 'Current velocity direction (°)';
%% plot B
%utest = reshape(-ut(startx,starty,:,1),length(depth),1);
%vtest = reshape(-vt(startx,starty,:,1),length(depth),1);
%hyptest = hypot(utest,vtest);
%figure
%plot(hyptest,depth)
%xlabel('Current velocity magnitude (m/s)')
%ylabel('Depth (m)')
%% Plot C
%utest = reshape(-ut(startx,starty,:,1),length(depth),1);
%vtest = reshape(-vt(startx,starty,:,1),length(depth),1);
%hyptest = hypot(utest,vtest);
%figure
%plot(utest,depth)
%xlabel('Northward current velocity (m/s)')
%ylabel('Depth (m)')
%yline(z_column) %change to z_column in curpred
%figure
%plot(vtest,depth)
%xlabel('Eastward current velocity (m/s)')
%ylabel('Depth (m)')
%yline(z_column) %change to z_column in curpred
%figure
%plot(hyptest(1:22),depth(1:22))
%xlabel('Current velocity magnitude (m/s)')
%ylabel('Depth (m)')
%yline(z_column(1:49)) %change to z_column in curpred
%% Plot D
%figure
%plot(xm,z_column)
%xlabel('Total eastward movement (m)')
%ylabel('Depth (m)')
%figure
%plot(ym,z_column)
%xlabel('Total northward movement (m)')
%ylabel('Depth (m)')
%figure
%plot(xm,ym);
%grid on
%set(gca,'XAxisLocation','origin','YAxisLocation', 'origin');
%xlabel('Eastward movement (m)')
%ylabel('Northward movement (m)')
%% Plot E, path of particle
figure
plot(xlon,ylat)
grid on
xlabel('Longitude (°)')
ylabel('Latitude (°)')
hold
scatter(xlon(1),ylat(1),"filled",'o')
scatter(xlon(end),ylat(end),"filled",'d')
text(xlon(1)-1.15,ylat(1),'Sediment Trap')
text(xlon(end)+0.15,ylat(end)+0.05,'Provenance point')
text(xlon(end)+0.03,ylat(end)-0.03,[num2str(lengt(end))])
end
