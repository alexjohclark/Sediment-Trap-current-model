%% xyz of a particle following currents until trap depth
function te = temppred2(lon,lat,depth,mixlay,t,ut,vt,temp,long,lati,mini,maxi,S,Tdep,nrsamp,day)


%S = 100; Tdep = 1996; day = 8.50; mini =20; maxi=80; long = 58.80; lati = 17.40; %input parameters
%%
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

%% Temperature at location and for maximum production depth
% Reshape into lon-lat vs depth and time
startx = discretize(long,lon);
starty = discretize(lati,lat);

%% Temperature
% Temperature range for maximum production depths
temptot = zeros(length(depth),length(t),length(lon),length(lat));
temp = flip(flip(temp,1),2);
for i = 1:length(lon)
    for ii = 1:length(lat)
        temptot(:,:,i,ii) = reshape(temp(i,ii,:,:),length(depth),height(t)); %depth = rows, time is columns, last two are lon and lat
    end
end

depthmin = discretize(mini,depth);
depthmax = discretize(maxi,depth);

%% Above mixed layer depth capture cone for open time of trap at sedi trap
mixlay(isnan(mixlay))=0;
for i = -1 : 1 : 1
    for ii = -1 : 1 : 1
        for iii = i + 2
            for iv = ii + 2
                for v = optrap : cltrap
                    mixdeptest(iii,iv,v) = discretize(mixlay(startx + i,starty + ii,v),depth);
                    mixdep(iii,iv,v) = mixlay(startx + i,starty + ii,v);
                    for vi = 1 : mixdeptest(iii,iv,v)
                    mixtemper(vi,v,iii,iv) = temptot(vi, v, startx + i,starty + ii);
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

%% Maximum production depth capture cone for open time of trap at sedi trap
for i = -1 : 1 : 1
    for ii = -1 : 1 : 1
        for iii = i + 2
            for iv = ii + 2
                for v = optrap : cltrap
                    for vi = depthmin : depthmax
                    maxtemper(vi,v,iii,iv) = temptot(vi, v, startx + i,starty + ii);
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

%% Temperature for whole photic zone 0-200 m capture cone at sedi trap
for i = -1 : 1 : 1
    for ii = -1 : 1 : 1
        for iii = i + 2
            for iv = ii + 2
                for v = optrap : cltrap
                    for vi = 1 : 22
                    photemper(vi,v,iii,iv) = temptot(vi, v, startx + i,starty + ii);
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

%% Sea Surface Temperature capture cone at sedi trap
for i = -1 : 1 : 1
    for ii = -1 : 1 : 1
        for iii = i + 2
            for iv = ii + 2
                for v = optrap : cltrap
                    for vi = 1 : 2
                    ssttemper(vi,v,iii,iv) = temptot(vi, v, startx + i,starty + ii);
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

%% Depth range of D47 temperature at sedi trap
D47_temptable = readtable("D47_sed_simpl.csv",'VariableNamingRule','preserve');
D47_temp = table2array(D47_temptable(:,'D47 Temp'));
D47_temp_SD = table2array(D47_temptable(:,'D47 Temp SD'));
%nrsamp = 7;
depthD47 = zeros(length(depth),ceil(dur),length(lon),length(lat));
for j = 1:length(depth)
    for ji = optrap:cltrap
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
mdepD47 = mean(depthD47range,2,'omitnan');
%sdepD47 = [std(depthD47range(1,:)) std(depthD47range(2,:))]; % Potentially for the temporal range 
%% Does D47 temperature fit with which temperature range at sedi trap
if mix_dep_temp + std_dep_temp <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && mix_dep_temp + std_dep_temp >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,1) = mix_dep_temp;
    D47_comp(2,1) = std_dep_temp;
    D47_comp(3,1) = 1;
    D47_comp(4,1) = long;
    D47_comp(5,1) = lati;
    D47_comp(6,1) = mean(mean(mean(mixdep,'omitnan'),'omitnan'),'omitnan');
elseif mix_dep_temp - std_dep_temp <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && mix_dep_temp - std_dep_temp >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,1) = mix_dep_temp;
    D47_comp(2,1) = std_dep_temp;
    D47_comp(3,1) = 1;
    D47_comp(4,1) = long;
    D47_comp(5,1) = lati;
    D47_comp(6,1) = mean(mean(mean(mixdep,'omitnan'),'omitnan'),'omitnan');
else
    D47_comp(1,1) = mix_dep_temp;
    D47_comp(2,1) = std_dep_temp;
    D47_comp(3,1) = 0;
    D47_comp(4,1) = long;
    D47_comp(5,1) = lati;
    D47_comp(6,1) = mean(mean(mean(mixdep,'omitnan'),'omitnan'),'omitnan');
end
if max_prod_temp + std_prod_temp <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && max_prod_temp + std_prod_temp >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,2) = max_prod_temp;
    D47_comp(2,2) = std_prod_temp;
    D47_comp(3,2) = 1;
    D47_comp(4,2) = long;
    D47_comp(5,2) = lati;
elseif max_prod_temp - std_prod_temp <= D47_temp(nrsamp) + D47_temp_SD(nrsamp) && max_prod_temp - std_prod_temp >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,2) = max_prod_temp;
    D47_comp(2,2) = std_prod_temp;
    D47_comp(3,2) = 1;
    D47_comp(4,2) = long;
    D47_comp(5,2) = lati;
else
    D47_comp(1,2) = max_prod_temp;
    D47_comp(2,2) = std_prod_temp;
    D47_comp(3,2) = 0;
    D47_comp(4,2) = long;
    D47_comp(5,2) = lati;
end
if max_prod_temp1 + std_prod_temp1 <= D47_temp(nrsamp) + D47_temp_SD(nrsamp)  && max_prod_temp1 + std_prod_temp1 >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,3) = max_prod_temp1;
    D47_comp(2,3) = std_prod_temp1;
    D47_comp(3,3) = 1;
    D47_comp(4,3) = long;
    D47_comp(5,3) = lati;
elseif max_prod_temp1 - std_prod_temp1 <= D47_temp(nrsamp) + D47_temp_SD(nrsamp)  && max_prod_temp1 - std_prod_temp1 >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,3) = max_prod_temp1;
    D47_comp(2,3) = std_prod_temp1;
    D47_comp(3,3) = 1;
    D47_comp(4,3) = long;
    D47_comp(5,3) = lati;
else
    D47_comp(1,3) = max_prod_temp1;
    D47_comp(2,3) = std_prod_temp1;
    D47_comp(3,3) = 0;
    D47_comp(4,3) = long;
    D47_comp(5,3) = lati;
end
if meantempsst + stdtempsst<= D47_temp(nrsamp) + D47_temp_SD(nrsamp)  && meantempsst + stdtempsst >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,4) = meantempsst;
    D47_comp(2,4) = stdtempsst;
    D47_comp(3,4) = 1;
    D47_comp(4,4) = long;
    D47_comp(5,4) = lati;
elseif meantempsst - stdtempsst <= D47_temp(nrsamp) + D47_temp_SD(nrsamp)  && meantempsst - stdtempsst >= D47_temp(nrsamp) - D47_temp_SD(nrsamp)
    D47_comp(1,4) = meantempsst;
    D47_comp(2,4) = stdtempsst;
    D47_comp(3,4) = 1;
    D47_comp(4,4) = long;
    D47_comp(5,4) = lati;
else
    D47_comp(1,4) = meantempsst;
    D47_comp(2,4) = stdtempsst;
    D47_comp(3,4) = 0;
    D47_comp(4,4) = long;
    D47_comp(5,4) = lati;
end
D47_comp(1,5) = D47_temp(nrsamp);
D47_comp(2,5) = D47_temp_SD(nrsamp);
D47_comp(3,5) = 1;
D47_comp(4,5) = mdepD47(1);
D47_comp(5,5) = mdepD47(2);

te = array2table(D47_comp,'VariableNames',{'Above mixed layer Temp (°C)','Maximum Production Temp (°C)', 'Photic Zone Temp (°C)', 'Sea Surface Temp (°C)', 'D47 Temp (°C)'});
%% For assembling and exporting into a matrix and table
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
writetable(D47_mix1,'ttemps_mixlay.xlsx');
writetable(D47_max1,'ttemps_maxpro.xlsx');
writetable(D47_phot1,'ttemps_phozo.xlsx');
writetable(D47_sst1,'ttemps_SST.xlsx');

end