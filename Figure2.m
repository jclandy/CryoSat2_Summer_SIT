%% Code to reproduce Figure 2 in the main paper

% Path to 'Sea Ice Thickness Product'
folder = '...\Sea Ice Thickness';

% Import files
cs2_files = dir(fullfile(folder,'*.nc'));
filename = fullfile(folder,cs2_files(1).name);

%% Plotting

% Download codes and relevant plotting packages to 'folder'
% panel-2.14 https://uk.mathworks.com/matlabcentral/fileexchange/20003-panel
% othercolor https://uk.mathworks.com/matlabcentral/fileexchange/30564-othercolor
% PyColormap4Matlab https://uk.mathworks.com/matlabcentral/fileexchange/68239-pycolormap4matlab
addpath(genpath(folder))

% Remove plotting warnings
warning('off','all')

%% CryoSat-2 data

% Load grid
lat_cs2 = ncread(filename,'Latitude');
lon_cs2 = ncread(filename,'Longitude');

% Load time vector [Y M D]
tvec_cs2 = datevec(ncread(filename,'Time'));
tvec_edges = [1 1; 1 16; 2 1; 2 15; 3 1; 3 16; 4 1; 4 16; 5 1; 5 16; 5 31; 6 15; 6 30; 7 16; 7 31; 8 16; 8 31; 9 15; 10 1; 10 16; 11 1; 11 15; 12 1; 12 16; 12 32];

% Load SIT
sit_cs2 = ncread(filename,'Sea_Ice_Thickness');
sit_unc_cs2 = ncread(filename,'Sea_Ice_Thickness_Uncertainty');

% Load SIC and ice type
sic_cs2 = ncread(filename,'Sea_Ice_Concentration');
ice_type_cs2 = ncread(filename,'Sea_Ice_Type');

%% Calculate sea ice volume from CryoSat-2 data

% LatLon to PolarStereo
[x_cs2,y_cs2]=polarstereo_fwd(lat_cs2,lon_cs2,6378137,0.08181919,70,0);

% Calculate SIV and uncertainty from each SIT and SIC grid
siv_cs2 = NaN(size(sit_cs2));
siv_unc_cs2 = NaN(size(sit_cs2));
for j = 1:size(sit_cs2,3)
    
    temp = sit_cs2(:,:,j);
    
    % Identify valid sea ice grid cells with missing data
    idg = find(~isnan(temp));
    idN = find(isnan(temp) & sic_cs2(:,:,j)>0.15);
    
    % Inverse distance weights
    [IDX,D] = knnsearch([x_cs2(idg) y_cs2(idg)],[x_cs2(idN) y_cs2(idN)],'K',25);
    w = (1./D)./sum(1./D,2);
    
    % Fill missing grid cells with inverse distance weighted interpolation
    temp_int = sum(temp(idg(IDX)).*w,2);
    temp(idN) = temp_int;
    temp(temp<0) = 0;
    
    % Sea ice volume from SIT*SIC*gridarea
    siv_cs2(:,:,j) = temp.*(sic_cs2(:,:,j))*(80e3^2); % 80 km grid cells
    
    % Uncertainty of missing grid cells from nearest neighbours
    temp = sit_unc_cs2(:,:,j);
    temp_int = sum(temp(idg(IDX)).*w,2);
    temp(idN) = temp_int;
    temp(temp<0) = 0;
    
    % Sea ice volume from SIT_uncertainty*SIC*gridarea
    siv_unc_cs2(:,:,j) = temp.*(sic_cs2(:,:,j))*(80e3^2);
        
end

%% Calculate pan-Arctic time series of CryoSat-2 sea ice volume

% Assume missing NSIDC ice type data in MIZ is first-year ice
ice_type2 = ice_type_cs2; ice_type2(isnan(ice_type2) & siv_cs2>0) = 0;

% Time series for total SIV
siv_total_cs2 = nansum(siv_cs2,[1 2]); siv_total_cs2 = siv_total_cs2(:);
siv_total_unc_cs2 = nansum(siv_unc_cs2,[1 2]); siv_total_unc_cs2 = siv_total_unc_cs2(:);

% Time series for FYI only
siv_FYI_cs2 = siv_cs2.*(1 - ice_type2);
siv_FYI_cs2 = nansum(siv_FYI_cs2,[1 2]); siv_FYI_cs2 = siv_FYI_cs2(:);
siv_FYI_unc_cs2 = siv_unc_cs2.*(1 - ice_type2);
siv_FYI_unc_cs2 = nansum(siv_FYI_unc_cs2,[1 2]); siv_FYI_unc_cs2 = siv_FYI_unc_cs2(:);

% Time series for MYI only
siv_MYI_cs2 = siv_cs2.*ice_type2;
siv_MYI_cs2 = nansum(siv_MYI_cs2,[1 2]); siv_MYI_cs2 = siv_MYI_cs2(:);
siv_MYI_unc_cs2 = siv_unc_cs2.*ice_type2;
siv_MYI_unc_cs2 = nansum(siv_MYI_unc_cs2,[1 2]); siv_MYI_unc_cs2 = siv_MYI_unc_cs2(:);

%% PIOMAS V2.1 data

% Download relevant datasets for period 2010-2020 to 'PIOMAS' directory in 'folder'
% (http://psc.apl.uw.edu/research/projects/arctic-sea-ice-volume-anomaly/data/model_grid)
% 'hiday' Sea ice thickness (Volume per unit Area) , daily mean
% https://pscfiles.apl.washington.edu/zhang/PIOMAS/data/v2.1/hiday/
% 'aiday' Sea ice concentration , daily mean
% https://pscfiles.apl.washington.edu/zhang/PIOMAS/data/v2.1/area/
addpath(genpath(folder))

% Grid parameters are provided in 'PIOMAS Grid' directory

% Import grid
fileid = fopen('Region.bin');
region_piomas = fread(fileid,[360 120]);
fclose(fileid);
fileid = fopen('Lat.bin');
lat_piomas = fread(fileid,[360 120],'single');
fclose(fileid);
fileid = fopen('Lon.bin');
lon_piomas = fread(fileid,[360 120],'single');
fclose(fileid);
fileid = fopen('GridArea.bin');
gridarea_piomas = fread(fileid,[360 120],'single'); % km^2
fclose(fileid);

% LatLon to PolarStereo
[x_piomas,y_piomas] = polarstereo_fwd(lat_piomas,lon_piomas,6378137,0.08181919,70,0);

%% Calculate PIOMAS SIC, SIT and SIV on same timescale as CS2 data

% PIOMAS SIT and SIC binary files
sitfiles = dir(fullfile(folder,'PIOMAS','Daily','hiday*'));
sicfiles = dir(fullfile(folder,'PIOMAS','Daily','aiday*'));

% Iterate through years averaging SIT and SIC on biweekly timescale
years = unique(tvec_cs2(:,1));
sit_piomas = NaN(size(lat_piomas,1),size(lat_piomas,2),size(sitfiles,1)*(size(tvec_edges,1)-1));
sic_piomas = NaN(size(lat_piomas,1),size(lat_piomas,2),size(sitfiles,1)*(size(tvec_edges,1)-1));
for i = 1:length(years)
    
    YD = 365;
    
    % Load daily sit (siv per unit area)
    fname = fullfile('Daily',sitfiles(i).name);
    fileID = fopen(fname);
    sit_full = fread(fileID,[360 120*YD],'single');
    fclose(fileID);
    
    sit_full = reshape(sit_full,[360 120 YD]);
    for j = 1:YD
        sitt = sit_full(:,:,j);
        sitt(region==15) = NaN; % Mask out land
        sit_full(:,:,j) = sitt;
    end
    
    % Load daily sic
    fname = fullfile('Daily',sicfiles(i).name);
    fileID = fopen(fname);
    sic_full = fread(fileID,[360 120*YD],'single');
    fclose(fileID);
    
    sic_full = reshape(sic_full,[360 120 YD]);
    for j = 1:YD
        sicc = sic_full(:,:,j);
        sicc(region==15) = NaN; % Mask out land
        sic_full(:,:,j) = sicc;
    end
        
    % Time
    piomas_days = datenum(ones(YD,1)*years(i),zeros(YD,1),(1:YD)');
    
    % Average over CS2 time intervals
    % Don't include zeros in thickness averaging to reproduce CS2 method of
    % filtering out non-ice returns
    sit_biweekly = NaN(size(sit_full,1),size(sit_full,2),size(tvec_edges,1)-1);
    sic_biweekly = NaN(size(sit_full,1),size(sit_full,2),size(tvec_edges,1)-1);
    for j = 1:size(tvec_edges,1)-1
        
        % All daily grids within a biweekly time interval
        idt = find(piomas_days>=datenum([years(i),tvec_edges(j,:)]) & piomas_days<datenum([years(i),tvec_edges(j+1,:)]));
        
        temp = sit_full(:,:,idt); temp(temp==0) = NaN;
        sit_mean = nanmean(temp,3);
        
        temp = sic_full(:,:,idt);
        sic_mean = nanmean(temp,3);
        
        sit_mean(isnan(sit_mean) & sic_mean>=0) = 0;
        
        sit_biweekly(:,:,j) = sit_mean;
        sic_biweekly(:,:,j) = sic_mean;
        
    end
        
    sit_piomas(:,:,(i-1)*(size(tvec_edges,1)-1)+1:i*(size(tvec_edges,1)-1)) = sit_biweekly;
    sic_piomas(:,:,(i-1)*(size(tvec_edges,1)-1)+1:i*(size(tvec_edges,1)-1)) = sic_biweekly;
end
    
% PIOMAS time vector
tvec_piomas = [repelem(years,size(tvec_edges,1)-1,1) repmat(unique(tvec_cs2(:,2:3),'rows'),size(years,1),1)];

% Calculate sea ice volume from PIOMAS
% Multiply SIT (i.e., SIV per unit area) by the grid cell area
siv_piomas = (sit_piomas*1e-3).*gridarea_piomas;

% Remove areas with SIC<15% to mirror the CS2 method
siv_piomas(sic_piomas<0.15) = NaN;

%% Prepare PIOMAS data for comapring to CS2

tvec_id = ismember(tvec_piomas,tvec_cs2(:,1:3),'rows');

% Mask areas not covered by CS2
k = boundary(x_cs2(:),y_cs2(:));
in = inpolygon(x_piomas(:),y_piomas(:),x_cs2(k),y_cs2(k));
cs2_mask = reshape(in, size(x_piomas));

% Resample NSIDC sea ice type data onto PIOMAS grid
ice_type2_piomas = NaN(size(sic_piomas(:,:,tvec_id)));
for j = 1:size(ice_type2,3)
    
    temp = ice_type2(:,:,j);
    [IDX,D] = knnsearch([x_cs2(:) y_cs2(:)],[x_piomas(:) y_piomas(:)],'K',1);
    ice_type2_piomas(:,:,j) = reshape(temp(IDX), size(x_piomas));    
    
end

%% Calculate pan-Arctic time series of PIOMAS sea ice volume

% Time series for total SIV
siv_all_piomas = siv_piomas(:,:,tvec_id); siv_all_piomas(cs2_mask.*ones(size(siv_all_piomas))<1) = NaN;
siv_total_piomas = nansum(siv_all_piomas,[1 2]); siv_total_piomas = siv_total_piomas(:);

% Time series for FYI only
siv_FYI_piomas = siv_all_piomas.*(1 - ice_type2_piomas);
siv_FYI_piomas = nansum(siv_FYI_piomas,[1 2]); siv_FYI_piomas = siv_FYI_piomas(:);

% Time series for MYI only
siv_MYI_piomas = siv_all_piomas.*ice_type2_piomas;
siv_MYI_piomas = nansum(siv_MYI_piomas,[1 2]); siv_MYI_piomas = siv_MYI_piomas(:);

%% Calculate time series of anomalies after removing the climatological seasonal cycles

% Climatology time vector
UU = unique(tvec(:,2:3),'rows');
[~,LocB] = ismember(tvec(:,2:3),UU,'rows');

% Seasonal variations from PIOMAS
PIOMAS_total_climatology = accumarray(LocB,siv_total_piomas*1e-3,size(UU(:,1)),@nanmean);
PIOMAS_total_SS = siv_total_piomas*1e-3 - PIOMAS_total_climatology(LocB);
PIOMAS_FYI_climatology = accumarray(LocB,siv_FYI_piomas*1e-3,size(UU(:,1)),@nanmean);
PIOMAS_FYI_SS = siv_FYI_piomas*1e-3 - PIOMAS_FYI_climatology(LocB);
PIOMAS_MYI_climatology = accumarray(LocB,siv_MYI_piomas*1e-3,size(UU(:,1)),@nanmean);
PIOMAS_MYI_SS = siv_MYI_piomas*1e-3 - PIOMAS_MYI_climatology(LocB);

% Seasonal variations from CryoSat-2
CS2_total_climatology = accumarray(LocB,siv_total_cs2*1e-9*1e-3,size(UU(:,1)),@nanmean);
CS2_total_SS = siv_total_cs2*1e-9*1e-3 - CS2_total_climatology(LocB);
CS2_FYI_climatology = accumarray(LocB,siv_FYI_cs2*1e-9*1e-3,size(UU(:,1)),@nanmean);
CS2_FYI_SS = siv_FYI_cs2*1e-9*1e-3 - CS2_FYI_climatology(LocB);
CS2_MYI_climatology = accumarray(LocB,siv_MYI_cs2*1e-9*1e-3,size(UU(:,1)),@nanmean);
CS2_MYI_SS = siv_MYI_cs2*1e-9*1e-3 - CS2_MYI_climatology(LocB);

% Anomaly correlation coefficients
corrcoef(CS2_total_SS,PIOMAS_total_SS)
corrcoef(CS2_FYI_SS,PIOMAS_FYI_SS)
corrcoef(CS2_MYI_SS,PIOMAS_MYI_SS)

%% Final time series plot

figure(99); clf(99);
set(gcf,'Position',[-2540.60000000000,-323.800000000000,1643.20000000000,1181.60000000000])
p = panel();
p.margin = [15 20 5 5];
p.pack({2.8/5 []},1);
p(2,1).pack(1,3);
p.de.margin = 20;

cmap = getPyPlot_cMap('PRGn',2);
winter = tvec(:,2)<=4 | tvec(:,2)>=10;
summer = tvec(:,2)>=5 & tvec(:,2)<=9;

p(1, 1).select();
cla
hold on
% Scale CS2 and PIOMAS ice volume to the same units
patch([datenum(tvec_cs2); flipud(datenum(tvec_cs2))],[siv_total_cs2 + siv_total_unc_cs2; flipud(siv_total_cs2 - siv_total_unc_cs2)]*1e-9*1e-3,[0.5 0.5 0.5],'facealpha',0.4,'edgecolor','none')
plot(datenum(tvec_cs2),siv_total_cs2*1e-9*1e-3,'k','linewidth',2)
patch([datenum(tvec_cs2); flipud(datenum(tvec_cs2))],[siv_FYI_cs2 + siv_FYI_unc_cs2; flipud(siv_FYI_cs2 - siv_FYI_unc_cs2)]*1e-9*1e-3,cmap(1,:),'facealpha',0.3,'edgecolor','none')
plot(datenum(tvec_cs2),siv_FYI_cs2*1e-9*1e-3,'color',cmap(1,:),'linewidth',1.5)
patch([datenum(tvec_cs2); flipud(datenum(tvec_cs2))],[siv_MYI_cs2 + siv_MYI_unc_cs2; flipud(siv_MYI_cs2 - siv_MYI_unc_cs2)]*1e-9*1e-3,cmap(2,:),'facealpha',0.3,'edgecolor','none')
plot(datenum(tvec_cs2),siv_MYI_cs2*1e-9*1e-3,'color',cmap(2,:),'linewidth',1.5)
plot(datenum(tvec_cs2),siv_total_piomas*1e-3,'k--','linewidth',2)
plot(datenum(tvec_cs2),siv_FYI_piomas*1e-3,'--','color',cmap(1,:),'linewidth',1.5)
plot(datenum(tvec_cs2),siv_MYI_piomas*1e-3,'--','color',cmap(2,:),'linewidth',1.5)
grid on
box on
datetick('x','yyyy')
xlim([734050 738157])
xlabel('Year')
ylabel('Sea Ice Volume [1000 km^3]')
legend('CryoSat-2 Total','','CryoSat-2 FYI','','CryoSat-2 MYI','','PIOMAS Total','PIOMAS FYI','PIOMAS MYI','Location','northwest')
text(datenum('09/15/2020'),29,'(a)','fontsize',20)
set(gca,'fontsize',14)
hold off  

p(2,1,1,1).select();
cla
hold on
scatter(CS2_total_SS(winter),PIOMAS_total_SS(winter),30,'k')
scatter(CS2_total_SS(summer),PIOMAS_total_SS(summer),30,'k','filled')
cfit = fit(CS2_total_SS,PIOMAS_total_SS,'poly1');
xx = min(CS2_total_SS):0.1:max(CS2_total_SS);
plot(xx,cfit(xx),'k','linewidth',2)
plot([-3 3],[-3 3],'k:','linewidth',1)
grid on
box on
xlim([-3 3])
ylim([-3 3])
xlabel('CryoSat-2 Sea Ice Volume Anomaly [1000 km^3]')
ylabel('PIOMAS Sea Ice Volume Anomaly [1000 km^3]')
legend('Total Winter','Total Summer','Best Fit','One-to-One Line','Location','northwest')
text(2.5,2.8,'(b)','fontsize',20)
set(gca,'fontsize',14)

p(2,1,1,2).select();
cla
hold on
scatter(CS2_FYI_SS(winter),PIOMAS_FYI_SS(winter),30,cmap(1,:))
scatter(CS2_FYI_SS(summer),PIOMAS_FYI_SS(summer),30,cmap(1,:),'filled')
cfit = fit(CS2_FYI_SS,PIOMAS_FYI_SS,'poly1');
xx = min(CS2_FYI_SS):0.1:max(CS2_FYI_SS);
plot(xx,cfit(xx),'color',cmap(1,:),'linewidth',2)
plot([-3 3],[-3 3],'k:','linewidth',1)
grid on
box on
xlim([-3 3])
ylim([-3 3])
xlabel('CryoSat-2 Sea Ice Volume Anomaly [1000 km^3]')
legend('FYI Winter','FYI Summer','Best Fit','One-to-One Line','Location','northwest')
text(2.5,2.8,'(c)','fontsize',20)
set(gca,'fontsize',14)

p(2,1,1,3).select();
cla
hold on
scatter(CS2_MYI_SS(winter),PIOMAS_MYI_SS(winter),30,cmap(2,:))
scatter(CS2_MYI_SS(summer),PIOMAS_MYI_SS(summer),30,cmap(2,:),'filled')
cfit = fit(CS2_MYI_SS,PIOMAS_MYI_SS,'poly1');
xx = min(CS2_MYI_SS):0.1:max(CS2_MYI_SS);
plot(xx,cfit(xx),'color',cmap(2,:),'linewidth',2)
plot([-3 3],[-3 3],'k:','linewidth',1)
grid on
box on
xlim([-3 3])
ylim([-3 3])
xlabel('CryoSat-2 Sea Ice Volume Anomaly [1000 km^3]')
legend('MYI Winter','MYI Summer','Best Fit','One-to-One Line','Location','northwest')
text(2.5,2.8,'(d)','fontsize',20)
set(gca,'fontsize',14)
