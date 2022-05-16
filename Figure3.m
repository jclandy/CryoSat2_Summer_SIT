%% Code to reproduce Figure 3 in the main paper

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

%% MASIE Arctic region mask

load('MASIE_region_mask.mat')

% Region Number IdentifierName
% 1	Beaufort Sea
% 2	Chukchi Sea
% 3	East Siberian Sea
% 4	Laptev Sea
% 5	Kara Sea
% 6	Barents Sea
% 7	Greenland Sea
% 8	Baffin Bay/Gulf of St. Lawrence
% 9	Canadian Archipelago
% 10	Hudson Bay
% 11	Central Arctic
% 12	Bering Sea
% 13	Baltic Sea
% 14	Sea of Okhotsk
% 15	Yellow Sea
% 16	Cook Inlet

%% Calculate pan-Arctic time series of biweekly CryoSat-2 sea ice volume and OSISAF sea ice extent

% Time series for total SIV
siv_sum = nansum(siv_cs2,[1 2]); siv_sum = siv_sum(:);
sivt = reshape(siv_cs2, [size(siv_cs2,1)*size(siv_cs2,2) size(siv_cs2,3)]);

% Time series for SIV in the Central Arctic
siv_sum_CA = nansum(sivt(region_mask_cs2(:)==11,:),1);

% Time series for total SIE
sie_sum = nansum(sic_cs2>0.15,[1 2])*80^2; sie_sum = sie_sum(:); % km^2

%% Download daily OSISAF sea ice concentration data

% Download OSI-SAF Daily Global Sea Ice Concentration climate data record
% (SMMR/SSMI/SSMIS) OSI-450 for 2010-2020 to 'SIC' directory in 'folder'
% https://osi-saf.eumetsat.int/products/osi-450

sic_files = dir(fullfile(folder,'SIC','*/*.nc'));

% SIC grid
lat_sic = double(ncread(fullfile(sic_files(1).folder,sic_files(1).name),'lat'));
lon_sic = double(ncread(fullfile(sic_files(1).folder,sic_files(1).name),'lon'));

% LatLon to PolarStereo
[x_sic,y_sic]=polarstereo_fwd(lat_sic,lon_sic,6378137,0.08181919,70,0);

% SIC time vector
tvec_sic = NaN(length(sic_files),3);
for j = 1:length(sic_files)
    if str2double(sic_files(j).name(32:35))<=2015
        tvec_sic(j,:) = [str2double(sic_files(j).name(32:35)) str2double(sic_files(j).name(36:37)) str2double(sic_files(j).name(38:39))];
    else
        tvec_sic(j,:) = [str2double(sic_files(j).name(33:36)) str2double(sic_files(j).name(37:38)) str2double(sic_files(j).name(39:40))];        
    end
end
sic_files = sic_files(tvec_sic(:,1)>=2010 & tvec_sic(:,1)<=2020);
tvec_sic = tvec_sic(tvec_sic(:,1)>=2010 & tvec_sic(:,1)<=2020,:);
t_sic = datenum(tvec_sic);

%% Resample to CryoSat-2 grid for a fair comparison

% Grid nearest neighbours
[IDX,D] = knnsearch([x_sic(:) y_sic(:)],[x_cs2(:) y_cs2(:)],'K',9);

% Iterate through days and resample grid
sic_full_daily = NaN(size(x_cs2,1),size(x_cs2,2),size(t_sic,1));
for j = 1:size(t_sic,1)
    
    sic_temp = double(ncread(fullfile(sic_files(j).folder,sic_files(j).name),'ice_conc'));
    
    % Resample from 25 km to 80 km CS2 grid
    D2 = D; D2(isnan(sic_temp(IDX))) = NaN;
    sic_temp = nanmean(sic_temp(IDX),2);
    sic_temp(nanmin(D2,[],2)>20e3) = NaN;
    sic_full_daily(:,:,j) = reshape(sic_temp, size(lat_cs2));

end

%% Calculate pan-Arctic time series of daily sea ice extent

% Time series for total daily SIE
sie_sum_daily = nansum(sic_full_daily>15,[1 2])*80^2; sie_sum_daily = sie_sum_daily(:); % km^2

%% Calculate lagged correlation matrices between ice extent and volume

% Identify where daily SIC observation intersects with the mid-point of a
% biweekly time interval
[~,ia] = ismember(tvec_cs2,tvec_sic,'rows');

% Day of year
doy = tvec_sic(tvec_sic==2010,:); doy(:,1) = [];

% Last viable day of daily SIC time series
idj_max = find(tvec_sic(:,1)==tvec_cs2(end,1) & tvec_sic(:,2)==tvec_cs2(end,2),1,'last') + 1;

% Correlations between SIV and SIE time series separated by up to 364 days
SIE_SIV_corr = NaN(365,365);
SIE_SIE_corr = NaN(365,365);
SIE_SIV_p = NaN(365,365);
SIE_SIE_p = NaN(365,365);
for i = 1:size(doy,1)
    
    % Find daily SIE in time series for i'th day of the year
    idx = find(tvec_sic(:,2)==doy(i,1) & tvec_sic(:,3)==doy(i,2));
    
    for j = 0:364
        
        % Lag time in days
        i_lag = i - j; i_lag = IF(i_lag>0, i_lag, size(doy,1)+i_lag);
        
        % All viable years at specified lag
        IN = (datenum(tvec_sic(idx,:))-j) >= datenum(tvec_cs2(1,:)) & idx<=(idj_max+j);
        
        % Indices of daily SIE that line up for specified lag
        idj = find(tvec_sic(:,2)==doy(i_lag,1) & tvec_sic(:,3)==doy(i_lag,2));
        idj = idj(datenum(tvec_sic(idj,:)) >= datenum(tvec_cs2(1,:)) & idj<=max(idx) & idj<=idj_max);
        
        [~,ib] = ismember(idj,ia);
        
        % If pairs of SIV and SIE are available at lag time
        if sum(ib)>0
            
            % SIE on target day
            SIE_temp = sie_sum_daily(idx);
            
            % SIE/SIV at lead time
            SIE_lags = sie_sum(ib);
            SIV_lags = siv_sum(ib) - siv_sum_CA(ib)';
            
            % SIV-SIE correlations (with or without detrending)
            lfit = fitlm(SIV_lags,SIE_temp(IN));
            % lfit = fitlm(detrend(SIV_lags),detrend(SIE_temp(IN)));
            SIE_SIV_corr(j+1,i) = sqrt(lfit.Rsquared.Ordinary);
            SIE_SIV_p(j+1,i) = coefTest(lfit);
            
            % SIE-SIE correlations (with or without detrending)
            lfit = fitlm(SIE_lags,SIE_temp(IN));
            % lfit = fitlm(detrend(SIE_lags),detrend(SIE_temp(IN)));
            SIE_SIE_corr(j+1,i) = sqrt(lfit.Rsquared.Ordinary);
            SIE_SIE_p(j+1,i) = coefTest(lfit);
            
        end
                 
    end   

end

%% Smooth correlation matrices to fill gaps

% Filtered matrices, July to Jan
idt = [152:365, 1:31];
tt = datenum(tvec_sic(tvec_sic==2010,:));
tt = tt(idt(1)) : (tt(idt(1)) + length(idt) - 1);
[x,y] = meshgrid(tt,0:364);

% Median filter with radius of 21 days
SIE_SIV_corr2 = nanmedfilt2(SIE_SIV_corr(:,idt),[41 41]);
SIE_SIE_corr2 = nanmedfilt2(SIE_SIE_corr(:,idt),[41 41]);
SIE_SIV_p2 = nanmedfilt2(SIE_SIV_p(:,idt),[41 41]);
SIE_SIE_p2 = nanmedfilt2(SIE_SIE_p(:,idt),[41 41]);

%% Plot final correlation matrices

figure(200); clf(200);
set(gcf,'Position',[-2470.20000000000,-35.8000000000000,2232.80000000000,777.600000000000])
p = panel();
p.margin = [25 20 15 15];
p.pack(1,3);
p.de.margin = 35;

cmap = flipud(othercolor('RdBu11',20));
cmap(9,:)=cmap(10,:); cmap(12,:)=cmap(11,:); % all corrs <0.2 in white
monthlab = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% SIV correlated with future SIE
p(1, 1).select();
cla
hold on
contourf(tt,0:364,SIE_SIV_corr2,-1:0.1:1,'linestyle','none')
contour(tt,0:364,tt-(0:364)'-tt(1)+1,-365:365/12:365,'linecolor',[0.5 0.5 0.5],'linewidth',0.25);
caxis([-1 1])
colormap(cmap)
datetick('x','mmm')
contour(tt,0:364,SIE_SIV_p2,[0.1 0.1],'linecolor','k','linewidth',1)
SIE_SIV_p2_PIOMAS = load(fullfile(folder,'PIOMAS Grid','SIE_SIV_p2_PIOMAS.mat')); % Load equiavlent significance matrix derived from PIOMAS data
contour(tt,0:364,SIE_SIV_p2_PIOMAS.SIE_SIV_p2,[0.1 0.1],'k:','linewidth',1) % PIOMAS p=0.1
stipple(x,y,SIE_SIV_corr2>=SIE_SIE_corr2 & abs(SIE_SIV_p2)<=0.1,'color','k','markersize',5)
xlim([min(tt) max(tt)+1])
yyaxis left
set(gca,'ycolor',[0.5 0.5 0.5])
ylim([0 365])
yticks(0:365/12:365)
yticklabels([{''} fliplr(monthlab(1:5)) fliplr(monthlab(6:12))])
ytickangle(35)
ylabel('Lead Time','color','k')
yyaxis right
set(gca,'ycolor',[0.5 0.5 0.5])
ylim([0 365])
yticks(1:365/12:364)
yticklabels([{''} monthlab(1) fliplr(monthlab(3:12))])
ytickangle(35)
hold off
set(gca,'fontsize',16)
xlabel('Target Day')
title('(a) Sea Ice Volume Correlated with Future Ice Extent')
box on

% SIE correlated with future SIE
p(1, 2).select();
cla
hold on
contourf(tt,0:364,SIE_SIE_corr2,-1:0.1:1,'linestyle','none')
contour(tt,0:364,tt-(0:364)'-tt(1)+1,-365:365/12:365,'linecolor',[0.5 0.5 0.5],'linewidth',0.25);
caxis([-1 1])
colormap(cmap)
datetick('x','mmm')
contour(tt,0:364,SIE_SIE_p2,[0.1 0.1],'linecolor','k','linewidth',1)
stipple(x,y,SIE_SIE_corr2>SIE_SIV_corr2 & abs(SIE_SIE_p2)<=0.1,'color','k','markersize',5)
xlim([min(tt) max(tt)+1])
yyaxis left
set(gca,'ycolor',[0.5 0.5 0.5])
ylim([0 365])
yticks(0:365/12:365)
yticklabels([{''} fliplr(monthlab(1:5)) fliplr(monthlab(6:12))])
ytickangle(35)
yyaxis right
set(gca,'ycolor',[0.5 0.5 0.5])
ylim([0 365])
yticks(1:365/12:364)
yticklabels([{''} monthlab(1) fliplr(monthlab(3:12))])
ytickangle(35)
hold off
set(gca,'fontsize',16)
xlabel('Target Day')
title('(b) Sea Ice Extent Correlated with Future Ice Extent')
box on
colorbar('north','fontsize',16,'ticks',-1:0.2:1)

% Mean September correlations
cmap = getPyPlot_cMap('PRGn',10);
p(1, 3).select();
cla
hold on
patch([0:364 fliplr(0:364)],[nanmean(SIE_SIV_corr2(:,doy(idt,1)==9),2) + nanstd(SIE_SIV_corr2(:,doy(idt,1)==9),[],2); flipud(nanmean(SIE_SIV_corr2(:,doy(idt,1)==9),2) - nanstd(SIE_SIV_corr2(:,doy(idt,1)==9),[],2))],cmap(2,:),'FaceAlpha',0.2,'edgecolor','none')
P = plot(0:364,nanmean(SIE_SIV_corr2(:,doy(idt,1)==9),2),'color',cmap(2,:),'linewidth',2);
plot([108 108],[0 1],':','color',cmap(2,:),'linewidth',1.5) % 108 days is the mean contour for p=0.1
patch([0:364 fliplr(0:364)],[nanmean(SIE_SIE_corr2(:,doy(idt,1)==9),2) + nanstd(SIE_SIE_corr2(:,doy(idt,1)==9),[],2); flipud(nanmean(SIE_SIE_corr2(:,doy(idt,1)==9),2) - nanstd(SIE_SIE_corr2(:,doy(idt,1)==9),[],2))],cmap(9,:),'FaceAlpha',0.2,'edgecolor','none')
Q = plot(0:364,nanmean(SIE_SIE_corr2(:,doy(idt,1)==9),2),'color',cmap(9,:),'linewidth',2);
plot([71 71],[0 1],':','color',cmap(9,:),'linewidth',1.5) % 71 days is the mean contour for p=0.1
T = patch([26 26 137 137],[0 1 1 0],'k','FaceAlpha',0.1,'edgecolor','none');
U = patch([249 249 364 364],[0 1 1 0],'r','FaceAlpha',0.1,'edgecolor','none');
ylabel('Correlation Coefficient')
xlabel('Lead Time')
xlim([0 365])
xticks(0:365/12:365)
xticklabels([fliplr(monthlab(1:9)) fliplr(monthlab(9:12))])
set(gca,'fontsize',16)
legend([P,Q,T,U],'SIV Correlated with Future SIE','SIE Correlated with Future SIE','SIV Anomaly Persistence Region','SIV Anomaly Reemergence Region')
title('(c) Mean September Correlation')
box on
grid on


