%% Code to convert CryoSat-2 summer radar freeboards to sea ice thickness

% Path to 'Developers Product'
folder = '...\Developers Product';

% Import files
cs2_files = dir(fullfile(folder,'*.nc'));

%% Plotting

% Download codes and relevant plotting packages to 'folder'
% panel-2.14 https://uk.mathworks.com/matlabcentral/fileexchange/20003-panel
% othercolor https://uk.mathworks.com/matlabcentral/fileexchange/30564-othercolor
% PyColormap4Matlab https://uk.mathworks.com/matlabcentral/fileexchange/68239-pycolormap4matlab
% ScientificColourMaps4 https://zenodo.org/record/2649252#.YjiRherP0uU
addpath(genpath(folder))

% Plot titles
Titles = {'May 1^{st}-15^{th}','May 16^{th}-31^{st}'...
,'June 1^{st}-15^{th}','June 16^{th}-30^{th}','July 1^{st}-15^{th}','July 16^{th}-31^{st}'...
,'August 1^{st}-15^{th}','August 16^{th}-31^{st}','September 1^{st}-15^{th}','September 16^{th}-30^{th}'};

% Remove plotting warnings
warning('off','all')

%% Cycle through gridded freeboard files

for i = 1:length(cs2_files)
    
    filename = fullfile(folder,cs2_files(i).name);
    
    % Load grid
    lat = ncread(filename,'Latitude');
    lon = ncread(filename,'Longitude');
    
    % Load radar freeboard
    rfb = ncread(filename,'Radar_Freeboard');
    rfb_err = ncread(filename,'Radar_Freeboard_Error'); % standard error
    N_leads = ncread(filename,'N_Leads'); % number of lead returns in grid cell
    
    %% Load and plot surface roughness from Lagrangian tracking
    
    % We obtain pan-Arctic sea ice surface roughness observations for
    % summer months by propagating CryoSat-2 estimates of ? from the 25-km
    % gridded Lognormal Altimeter Retracker Model (LARM) dataset forward
    % and backward from winter months, based on observations of the sea ice
    % drift. These roughness observations are assumed to represent the
    % standard deviation of the snow-sea ice interface.
    
    roughness_cs2 = ncread(filename,'Sea_Ice_Surface_Roughness');
    
    figure(100); clf(100);
    p = panel();
    p.pack(3,4);
    p.margin = 5;
    set(gcf,'Position',[51.1   56.2    987.2    711.2])
    for j = 1:10
        sigma = roughness_cs2(:,:,j);
        sigma(lat>88.2 | sigma>0.35) = NaN;
        [m,n] = ind2sub([4 3],j);
        p(n,m).select();
        ncpolarm('lat',72,'lon',0)
        pcolorm(lat,lon,sigma)
        caxis([0.05 0.35])
        colormap(flipud(othercolor('Spectral10',64)))
        textm(68,47,Titles{j},'fontsize',14,'fontweight','bold','HorizontalAlignment','right')
    end
    % Add colorbar
    p(n, m+1).select();
    cla
    set(gca, 'visible', 'off')
    caxis([0.05 0.35])
    colorbar('south','fontsize',14)
    text(0,0.36,'Sea Ice Surface Roughness [m]','fontsize',16)
    
    %% Load and plot melt pond fraction from Sentinel-3 OLCI
    
    % Remotely sensed observations of melt pond fraction are obtained from
    % the Sentinel-3 OLCI sensor through the University of Bremen
    % https://seaice.uni-bremen.de/melt-ponds/. This is a daily 12.5 km
    % pan-Arctic product based on the Version 1.5 algorithm of Istomina et
    % al. and covering the period between 2017 and 2020. Since cloud
    % cover can heavily obscure the coverage of daily observations and only
    % the final four years of our freeboard record had coinciding
    % measurements of f_p, we calculate a seasonal climatology of the f_p
    % observations that we could then apply to all years of our study
    % 2011-2020. Biweekly 80-km melt pond fraction fields are obtained from
    % the average of all cloud-free OLCI pixels between 2017 and 2020
    % within each two-week summer window and 80-km grid cell.
    
    mpf_clim_cs2 = ncread(filename,'Melt_Pond_Fraction');
    
    figure(101); clf(101);
    p = panel();
    p.pack(3,4);
    p.margin = 5;
    set(gcf,'Position',[51.1   56.2    987.2    711.2])
    for j = 1:10
        mpf = nanmedfilt2(mpf_clim_cs2(:,:,j),[3 3]);
        mpf(isnan(mpf_clim_cs2(:,:,j)) | lat>88.2 | mpf>0.7) = NaN;
        [m,n] = ind2sub([4 3],j);
        p(n,m).select();
        ncpolarm('lat',72,'lon',0)
        pcolorm(lat,lon,mpf)
        caxis([0 0.4])
        colormap(othercolor('Blues9',64))
        textm(68,47,Titles{j},'fontsize',14,'fontweight','bold','HorizontalAlignment','right')
    end
    % Add colorbar
    p(n, m+1).select();
    cla
    set(gca, 'visible', 'off')
    caxis([0 0.4])
    colorbar('south','fontsize',14)
    text(0.1,0.36,'Melt Pond Fraction','fontsize',16)

    %% Load and plot EM bias correction
    
    % The EM range bias correction is calculated from the function shown in
    % Figure S3 using estimates of ? from CryoSat-2 and f_p from Sentinel-3
    % OLCI. 
    
    rfb_bias_corr_cs2 = ncread(filename,'Radar_Freeboard_Bias_Correction');
    
    % This correction is not applicable when a significant snowpack is
    % present on the sea ice surface, so that melt pond coverage would be
    % limited. Therefore, we do not apply the correction when snow depth
    % h_s ? 60 cm (see below) and reduce the correction
    % linearly as a function of h_s between 0 and 60 cm.
    
    snow_depth = ncread(filename,'Snow_Depth');
    sd_redu = (1 - snow_depth/0.6); sd_redu(sd_redu<0) = 0;
    
    figure(102); clf(102);
    p = panel();
    p.pack(3,4);
    p.margin = 5;
    set(gcf,'Position',[51.1   56.2    987.2    711.2])
    for j = 1:10
        rfb_corr = -nanmedfilt2(rfb_bias_corr_cs2(:,:,j).*sd_redu(:,:,j),[3 3]);
        rfb_corr(isnan(rfb_bias_corr_cs2(:,:,j)) | lat>88.2 | rfb_corr>0.7) = NaN;
        [m,n] = ind2sub([4 3],j);
        p(n,m).select();
        ncpolarm('lat',72,'lon',0)
        pcolorm(lat,lon,rfb_corr)
        caxis([-0.2 0.2])
        colormap(othercolor('RdBu9',64))
        textm(68,47,Titles{j},'fontsize',14,'fontweight','bold','HorizontalAlignment','right')
    end
    % Add colorbar
    p(n, m+1).select();
    cla
    set(gca, 'visible', 'off')
    caxis([-0.2 0.2])
    colorbar('south','fontsize',14)
    text(0.02,0.36,'EM Bias Correction [m]','fontsize',16)

    %% Remove bias from radar freeboards and plot
    
    % Correction removed from rfb
    rfb_bias_corr_cs2(isnan(rfb_bias_corr_cs2)) = 0;
    rfb_corrected = rfb - rfb_bias_corr_cs2.*sd_redu;
    
    % Load OSISAF sea ice concentration
    sic = ncread(filename,'Sea_Ice_Concentration_OSISAF');
    
    figure(103); clf(103);
    p = panel();
    p.pack(3,4);
    p.margin = 5;
    set(gcf,'Position',[51.1   56.2    987.2    711.2])
    for j = 1:10
        rfbt = rfb_corrected(:,:,j);
        rfbt(isnan(rfbt) | lat>88.2 | rfbt>0.5) = NaN;
        [m,n] = ind2sub([4 3],j);
        p(n,m).select();
        ncpolarm('lat',72,'lon',0)
        contourfm(lat,lon,rfbt,0:0.01:0.3,'LineStyle','none')
        contourm(lat,lon,sic(:,:,j),[15 15],'linecolor',[0.8 0.8 0.8],'Fill','off','linewidth',0.1)
        caxis([0 0.3])
        colormap(flipud(othercolor('Spectral10',64)))
        textm(68,47,Titles{j},'fontsize',14,'fontweight','bold','HorizontalAlignment','right')
    end
    % Add colorbar
    p(n, m+1).select();
    cla
    set(gca, 'visible', 'off')
    caxis([0 0.3])
    colorbar('south','fontsize',14)
    text(0.06,0.36,'Radar Freeboard [m]','fontsize',16)
    
    %% Load and plot snow depth
    
    % Snow load (depth and density) estimates are obtained from the
    % Lagrangian snow evolution scheme SnowModel-LG. This scheme uses
    % the MERRA2 atmospheric reanalysis and NSIDC Polar Pathfinder ice
    % motion observations to simulate the accumulation of snow on Arctic
    % sea ice between September and April, while also modelling snowpack
    % metamorphism and melt between May and August.
    
    snow_depth = ncread(filename,'Snow_Depth');
    snow_density = ncread(filename,'Snow_Density');
    
    figure(104); clf(104);
    p = panel();
    p.pack(3,4);
    p.margin = 5;
    set(gcf,'Position',[51.1   56.2    987.2    711.2])
    for j = 1:10
        sd = snow_depth(:,:,j);
        sd(isnan(sd) | lat>88.2) = NaN;
        [m,n] = ind2sub([4 3],j);
        p(n,m).select();
        ncpolarm('lat',72,'lon',0)
        pcolorm(lat,lon,sd)
        contourm(lat,lon,sic(:,:,j),[15 15],'linecolor',[0.8 0.8 0.8],'Fill','off','linewidth',0.1)
        caxis([0 0.3])
        colormap(othercolor('Blues9',64))
        textm(68,47,Titles{j},'fontsize',14,'fontweight','bold','HorizontalAlignment','right')
    end
    % Add colorbar
    p(n, m+1).select();
    cla
    set(gca, 'visible', 'off')
    caxis([0 0.35])
    colorbar('south','fontsize',14)
    text(0.18,0.36,'Snow Depth [m]','fontsize',14)
    
    %% Convert bias-corrected radar freeboard to sea ice thickness
    
    % Load NSIDC sea ice type
    ice_type = ncread(filename,'Sea_Ice_Type'); % 0=FYI, 1=MYI
    
    % Sea ice density (fyi = 917, myi = 882 kg/m3)
    ocean_water_density = 1024;
    sea_ice_density = 916.7 - ice_type*34.7;
    
    % Normalized snow penetration depth (fraction of snow depth)
    D_pen = 0.9;
    
    % As a first approximation we assume the Ku-band radar penetrates a
    % constant 90% of the snow cover wherever snow is present between May
    % and September, which produces a largely consistent transition derived
    % sea ice thickness between April and May, and between September and
    % October. However, the assumed Ku-band radar penetration depth into
    % snow during the Arctic melting season does impact the estimated sea
    % ice thickness (see Supplementary Information Section B) and should
    % therefore be the subject of further study.
    
    % Sea ice freeboard from radar freeboard
    fb = rfb_corrected + D_pen*snow_depth.*(1./((1 + 0.51*snow_density*1e-3).^(-1.5)) - 1);
    
    % Convert to sea ice thickness
    sit = fb.*(ocean_water_density./(ocean_water_density-sea_ice_density)) + D_pen*snow_depth.*(snow_density./(ocean_water_density-sea_ice_density)) - (1 - D_pen)*snow_depth.*((ocean_water_density-snow_density)./(ocean_water_density - sea_ice_density));
    
    % Load sea ice thickness uncertainty
    sit_unc = ncread(filename,'Sea_Ice_Thickness_Uncertainty');
    
    % Plot SIT
    
    figure(105); clf(105);
    p = panel();
    p.pack(3,4);
    p.margin = 5;
    set(gcf,'Position',[51.1   56.2    987.2    711.2])
    for j = 1:10
        sitt = sit(:,:,j);
        sitt(isnan(sitt) | lat>88.2) = NaN;
        [m,n] = ind2sub([4 3],j);
        p(n,m).select();
        ncpolarm('lat',72,'lon',0)
        pcolorm(lat,lon,sitt)
        contourm(lat,lon,sic(:,:,j),[15 15],'linecolor',[0.8 0.8 0.8],'Fill','off','linewidth',0.1)
        caxis([0 4])
        colormap(flipud(othercolor('Spectral10',64)))
        textm(68,47,Titles{j},'fontsize',14,'fontweight','bold','HorizontalAlignment','right')
    end
    % Add colorbar
    p(n, m+1).select();
    cla
    set(gca, 'visible', 'off')
    caxis([0 4])
    colorbar('south','fontsize',14)
    text(0.1,0.36,'Sea Ice Thickness [m]','fontsize',14)
    
    % Plot SIT uncertainty
    
    figure(106); clf(106);
    p = panel();
    p.pack(3,4);
    p.margin = 5;
    set(gcf,'Position',[51.1   56.2    987.2    711.2])
    load('lajolla')
    for j = 1:10
        sitt = sit_unc(:,:,j);
        sitt(isnan(sitt) | lat>88.2) = NaN;
        [m,n] = ind2sub([4 3],j);
        p(n,m).select();
        ncpolarm('lat',72,'lon',0)
        pcolorm(lat,lon,sitt)
        contourm(lat,lon,sic(:,:,j),[15 15],'linecolor',[0.8 0.8 0.8],'Fill','off','linewidth',0.1)
        caxis([0 0.8])
        colormap(lajolla)
        textm(68,47,Titles{j},'fontsize',14,'fontweight','bold','HorizontalAlignment','right')
    end
    % Add colorbar
    p(n, m+1).select();
    cla
    set(gca, 'visible', 'off')
    caxis([0 0.8])
    colorbar('south','fontsize',14)
    text(-0.04,0.36,'Ice Thickness Uncertainty [m]','fontsize',14)
    
    %% Interpolate sea ice thickness in regions of missing data
    
    % Waveform parameters for CryoSat-2 samples classified as sea ice
    sig0_ice = ncread(filename,'Sigma_0_Ice');
    N_ice = ncread(filename,'N_Ice');
    PP_ice = ncread(filename,'Pulse_Peakiness_Ice');
    RP_ice = ncread(filename,'RIP_Peakiness_Ice');
    
    % Add thin ice in MIZ based on waveform parameters and SIC
    thin_ice = NaN(size(sit));
    for j = 1:size(sit,3)
        
        % Define constant thickness of unretrievable thin ice
        % as lowest 5th percentile thickness
        sitt = sit(:,:,j);
        ti = prctile(sitt(:),5);
        
        % Define MIZ
        miz = sic(:,:,j)>15 & sic(:,:,j)<60;
        
        % Valid specular returns (i.e. smooth level melting ice) 
        ids = sig0_ice(:,:,j)>44 | RP_ice(:,:,j)>25 | PP_ice(:,:,j)>0.3;
        
        % Valid grid cells for thin ice estimate
        idv = isnan(sit(:,:,j)) & miz & ids & ~isnan(snow_depth(:,:,j));
        
        thin_ice(:,:,j) = idv*ti;
        
    end
    thin_ice(thin_ice==0) = NaN;
    sit_int = nansum(cat(4,sit,thin_ice),4);
    sit_int(sit_int==0) = NaN;
    
    % Linear interpolation to no more than one grid cell from available
    % data
    for j = 1:size(sit,3)
        
        sitt = sit_int(:,:,j);
        
        % Using 2D median filter
        sit2 = nanmedfilt2(sitt,[3 3]);
        
        % Remove areas of valid data and data where SIC<15%
        sit2(~isnan(sitt)) = sitt(~isnan(sitt));
        sit2(sic(:,:,j)<=15 | isnan(sic(:,:,j))) = NaN;
        sit2(lat>88) = NaN;
        
        sit_int(:,:,j) = sit2;
        
    end
    
    % Plot SIT with interpolation
    
    figure(107); clf(107);
    p = panel();
    p.pack(3,4);
    p.margin = 5;
    set(gcf,'Position',[51.1   56.2    987.2    711.2])
    for j = 1:10
        sitt = sit_int(:,:,j);
        sitt(isnan(sitt) | lat>88.2) = NaN;
        [m,n] = ind2sub([4 3],j);
        p(n,m).select();
        ncpolarm('lat',72,'lon',0)
        pcolorm(lat,lon,sitt)
        contourm(lat,lon,sic(:,:,j),[15 15],'linecolor',[0.8 0.8 0.8],'Fill','off','linewidth',0.1)
        caxis([0 4])
        colormap(flipud(othercolor('Spectral10',64)))
        textm(68,47,Titles{j},'fontsize',14,'fontweight','bold','HorizontalAlignment','right')
    end
    % Add colorbar
    p(n, m+1).select();
    cla
    set(gca, 'visible', 'off')
    caxis([0 4])
    colorbar('south','fontsize',14)
    text(0.1,0.36,'Sea Ice Thickness [m]','fontsize',14)
    
end

