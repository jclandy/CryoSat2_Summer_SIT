%% Code to reproduce Figure 1 in the main paper

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

% Plot titles
Titles = {'January 1^{st}-15^{th}','January 16^{th}-31^{st}','February 1^{st}-14^{th}','February 15^{th}-28^{th}'...
    ,'March 1^{st}-15^{th}','March 16^{th}-31^{st}','April 1^{st}-15^{th}','April 16^{th}-30^{th}','May 1^{st}-15^{th}','May 16^{th}-31^{st}'...
    ,'June 1^{st}-15^{th}','June 16^{th}-30^{th}','July 1^{st}-15^{th}','July 16^{th}-31^{st}'...
    ,'August 1^{st}-15^{th}','August 16^{th}-31^{st}','September 1^{st}-15^{th}','September 16^{th}-30^{th}'...
    ,'October 1^{st}-15^{th}','October 16^{th}-31^{st}','November 1^{st}-15^{th}','November 16^{th}-30^{th}','December 1^{st}-15^{th}','December 16^{th}-31^{st}'};

% Remove plotting warnings
warning('off','all')

%% Choose year and import data

year = 2016;

% Load grid
lat = ncread(filename,'Latitude');
lon = ncread(filename,'Longitude');

% Load time vector [Y M D]
tvec = datevec(ncread(filename,'Time'));

% Load SIT
sit = ncread(filename,'Sea_Ice_Thickness');
sit_unc = ncread(filename,'Sea_Ice_Thickness_Uncertainty');

% Load SIC and ice type
sic = ncread(filename,'Sea_Ice_Concentration');
ice_type = ncread(filename,'Sea_Ice_Type');


%% Plot chosen year

figure(99); clf(99);
set(gcf,'Position',[-2559,-535,885.600000000000,1317.60000000000])
p = panel();
p.margin = [3 3 10 3];
p.de.margin = 3;
p.pack(6,4);

% Add subplots
idt = find(tvec(:,1)==year);
for i = 1:6
    for j = 1:4
        p(i, j).select();
        ncpolarm('lat',68,'lon',0);
        temp2 = double(sic(:,:,idt((i-1)*4 + j))>0.15);
        temp2(temp2==0 | lat>88) = NaN;
        geoshow(lat,lon,temp2*-0.01,'DisplayType','Surface')
        
        temp = sit(:,:,idt((i-1)*4 + j));
        temp2 = inpaint_nans(temp,1); temp(lat>88 & isnan(temp)) = temp2(lat>88 & isnan(temp));
        temp(lat>88.5 | temp==0 | sic(:,:,idt((i-1)*4 + j))<0.15) = NaN;
        contourfm(lat,lon,temp,0:0.1:4,'LineStyle','none')
        contourm(lat,lon,sic(:,:,idt((i-1)*4 + j)),[0.15 0.15],'k','Fill','off','linewidth',0.5)

        caxis([-0.01 4])
        cmp = flipud(othercolor('Spectral10',64));
        % cmp = getPyPlot_cMap('viridis',64); % gist_earth viridis 
        colormap([0.8 0.8 0.8; cmp])
        textm(61,48.5,Titles{(i-1)*4 + j},'fontsize',12,'fontweight','bold','HorizontalAlignment','right','backgroundcolor','white','margin',0.5)
    end
end
% Add colorbar
p(i, j).select();
caxis([0 4])
C = colorbar('west','fontsize',14,'fontweight','bold','ticks',0:1:4); pause(0.1);
set(C,'position',[0.968,0.016,0.012,0.142])
set(C,'YAxisLocation','right')


