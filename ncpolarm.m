function [maph]=ncpolarm(varargin)

% NCPOLARM: Square or rectangular polar stereographic map with land and grid
%
% function [maph]=ncpolarm('lat',lat,'lon',lon,'grid','label','asm','seaice')
%
% This function generates a polar stereographic map for either the 
% Arctic [default] or Antarctic on a square-framed map (the default
% stereographic map in matlab's mapping toolbox produces a circular plot,
% but this is not always the the most space-efficient method for 
% displaying oceanographic and meteorological information at the 
% poles). The are also several special Arctic preset maps available in 
% this function for displaying information about the Arctic System.
% 
% Several optional arguments may be used:
%
% Input:
% 'lat'   - Provides the maximum equatorward extent the map
%           and is proceeded by a scalar latitude.  Positive
%           latitudes produce an Arctic polar stereographic
%           map and negative latitudes produce an Antarctic 
%           polar stereographic map.  The default is 60 degrees
%           north, giving a default Arctic map.
%
% lat     - Actual latitude value.
%
% 'lon'   - Central meridian in the lower (upper) portion for the
%           northern (southern) hemisphere. The default is the
%           Greenwich Meridian.
%
% lon     - Actual value of lon.
%
% 'noland'- Stop land mask being plotted.
%
% 'grid'  - Overlays a grid on the map.
%
% 'label' - Adds latitude and longitude labels to the grid.
%
% 'asm'   - Arctic System map projection designed to include
%           areas that may be considered part of the regional
%           "Arctic System". This option does not work with 'seaice'. 
%
% 'seaice'- Plot most efficient rectangular map for Arctic
%           sea ice. This option does not work if 'asm' is selected.
%
% Output:
% This function generates a handle for the map:  maph
%
% An Azimuthal Stereographic Projection map projection is
% used.  Typing 'help stereo' in matlab will provide more
% information about this projection.  To see the distortion 
% provided by this map, use the tissot function.
%
% Examples:
% ncpolarm('asm') 
% ncpolarm('seaice')
% ncpolarm('lat',40,'lon',-40,'grid','label')
% ncpolarm('lat',-60,'grid')
%
% When printing maps with Postscript or EPSC, it is 
% advisable to use the -zbuffer option to avoid
% thin white filaments showing up as part of the land 
% filling. This is an artifact from within Matlab's
% Mapping toolbox.
%
% Written by Andrew Roberts
% Department of Oceanography, Naval Postgraduate School, 2011
% Tested using MATLAB Version 7.11.0.584 (R2010b)
%
% $Id: ncpolarm.m 256 2011-02-12 02:24:48Z aroberts $

% map handle
global maph

% Check that the mapping toolbox is installed
h=ver('map') ; 
if length(h)==0 ; error('Mapping toolbox not installed') ; end

% defaults
centralmeridian=0 ; 
equatorextent=60 ; 
parallelstop=80;
grid=0;
label=0;
land=1;
asm=0;
seaice=0;

% edge, grid and label color
gridcolor=[0.400 0.400 0.400] ; 

if nargin == 0
 disp('Default Arctic stereographic map');
else
 i=0;
 while i < nargin
  i=i+1;
  switch lower(varargin{i})
   case 'seaice'
	   seaice=1;
   case 'asm'
	   asm=1;
   case 'lat'
  	   if(i+1>nargin || ischar(varargin{i+1}))
		   error(['Missing latitude argument for ',varargin{i}])
  	   end
	   i=i+1;
	   equatorextent=varargin{i};
   case 'lon'
  	   if(i+1>nargin || ischar(varargin{i+1}))
		   error(['Missing longitiude argument for ',varargin{i}])
  	   end
	   i=i+1;
	   centralmeridian=varargin{i};
   case 'grid'
	   grid=1;
   case 'label'
	   label=1;
   case 'noland'
	   land=0;
   otherwise
	   error(['Option ',varargin{i},' is incorrect'])
  end
 end
end

if seaice & asm
 error('Cannot choose both seaice and asm as options in ncpolarm')
end

if asm % create preset arctic system map
 equatorextent=40;
 centralmeridian=-45;
elseif seaice % create preset arctic system map
 equatorextent=40;
 centralmeridian=-45;
end

% make polar map southern hemisphere if southern extent is zero
if(abs(equatorextent)>89)
	error('Latitude extent of map too close to a pole.');
elseif(equatorextent<0)
	mult=-1;
	disp('Southern Hemisphere')
else
	mult=1;
	disp('Northern Hemisphere')
end

% clear current axes
cla; 

% set default map properties
defaultm; 

% map projection
maph=axesm('MapProjection','stereo',...
  'AngleUnits','degrees',...
  'Aspect','normal',...
  'FalseNorthing',0,...
  'FalseEasting',0,...
  'MapLatLimit',sort([0 mult*90]),...
  'Geoid',[1 0],...
  'Origin',[mult*90 centralmeridian 0],...
  'Scalefactor',1,...
  'Frame','off',...
  'FFill',2000,...
  'FLatLimit',[-mult*Inf mult*90],...
  'FEdgeColor',gridcolor,...
  'FFaceColor','white',...
  'FLineWidth',0.25);

% add land to plot
if (land==1)
 land = shaperead('landareas.shp', 'UseGeoCoords', true);
 h=getm(maph);
 if strmatch(h.mapprojection,'globe') ;
  hl=linem([land.Lat],[land.Lon],'Color',0.45*[1 1 1]);
 else
  hl=geoshow(maph, land, 'FaceColor',[0.80 0.80 0.90],'EdgeColor',0.45*[1 1 1]);
 end
end

% map grid
if(grid==1) 
 setm(maph,'Grid','on',...
  'Galtitude',Inf,...
  'GColor',gridcolor,...
  'GLinestyle',':',...
  'Glinewidth',0.5,...
  'MLinefill',2000,...
  'MLineLimit',[mult*parallelstop 0],...
  'MLineException',[-90 0 90 180],...
  'MLineLocation',[0 30 60 90 120 150 180 -30 -60 -90 -120 -150],...
  'MLineVisible','on',...
  'PLineException',[],...
  'PLineFill',2000,...
  'PLineLimit',[-180 180],...
  'PLineLocation',[10],...
  'PLineVisible','on')
end



% grid label
if(label==1)
 if(mult>0)
  MLabelLocation=[0 180];
  PLabelMeridian=97;
 else
  MLabelLocation=[0 180];
  PLabelMeridian=135;
 end
 if asm;
  Mpos=37.5;
 elseif seaice;
  Mpos=46.5;
 else
  Mpos=equatorextent+mult*(90-mult*equatorextent)/10;
 end
 if abs(equatorextent)<40
  PLabelLocation=mult*[0 30 60];
 else
  PLabelLocation=mult*[40 50 60 70];
 end
 PLabelLocation=PLabelLocation(abs(PLabelLocation)>=abs(equatorextent));

 setm(maph,'Fontangle','normal',...
  'FontColor',gridcolor,...
  'Fontname','helvetica',...
  'FontSize',8,...
  'FontUnits','points',...
  'FontWeight','normal',...
  'LabelFormat','none',...
  'LabelRotation','off',...
  'MeridianLabel','on',...
  'MLabelLocation',MLabelLocation,...
  'MLabelParallel',Mpos,...
  'MLabelRound',0,...
  'ParallelLabel','on',...
  'PLabelLocation',PLabelLocation,...
  'PLabelMeridian',PLabelMeridian,...
  'PLabelRound',0);
end

% fit map tightly to axes
tightmap 
xlim=get(gca,'xlim');
ylim=get(gca,'ylim');
[lat,leftlon]=minvtran(xlim(1),ylim(1)+0.5*diff(ylim));
[lat,rightlon]=minvtran(xlim(2),ylim(1)+0.5*diff(ylim));
[lat,bottomlon]=minvtran(xlim(1)+0.5*diff(xlim),ylim(1));
[lat,toplon]=minvtran(xlim(1)+0.5*diff(xlim),ylim(2));
if asm
 [x1,y] = mfwdtran(equatorextent+7,leftlon);
 [x2,y] = mfwdtran(equatorextent+9,rightlon);
 [x,y1] = mfwdtran(equatorextent+3,bottomlon);
 [x,y2] = mfwdtran(equatorextent,toplon);
elseif seaice
 [x1,y] = mfwdtran(equatorextent+14,leftlon);
 [x2,y] = mfwdtran(equatorextent+16,rightlon);
 [x,y1] = mfwdtran(equatorextent+3,bottomlon);
 [x,y2] = mfwdtran(equatorextent,toplon);
else
 [x1,y] = mfwdtran(equatorextent,leftlon);
 [x2,y] = mfwdtran(equatorextent,rightlon);
 [x,y1] = mfwdtran(equatorextent,bottomlon);
 [x,y2] = mfwdtran(equatorextent,toplon);
end

set(gca,'Xlim',[x1 x2]); set(gca,'Ylim',[y1 y2]);

if nargout==0 ; clear maph; end

% Remove meridians and parallels from the legend (if one is used)
h1=handlem('parallel');
set(h1,'Clipping','on');
hCGroup=hggroup;
set(h1,'Parent',hCGroup)
set(get(get(hCGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2=handlem('meridian');
set(h2,'Clipping','on');
hSGroup=hggroup;
set(h2,'Parent',hSGroup)
set(get(get(hSGroup,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

