function plotNorthSeaBathymetricMap(pathOsparRegionsShapefile,...
    pathOsparSubregionsShapefile,pathOsparMpasShapefile,...
    pathUkOffshoreGcsLicensesShapefile,pathAreaStudyShapefile,...
    HPLC,filenameSedimentCoreData,filenameCefasCruiseData)

% PLOTNORTHSEABATHYMETRICMAP Plots a bathymetric map of the North Sea 
%   showing the locations on in situ data for this project as well as the
%   area for which satellite remote-sensing data have been extracted.
%
%   INPUT:
%       pathOsparRegionsShapefile          - shapefile with OSPAR regions
%       pathOsparSubregionsShapefile       - shapefile with OSPAR subregions
%       pathOsparMpasShapefile             - shapefile with OSPAR MPAs
%       pathUkOffshoreGcsLicensesShapefile - shapefile with NSTA's UKCS offshore licensed sites
%       pathAreaStudyShapefile             - shapefile with our area of study
%       HPLC                               - HPLC data filtered
%       filenameSedimentCoreData           - excel file with Malini's core samples locations
%       filenameCefasCruiseData            - excel file with cruise sample locations
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 7 May 2024
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Define regions and locations that will feature in the map

% Load the different geographic features and point locations for the map.
% Sources of the data:
%   OSPAR regions: https://odims.ospar.org/en/submissions/ospar_regions_2017_01/
%   OSPAR subregions: https://gis.ices.dk/sf/index.html?widget=station&program=OSPAR#
%   OSPAR MPAs: https://odims.ospar.org/en/maps/map-marine-protected-areas/
%   UK GCS sites: https://opendata-nstauthority.hub.arcgis.com/datasets/514393a0d3374847ab4251293781c6e6_0/explore
%   SmartBuoy network coordinates: https://data.cefas.co.uk/view/66
%   HPLC: https://data.cefas.co.uk/view/53

% OSPAR regions
osparReg = m_shaperead(pathOsparRegionsShapefile); % lat/lon coordinates
osparRegInfo = shapeinfo(pathOsparRegionsShapefile);
osparRegProj = osparRegInfo.CoordinateReferenceSystem;

% OSPAR subregions
osparSubreg = m_shaperead(pathOsparSubregionsShapefile); % map coordinates
osparSubregInfo = shapeinfo(pathOsparSubregionsShapefile);
osparSubregProj = osparSubregInfo.CoordinateReferenceSystem;

% OSPAR MPAs
osparMpas = m_shaperead(pathOsparMpasShapefile); % map coordinates
osparMpasInfo = shapeinfo(pathOsparMpasShapefile);
osparMpasProj = osparMpasInfo.CoordinateReferenceSystem;

% UK Continental Shelf Carbon Storage sites
ukcsGcsSites = m_shaperead(pathUkOffshoreGcsLicensesShapefile); % lat/lon coordinates

% Area of study around the Endurance
areaStudy = m_shaperead(pathAreaStudyShapefile); % lat/lon coordinates

% HPLC locations
hplc.lat = HPLC.Latitude;
hplc.lon = HPLC.Longitude;

% Selected buoys from the SmartBuoy network.
% Dowsing, North Dogger and Oyster Ground form a triangle around the
% Endurance site. They are no longer operational
buoys.lon = [1.056,  2.280,  4.039];
buoys.lat = [53.531, 55.683, 54.414];
buoys.name = {'Dowsing','North Dogger','Oyster Ground'};

% Malini core samples
TM = readtable(fullfile('.','data','raw','BGS',filenameSedimentCoreData));
cores.lat = TM.Y;
cores.lon = TM.X;

% Cruise
TC = readtable(fullfile('.','data','raw','CEFAS_Oxford_shipboard_survey',filenameCefasCruiseData),...
    'Sheet','Results','Range','A3:P9','ReadRowNames',true);
cruise.lat = TC.Latitude;
cruise.lon = TC.Longitude;

%% Figure

% Colour definitions (potential colours that I could use)
% Divide by 255 to convert to an RGB triplet between 0 and 1
mycolpalette.red            = [219 11 11]./255;
mycolpalette.pink           = [232 190 255]./255;
mycolpalette.beige          = [228 222 200]./255;
mycolpalette.blue           = [69 114 196]./255;
mycolpalette.purple         = [88 83 137]./255;
mycolpalette.banana         = [254 255 115]./255;
mycolpalette.sand           = [255 214 140]./255;
mycolpalette.orange         = [235 122 52]./255;
mycolpalette.green          = [0.7 1.0 0.7];
mycolpalette.darkgrey       = [0.4 0.4 0.4];
mycolpalette.lightgrey      = [0.7 0.7 0.7];
mycolpalette.deepblue       = [8 81 156]./255;
mycolpalette.salmon         = [255 168 127]./255;
mycolpalette.greencontinent = [200 222 175]./255;
mycolpalette.brown          = [168 56 0]./255;
mycolpalette.bordergrey     = [93 128 147]./255;

% Colour choices
colourContinent = mycolpalette.greencontinent;
colourCountryBorders = mycolpalette.darkgrey;
colourDeepBathymetry = mycolpalette.deepblue;
colourGcsSites = mycolpalette.banana;
colourEnduranceGcsSite = mycolpalette.brown;
colourOsparBorders = mycolpalette.bordergrey;

% Define the UK Shelf Sea boundary box
northsea.lon = [-10.5 8.5]; 
northsea.lat = [48.2 61.5];
northsea.lonticks = [-9 -6 -3 0 3 6];
northsea.latticks = [48 51 54 57 60];

% Define bathymetry layers/contours
bathyContours = [-200:50:-50 -30 -10];

% Options for labelling
isAddLabelsManually = 1; % 1=add labels manually using PPT text boxes, 0=use Matlab


figure()
hfig = set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.50],'Color','w');
m_proj('Equidistant','long',northsea.lon,'lat',northsea.lat); 

% Bathymetry
[CS,CH] = m_etopo2('contourf',bathyContours,'edgecolor','none');

% Coastline & borders
m_gshhs('i','patch',colourContinent,'EdgeColor','none'); % 'i' for intermediate resolution coastaline
m_gshhs('ib','color',colourCountryBorders) % 'b' for country borders

% In order to add colour to the bathymetric regions below the maximum depth 
% level defined by the bathymetric contours, we will add a background
% colour to the figure axes
hax = gca;
set(hax,'Color',colourDeepBathymetry)
hax.Position(1) = hax.Position(1) - 0.14; % move the axes to the left to make space for the inset figure

% Add grid labels (keep those two lines separated)
m_grid('xtick',northsea.lonticks,'ytick',northsea.latticks,'tickdir','out','yaxislocation','left','FontSize',11)
m_grid('xticklabels',[],'yticklabels',[],'ticklen',0,'linestyle','none')

% for iMpa = 1:length(osparMpas.ncst)
%     % Project x-y map coordinates to latitude/longitude coordinates
%     [lat,lon] = projinv(osparMpasProj,osparMpas.ncst{iMpa}(:,1),osparMpas.ncst{iMpa}(:,2));
%     if (iMpa == 1)
%         h7 = m_line(lon,lat,'linewi',1,'color','r'); 
%     else
%         m_line(lon,lat,'linewi',1,'color','r');
%     end
%     hold on
% end

% Add offshore GCS sites
for iGcsSite = 1:length(ukcsGcsSites.ncst)
    if (iGcsSite == 1) % Endurance is the first one!
        h1 = m_patch(ukcsGcsSites.ncst{iGcsSite}(:,1),ukcsGcsSites.ncst{iGcsSite}(:,2),...
            colourEnduranceGcsSite,'EdgeColor',colourCountryBorders);
    elseif (iGcsSite == 2)
        h2 = m_patch(ukcsGcsSites.ncst{iGcsSite}(:,1),ukcsGcsSites.ncst{iGcsSite}(:,2),...
            colourGcsSites,'EdgeColor',colourCountryBorders);
    elseif (iGcsSite > 2)
        m_patch(ukcsGcsSites.ncst{iGcsSite}(:,1),ukcsGcsSites.ncst{iGcsSite}(:,2),...
            colourGcsSites,'EdgeColor',colourCountryBorders);
    end
    hold on
end

% Add OSPAR subregions -select only offshore ones
for iSubreg = [1, 2, 4, 6:9, 14, 18, 20, 23, 26:28, 30, 33, 41, 42, 48:50] %1:length(osparSubreg.ncst) 
    % Project x-y map coordinates to latitude/longitude coordinates
    [lat,lon] = projinv(osparSubregProj,osparSubreg.ncst{iSubreg}(:,1),osparSubreg.ncst{iSubreg}(:,2));
    if (iSubreg == 1)
        h3 = m_line(lon,lat,'LineWidth',1,'Color',colourOsparBorders);
    else
        m_line(lon,lat,'LineWidth',1,'Color',colourOsparBorders); 
    end
    hold on
end

% % % % Add OSPAR regions
% % % for iReg = 1:length(osparReg.ncst) 
% % %     m_line(osparReg.ncst{iReg}(:,2),osparReg.ncst{iReg}(:,1),'LineWidth',1,'Color',colourOsparBorders);
% % %     hold on
% % % end

% Add HPLC locations
for iHplcPoint = 1:length(hplc.lat)
    if (iHplcPoint == 1)
        h4 = m_scatter(hplc.lon,hplc.lat,10,'r','filled');
    else 
        m_scatter(hplc.lon,hplc.lat,10,'r','filled');
    end
    hold on
end

% Add sediment core locations (blank)
% h5 = m_scatter(NaN,NaN,40,'k','x','LineWidth',1.5); hold on;

% Add SmartBuoy locations
buoyFeatures = {'Marker','o','Color','k','LineWidth',2,...
    'linest','none','MarkerFaceColor','w','clip','point'};
h6 = m_line(buoys.lon(1),buoys.lat(1),buoyFeatures{:}); hold on
h7 = m_line(buoys.lon(2),buoys.lat(2),buoyFeatures{:}); hold on
h8 = m_line(buoys.lon(3),buoys.lat(3),buoyFeatures{:}); hold on
% h6 = m_scatter(buoys.lon(1),buoys.lat(1),40,'k','filled'); hold on
% h7 = m_scatter(buoys.lon(2),buoys.lat(2),40,'k','filled'); hold on
% h8 = m_scatter(buoys.lon(3),buoys.lat(3),40,'k','filled'); hold on

% Add cruise locations (blank)
h9 = m_scatter(NaN,NaN,40,'k','x','LineWidth',1.5); hold on;

% Add Endurance boundary box
h10 = m_line(areaStudy.ncst{1}(:,1),areaStudy.ncst{1}(:,2),'LineWidth',2,'Color','k');

% Add bathymetry legend
colormap(m_colmap('blue')); 
caxis([-250 -10]); % set the limit higher than the deepest contour so that the darkest colour is not used in the deepest contour of choice
% [haxc,h] = m_contfbar([0.60 0.82],0.07,CS,CH,'levels',bathyContours,'endpiece','no','axfrac',0.025); % last number is thickness of bathymetry legend
[haxc,h] = m_contfbar(0.11,[0.0 0.45],CS,CH,'endpiece','no','axfrac',0.025,'FontSize',11); % last number is thickness of bathymetry legend
if (isAddLabelsManually == 0)
    title(haxc,{'Depth (m)'})
end

% Add text
if (isAddLabelsManually == 0)
    m_text(3.8,54,sprintf('%s\n%s',["Southern";"North Sea"]),... % two lines of text, centred
        'color','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top');
    m_text(0.2,56,sprintf('%s\n%s',["Northern";"North Sea"]),...
        'color','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top');
    m_text(2.2,55,sprintf('%s\n%s',["Dogger";"Bank"]),...
        'color','w','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top');
end
if (isAddLabelsManually == 0)
    m_text(-1,52.5,sprintf('%s',["GBR"]),... 
        'color','k','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top');
    m_text(3,52,sprintf('%s',["NLD"]),...
        'color','k','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top');
    m_text(3,52,sprintf('%s',["BEL"]),...
        'color','k','FontSize',12,'HorizontalAlignment','center','VerticalAlignment','top');
end
if (isAddLabelsManually == 0)
    for iBuoy = 1:4
        [~,ln,lt]=m_lldist([1 buoys.lon(iBuoy)],[1 buoys.lat(iBuoy)],40); 
        m_text(ln(end)-0.55,lt(end)-0.33,sprintf('%s',buoys.name{iBuoy}),'FontSize',12);
    end
end

% Add inset figure: create a new pair of axes inside current figure
axes('Position',[0.59 0.541 0.35 0.39])
m_proj('Equidistant','long',[-0.5 2.5],'lat',[53.4 55]); % longitudes have to come in 0.5 values to fill in the box
[CS,CH] = m_etopo2('contourf',bathyContours,'Edgecolor','none'); 
m_gshhs_i('patch',colourContinent,'edgecolor','none');
for iGcsSite = 1:length(ukcsGcsSites.ncst)
    if (iGcsSite == 1) % Endurance is the first one!
        m_patch(ukcsGcsSites.ncst{iGcsSite}(:,1),ukcsGcsSites.ncst{iGcsSite}(:,2),...
            colourEnduranceGcsSite,'EdgeColor',colourCountryBorders);
    elseif (iGcsSite > 1)
        m_patch(ukcsGcsSites.ncst{iGcsSite}(:,1),ukcsGcsSites.ncst{iGcsSite}(:,2),...
            colourGcsSites,'EdgeColor',colourCountryBorders);
    end
    hold on
end
for iSubreg = [1 27 49] % OSPAR subregions around Endurance
    [lat,lon] = projinv(osparSubregProj,osparSubreg.ncst{iSubreg}(:,1),osparSubreg.ncst{iSubreg}(:,2));
    m_line(lon,lat,'LineWidth',1,'color',colourOsparBorders); 
    hold on
end
% for iCore = 1:length(cores.lat)
%     m_scatter(cores.lon,cores.lat,60,'k','x','LineWidth',1.2); 
%     hold on
% end

% Add cruise locations
for iStat = 1:length(cruise.lat)
    m_scatter(cruise.lon,cruise.lat,40,'k','x','LineWidth',1.5); hold on;
end

% Find and add scatters inside my inset bounding box
lat_range = [53.4, 55]; 
lon_range = [-0.5, 2.5]; 
lat_indices = hplc.lat >= lat_range(1) & hplc.lat <= lat_range(2);
lon_indices = hplc.lon >= lon_range(1) & hplc.lon <= lon_range(2);
combined_indices = lat_indices & lon_indices;
lats_to_plot = hplc.lat(combined_indices);
lons_to_plot = hplc.lon(combined_indices);
for iHplcPoint = 1:length(lons_to_plot)
    m_scatter(lons_to_plot,lats_to_plot,20,'r','filled');
    hold on
end

m_line(areaStudy.ncst{1}(:,1),areaStudy.ncst{1}(:,2),'LineWidth',2,'color','k');
colormap(m_colmap('blue')); 
caxis([-250 -10]);
m_grid('xtick',([-2 -1 0 1 2]),'ytick',([53 54 55]),'tickdir','out','yaxislocation','right','FontSize',10)
m_grid('xticklabels',[],'yticklabels',[],'ticklen',0,'linestyle','none')
m_ruler([0.05 0.50],0.090,2,'tickdir','out','ticklen',[.007 .007]);
axis tight
box on

% Add legend
lg = legend([h2(1),h1(1),h10(1),h3(1),h6(1),h4(1),h9(1)],...
    'UK GCS licensed sites (spring 2024)',...
    'Endurance GCS site',...
    'Endurance study area',...
    'OSPAR subregions',...
    'CEFAS SmartBuoy network',...
    'CEFAS HPLC sampling sites',...
    'Oxford-CEFAS cruise locations');
lg.Position(1) = 0.64; lg.Position(2) = 0.12; 
lg.Orientation = 'vertical';
lg.FontSize = 11; 
lg.ItemTokenSize = [20,10];
set(lg,'Box','on','Color','w') % set background colour so that it's not painted blue

exportgraphics(gcf,fullfile('.','figures','bathym_map.png'),'Resolution',600)

end % plotNorthSeaBathymetricMap
