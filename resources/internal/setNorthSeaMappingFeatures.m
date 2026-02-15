function setNorthSeaMappingFeatures(dataTable,dataColName,lonColName,latColName,...
    cmin,cmax,pathEnduranceShapefile,pathAreaStudyShapefile,titleString,titleLocation)

% SETNORTHSEAMAPPINGFEATURES Plot data on a North Sea map.
%
%   INPUT: 
%       dataTable              - data table
%       dataColName            - name of the column containing the data
%       lonColName             - name of the column containing longitude values
%       latColName             - name of the column containing latitude values
%       cmin                   - minimum data value for the colour bar
%       cmax                   - maximum data value for the colour bar
%       pathEnduranceShapefile - shapefile with the Endurance GCS site
%       pathAreaStudyShapefile - shapefile with our area of study
%       titleString            - title string
%       titleLocation          - position of the title on the figure
% 
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 24 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Colour definitions (potential colours that I could use)

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
colourCountryBorders = mycolpalette.darkgrey;
colourEnduranceGcsSite = mycolpalette.brown;

%% Define the UK Shelf Sea boundary box

northsea.lon = [-10.5 8.5]; 
northsea.lat = [48.2 61.5];
northsea.lonticks = [-9 -6 -3 0 3 6];
northsea.latticks = [48 51 54 57 60];

%% Define bathymetry layers/contours

bathyContours = [0 -25 -50 -100 -250 -500 -1000];

%% Figure settings

% Define projection
m_proj('Equidistant','long',northsea.lon,'lat',northsea.lat); 
%m_proj('Lambert','long',northsea.lon,'lat',northsea.lat,'rectbox','on'); 

% Add bathymetry
m_elev('contour',bathyContours,'EdgeColor',[0.7 0.7 0.7])

% Add the data
[X,Y] = m_ll2xy(dataTable.(lonColName),dataTable.(latColName)); % transform coord (lon/lat) to map coord (x/y) to use with a function that is not m_
scatter(X,Y,20,dataTable.(dataColName),'filled')
shading flat; colormap(cmocean('haline'))
caxis([cmin cmax])

% Add the continents
m_gshhs('i','patch',[0.7 0.7 0.7],'EdgeColor','k'); % 'i' for intermediate resolution coastaline

% Add grid labels (keep those two lines separated)
m_grid('xtick',northsea.lonticks,'ytick',northsea.latticks,'tickdir','out','yaxislocation','left','FontSize',11)
m_grid('xticklabels',[],'yticklabels',[],'ticklen',0,'linestyle','none')

% Add Endurance site
enduranceSite = m_shaperead(pathEnduranceShapefile); % lat/lon coordinates
m_patch(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),...
    colourEnduranceGcsSite,'EdgeColor',colourCountryBorders);

% Add area of study
areaStudy = m_shaperead(pathAreaStudyShapefile); % lat/lon coordinates
m_line(areaStudy.ncst{1}(:,1),areaStudy.ncst{1}(:,2),'LineWidth',2,'color','k');

% Title
title(titleString)
axfigure = gca;
axfigure.TitleHorizontalAlignment = titleLocation;

set(gca,'FontSize',12)   

end % setNorthSeaMappingFeatures
