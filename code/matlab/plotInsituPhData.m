function plotInsituPhData(pathEnduranceShapefile,pathAreaStudyShapefile,PH)

% PLOTINSITUPHDATA Plot the in situ pH data.
%
%   INPUT:
%       pathEnduranceShapefile - shapefile with the Endurance GCS site
%       pathAreaStudyShapefile - shapefile with our area of study
%       PH                     - PH data filtered
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 25 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Definitions

cmin = 7.8; %min(P.pH);
cmax = 8.2; %max(P.pH);
c = [7.8 8 8.2];

%% Map locations

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.40],'Color','w') 

setNorthSeaMappingFeatures(PH,'pH','Longitude_degE','Latitude_degN',...
    cmin,cmax,pathEnduranceShapefile,pathAreaStudyShapefile,...
    sprintf('pH from SSB and UKOA programs, \\itN\\rm \\bf= %d',height(PH)),'center')

cb = colorbar('Location','southoutside');
cb.Label.String = 'pH';
cb.YTick = c;
cb.YTickLabel = c;

exportgraphics(gcf,fullfile('.','figures','pH_map_all.png'),'Resolution',600)

end % plotInsituPhData
