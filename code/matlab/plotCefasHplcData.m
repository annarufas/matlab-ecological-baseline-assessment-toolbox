function plotCefasHplcData(pathEnduranceShapefile,pathAreaStudyShapefile,HPLC)

% PLOTHPLCDATA Create various plots to visualise the HPLC dataset.
%
%   INPUT:
%       pathEnduranceShapefile - shapefile with the Endurance GCS site
%       pathAreaStudyShapefile - shapefile with our area of study
%       HPLC                   - HPLC data filtered
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 7 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Definitions

seasons = {'Winter','Spring','Summer','Autumn'};
cmin = log10(0.5); %log10(min(P.TChlA_ug_L));
cmax = log10(50); %log10(max(P.TChlA_ug_L));
c = [0.5 5 50];

%% Map locations

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.40],'Color','w') 

setNorthSeaMappingFeatures(HPLC,'TChlA_ug_L','Longitude','Latitude',...
    cmin,cmax,pathEnduranceShapefile,pathAreaStudyShapefile,...
    sprintf('CEFAS HPLC (0–30 m), \\itN\\rm \\bf= %d',height(HPLC)),'center')

cb = colorbar('Location','southoutside');
cb.Label.String = 'TChla (mg m^{-3})';
cb.YTick = log10(c);
cb.YTickLabel = c;

exportgraphics(gcf,fullfile('.','figures','HPLC_map_all.png'),'Resolution',600)

%% Map locations by season

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.60],'Color','w') 
haxis = zeros(length(seasons),1);

for iSeason = 1:length(seasons)
    
    haxis(iSeason) = subaxis(2,2,iSeason,'Spacing',0.020,'Padding',0.030,'Margin',0.07);
    ax(iSeason).pos = get(haxis(iSeason),'Position');
    
    % Get the current season data and set the plot title
    seasonName = seasons{iSeason};
    thisSeasonHplcData = HPLC(HPLC.season == seasonName, :);
    titleString = sprintf('%s, \\itN\\rm \\bf= %d', seasonName, height(thisSeasonHplcData));

    % Plot the data for the current season
    setNorthSeaMappingFeatures(thisSeasonHplcData,'TChlA_ug_L','Longitude','Latitude',...
        cmin,cmax,pathEnduranceShapefile,pathAreaStudyShapefile,...
        titleString,'left')

    hold on

end

ax(2).pos(1) = ax(2).pos(1) - 0.08;
ax(4).pos(1) = ax(4).pos(1) - 0.08;
for iSubplot = 1:4
    % Shift subplots upwards to accommodate legend
    ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.05;
    set(haxis(iSubplot),'Position',ax(iSubplot).pos) 
end
       
% Common legend to all subplots
cb = colorbar('Location','southoutside');
cb.Label.String = 'TChla (mg m^{-3})';
cb.YTick = log10(c);
cb.YTickLabel = c;
cb.Position(1) = cb.Position(1) - 0.405;
cb.Position(2) = cb.Position(2) - 0.12;
cb.Position(3) = 0.635; % length
cb.Position(4) = 0.025; % width

exportgraphics(gcf,fullfile('.','figures','HPLC_map_byseason.png'),'Resolution',600)  

%% Histogram with all the data

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.29],'Color','w') 

% Main plot
histogram(HPLC.TChlA_ug_L(:),42,'BinLimits',[0,42],'FaceColor', [0.5, 0.5, 0.5])
grid on;
mainAxes = gca;
mainAxes.XGrid = 'on'; 
mainAxes.YGrid = 'on'; 
mainAxes.MinorGridAlpha = 0; % set alpha (transparency) to 0 for secondary grid lines
xlabel('Chlorophyll a (mg m^{-3})');
ylabel('Frequency');
title('Distribution of HPLC observations from CEFAS');

% Inset plot
insetAxes = axes('Position', [0.4, 0.43, 0.47, 0.45]); 
histogram(insetAxes,HPLC.TChlA_ug_L(:),9,'BinLimits',[1e-1,1],'FaceColor', [0.5, 0.5, 0.5])
grid on;
insetAxes.XGrid = 'on';  
insetAxes.YGrid = 'on';  
insetAxes.MinorGridAlpha = 0;
insetAxes.XAxis.FontSize = 8;
insetAxes.YAxis.FontSize = 8;
xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
ylabel('Frequency','FontSize',8)

exportgraphics(gcf,fullfile('.','figures','HPLC_histogram_all.png'),'Resolution',600)

%% Histogram with the data by season

colourScheme = brewermap(length(seasons),'*Spectral');

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.50],'Color','w') 
haxis = zeros(2,2);

for iSeason = 1:length(seasons)
    
    haxis(iSeason) = subaxis(2,2,iSeason,'Spacing',0.018,'Padding',0.020,'Margin',0.07);
    
    ax(iSeason).pos = get(haxis(iSeason),'Position');
    if (iSeason == 1 || iSeason == 2)
        ax(iSeason).pos(2) = ax(iSeason).pos(2) + 0.05;
    end
    ax(iSeason).pos(1) = ax(iSeason).pos(1) + 0.01;
    set(haxis(iSeason),'Position',ax(iSeason).pos) 
 
    switch seasons{iSeason}
        case 'Winter'
            thisSeasonHplcData = HPLC.TChlA_ug_L(HPLC.season == "Winter");
        case 'Spring'
            thisSeasonHplcData = HPLC.TChlA_ug_L(HPLC.season == "Spring");
        case 'Summer'
            thisSeasonHplcData = HPLC.TChlA_ug_L(HPLC.season == "Summer");
        case 'Autumn'
            thisSeasonHplcData = HPLC.TChlA_ug_L(HPLC.season == "Autumn");
    end
    
    % Main plot
    
    % Plot all data
    histogram(haxis(iSeason),HPLC.TChlA_ug_L(:),42,'BinLimits',[0,42],'FaceColor',[0.7,0.7,0.7],'EdgeColor','none')
    hold on
    % Plot seasonal data
    histogram(haxis(iSeason),thisSeasonHplcData,42,'BinLimits',[0,42],'FaceColor',colourScheme(iSeason,:))
    hold on
    ylim([0 400])
    grid on;
    mainSeasonalAxes = gca;
    mainSeasonalAxes.XGrid = 'on';  
    mainSeasonalAxes.YGrid = 'on';  
    mainSeasonalAxes.MinorGridAlpha = 0;  

    if (iSeason == 1 || iSeason == 3)
        ylabel('Frequency');
    end
    if (iSeason == 3 || iSeason == 4)
        xlabel('Chlorophyll a (mg m^{-3})');
    end
    title(seasons{iSeason});
    
    % Inset plot
    
    insetSeasonalAxes = axes('Position', [ax(iSeason).pos(1) + ax(iSeason).pos(3)*0.40,... 
                                          ax(iSeason).pos(2) + ax(iSeason).pos(4)*0.43,... 
                                          0.52*ax(iSeason).pos(3),... 
                                          0.49*ax(iSeason).pos(4)]);
    
    histogram(insetSeasonalAxes,thisSeasonHplcData,14,'BinLimits',[1e-1,1.5],'FaceColor',colourScheme(iSeason,:))
    hold on
    ylim([0 60])
    grid on;
    insetSeasonalAxes.XGrid = 'on';  
    insetSeasonalAxes.YGrid = 'on'; 
    insetSeasonalAxes.MinorGridAlpha = 0;  
    insetSeasonalAxes.XAxis.FontSize = 8;
    insetSeasonalAxes.YAxis.FontSize = 8;
    xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
    ylabel('Frequency','FontSize',8) 
    
end
hold off

exportgraphics(gcf,fullfile('.','figures','HPLC_histogram_byseason.png'),'Resolution',600)

end % plotCefasHplcData
