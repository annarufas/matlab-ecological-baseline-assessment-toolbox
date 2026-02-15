% function plotCefasSmartBuoyData(SB)

% PLOTCEFASSMARTBUOYDATA Create various plots to visualise the SmartBuoy dataset.
%
%   INPUT:
%       SB - SmartBuoy dataset
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 19 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Preamble

buoys = categories(SB.buoy);
nBuoys = numel(buoys);

tracers = {'CHLOR','FTU','TOXN','SILICA','O2CONC','TEMP','SAL'};
labelTracers = {'Chl a (\mug/L)','Turbidity (FTU)','Nitrate (\mumol/L)',...
    'Silica (\mumol/L)','Oxygen (mg/L)','Temperature (ºC)','Salinity (PSU)'}; % 'FTU' = 'Formazin Turbidity Units'
nTracers = length(tracers);

% Notice that some buoys have been measuring at different depths over time.
% We will use these unique depths to plot the data (each buoy has a
% distinct depth of deployment)
depthCategories = unique(SB.depth);
nDepthCategories = numel(depthCategories);
% DOWSING: 1 m and 15 m 
% NTHDOGGER: 25 m and 31 m 
% OYSTER: 35 m
% OYSTERML: 45 m

%% Plot time series

colourCategories = brewermap(nDepthCategories,'*Paired');
labelLegend = {'Dowsing, 1 m','Dowsing, 15 m','North Dogger, 25 m',...
    'North Dogger, 31 m','Oyster Grounds, 35 m','Oyster Grounds, 45 m'};

startDate = min(SB.dateTime);
endDate = max(SB.dateTime);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.85],'Color','w') 
haxis = zeros(nTracers,1);

for iTracer = 1:nTracers
    
    haxis(iTracer) = subaxis(nTracers,1,iTracer,'Spacing',0.01,'Padding',0.01,'Margin', 0.10);
    ax(iTracer).pos = get(haxis(iTracer),'Position');
    ax(iTracer).pos(2) = ax(iTracer).pos(2) + 0.03; 
    set(haxis(iTracer),'Position',ax(iTracer).pos) 
    
    for iDepthCategory = 1:nDepthCategories
        thisBuoyTracer = find(SB.variable == tracers(iTracer) & SB.depth == depthCategories(iDepthCategory));
        scatter(haxis(iTracer),SB.dateTime(thisBuoyTracer),SB.value(thisBuoyTracer),...
            6,colourCategories(iDepthCategory,:),'filled'); hold on;
    end
    hold off;
    
    xlim([startDate endDate])    
%     xtickangle(45)

    % Log scale for turbidity
    if (iTracer == 2)
        set(gca,'YScale','log')
        ylim([0 1000])
        yticks([1e-1 1 10 1e2 ])
        yticklabels({'0.1','1','10','100'})
    end
    ylh = ylabel(labelTracers(iTracer),'FontSize',11);
    ylh.Position(1) = -1300;

    if (iTracer == 1)
        title('CEFAS SmartBuoy dataset','FontSize',12)
    end
        
    if (iTracer == nTracers)
        xlabel('Time (date)','FontSize',11);
        lg = legend(labelLegend(:));  
        lg.Orientation = 'vertical';
        lg.Position(1) = 0.35; lg.Position(2) = 0;
        lg.ItemTokenSize = [15,1];
        lg.FontSize = 11;
        lg.NumColumns = 3;
        set(lg,'Box','off')
    end        

    box on

end 

hold off

exportgraphics(gcf,fullfile('.','figures','SmartBuoy_CEFAS_timeseries_all.png'),'Resolution',600)

clear ax

%% Inspect chlorophyll and turbidity closer

iSubplot = 1;

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.50 0.50],'Color','w') 

for iDepthCategory = 1:nDepthCategories
    
    % Get the data
    thisBuoyChla = find(SB.variable == "CHLOR" & SB.depth == depthCategories(iDepthCategory));    
    thisBuoyTurb = find(SB.variable == "FTU" & SB.depth == depthCategories(iDepthCategory));
    chlaTime = SB.dateTime(thisBuoyChla);
    turbTime = SB.dateTime(thisBuoyTurb);
    chlaVals = SB.value(thisBuoyChla);
    turbVals = SB.value(thisBuoyTurb);
    
    if (~isempty(thisBuoyChla))
        
        haxis(iSubplot) = subaxis(3,1,iSubplot,'Spacing',0.015,'Padding',0.015,'Margin', 0.10);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        if (iSubplot == 1)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.08; 
        elseif (iSubplot == 2)
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.03; 
        else
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.02; 
        end    
        set(haxis(iSubplot),'Position',ax(iSubplot).pos) 
        
        yyaxis left
%         plot(haxis(iSubplot),chlaTime,chlaVals,'-','Color','g','LineWidth',1); hold on;
        scatter(haxis(iSubplot),chlaTime,chlaVals,6,'g','filled'); hold on;
        ylabel('Chla (mg m^{-3}','FontSize',11)
        
        yyaxis right
        scatter(haxis(iSubplot),turbTime,turbVals,6,'r','filled'); hold on;
%         plot(haxis(iSubplot),turbTime,turbVals,'-','Color','r','LineWidth',1); hold on;
        ylabel('Turbidity (FTU)','FontSize',11)
        
        xlim([chlaTime(1) chlaTime(end)])
        if (iSubplot == 3)
            xlabel('Time (date)','FontSize',11)
        end
        title(labelLegend(iDepthCategory),'FontSize',11)
%         ax = gca;
%         ax.YAxis(1).Color = 'k';
%         ax.YAxis(2).Color = 'k';
%         clear ax
%         clear pos
        grid on
        
        iSubplot = iSubplot + 1;
    end
% 
%     if (iDepthCategory == nDepthCategories)
%         lg = legend(labelLegend(:));  
%         lg.Orientation = 'vertical';
%         lg.Position(1) = 0.84; lg.Position(2) = 0.72;
%         lg.ItemTokenSize = [15,1];
%         lg.FontSize = 11;
%         set(lg,'Box','off')
%     end
end
hold off;

% ylim([0 1])
% grid on;
% ax = gca;
% newPosition = ax.Position;
% newPosition(1) = newPosition(1) - 0.08; % Adjust the value based on your needs
% newPosition(2) = newPosition(2) + 0.02; 
% ax.Position = newPosition;

% lg = legend(labelLegend(:));  
% lg.Orientation = 'vertical';
% lg.Position(1) = 0.84; lg.Position(2) = 0.72;
% lg.ItemTokenSize = [15,1];
% lg.FontSize = 11;
% set(lg,'Box','off')
box on

exportgraphics(gcf,fullfile('.','figures','SmartBuoy_CEFAS_timeseries_chla.png'),'Resolution',600) 
clear ax

%% Histogram with all the data

% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.29],'Color','w') 
% 
% % Main plot
% histogram(SB.value(SB.variable == "CHLOR"),50,'BinLimits',[0,20],'FaceColor', [0.5, 0.5, 0.5])
% ylim([1 6e4])
% set(gca,'YScale','log');
% grid on;
% mainAxes = gca;
% mainAxes.XGrid = 'on';  
% mainAxes.YGrid = 'on';  
% mainAxes.MinorGridAlpha = 0;
% % yticks([1, 1e1, 1e2, 1e3]);
% % yticklabels({'1','10','10^{2}','10^{3}'});
% xlabel('Chlorophyll a (mg m^{-3})');
% ylabel('Frequency');
% title('Distribution of SmartBuoy observations from CEFAS');
% 
% % Inset plot
% insetAxes = axes('Position', [0.52, 0.43, 0.36, 0.45]); % define the position of the inset plot
% histogram(insetAxes,SB.value(SB.variable == "CHLOR"),19,'BinLimits',[1e-1,2],'FaceColor', [0.5, 0.5, 0.5])
% ylim([1 4e4])
% set(insetAxes,'YScale','log');
% grid on;
% insetAxes.XGrid = 'on'; 
% insetAxes.YGrid = 'on'; 
% insetAxes.MinorGridAlpha = 0;  
% insetAxes.XAxis.FontSize = 8;
% insetAxes.YAxis.FontSize = 8;
% yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5]);
% yticklabels({'1','10','10^{2}','10^{3}','10^{4}','10^{5}'});
% xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
% ylabel('Frequency','FontSize',8)
% 
% exportgraphics(gcf,fullfile('.','figures','SmartBuoy_CEFAS_histogram.png'),'Resolution',600)

%% Histogram by season

seasons = {'Winter','Spring','Summer','Autumn'};
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
            seasonData = SB.value(SB.season == "Winter" & SB.variable == "CHLOR"); 
        case 'Spring'
            seasonData = SB.value(SB.season == "Spring" & SB.variable == "CHLOR"); 
        case 'Summer'
            seasonData = SB.value(SB.season == "Summer" & SB.variable == "CHLOR"); 
        case 'Autumn'
            seasonData = SB.value(SB.season == "Autumn" & SB.variable == "CHLOR"); 
    end
    
    % Main plot
    
    % Plot all data
    histogram(haxis(iSeason),SB.value(SB.variable == "CHLOR"),50,'BinLimits',[0,20],'FaceColor',[0.7,0.7,0.7],'EdgeColor','none')
    hold on
    % Plot seasonal data
    histogram(haxis(iSeason),seasonData,50,'BinLimits',[0,20],'FaceColor',colourScheme(iSeason,:))
    hold on
    ylim([1 6e4])
    set(gca,'YScale','log');
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
    insetSeasonalAxes = axes('Position', [ax(iSeason).pos(1) + ax(iSeason).pos(3)*0.45,... 
                                          ax(iSeason).pos(2) + ax(iSeason).pos(4)*0.47,... 
                                          0.52*ax(iSeason).pos(3),... 
                                          0.49*ax(iSeason).pos(4)]);
    
    histogram(insetSeasonalAxes,seasonData,19,'BinLimits',[1e-1,2],'FaceColor',colourScheme(iSeason,:))
    hold on
    ylim([1 1e4])
    set(gca,'YScale','log');
    grid on;
    insetSeasonalAxes.XGrid = 'on';
    insetSeasonalAxes.YGrid = 'on';
    insetSeasonalAxes.MinorGridAlpha = 0;
    insetSeasonalAxes.XAxis.FontSize = 8;
    insetSeasonalAxes.YAxis.FontSize = 8;
    yticks([1, 1e1, 1e2, 1e3, 1e4]);
    yticklabels({'1','10','10^{2}','10^{3}','10^{4}'});
    xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
    ylabel('Frequency','FontSize',8)  
    
end
hold off

exportgraphics(gcf,fullfile('.','figures','SmartBuoy_CEFAS_histogram_byseason.png'),'Resolution',600)
clear ax

% end % plotCefasSmartBuoyData
