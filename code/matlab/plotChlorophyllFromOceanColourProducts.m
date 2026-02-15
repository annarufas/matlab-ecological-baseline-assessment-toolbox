% function plotChlorophyllFromOceanColourProducts(AScmems,ASnasa,ASoccci,...
%     pathAreaStudyShapefile,pathEnduranceShapefile)

% PLOTCHLOROPHYLLFROMOCEANCOLOURPRODUCTS Explore chlorophyll from various
% ocean colour products downloaded (CMEMS, NASA and OC-CCI)
%
%   INPUT:
%       AScmems                - table with time-series data from CMEMS for our area of study
%       ASnasa                 - table with time-series data from NASA for our area of study
%       ASoccci                - table with time-series data from OC-CCI for our area of study
%       pathAreaStudyShapefile - shapefile with our area of study
%       pathEnduranceShapefile - shapefile with the Endurance GCS site
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

%% Presets

% Extract chlorphyll from OC products
[chlaStruct,prodLabel,prodName,nChlaProducts] =... 
    extractChlorophyllFromOceanColourProducts(AScmems,ASnasa,ASoccci);

% Cap calculations to this date
targetDateTimeEnd = datetime(2023,12,31);

%% Compute average data for our area of study

for iProduct = 1:numel(prodName)
            
    thisProductName = prodName{iProduct};

    % Filter to include only dates up to and including targetDateTimeEnd
    time = chlaStruct.(thisProductName).time;
    filteredTimeVector = time(time <= targetDateTimeEnd);
    nTimeSteps = length(filteredTimeVector);
    time = time(1:nTimeSteps);
    chlaVar = chlaStruct.(thisProductName).dataset(:,:,1:nTimeSteps);

    % Calc scene-averaged data
    if (isfield(chlaStruct.(thisProductName),'depth') == 0) % if the dataset has no depth dependency
        [chlaSceneAverage,~] = calcSceneAverage(chlaVar,time);
    else % if the dataset has depth-dependency, compute the pondered depth average 
        nDepthLevels = numel(plkStruct.(thisProductName).depth);
        chlaSceneAverage = calcScenePonderedDepthAverage(chlaVar,time,nDepthLevels,nTimeSteps);       
    end

    % Add fields to the structure array
    chlaStruct.(thisProductName) = setfield(chlaStruct.(thisProductName),'CHL_time',time);
    chlaStruct.(thisProductName) = setfield(chlaStruct.(thisProductName),'CHL_dailyMean',chlaSceneAverage);

end

%% Plot all satellite data together in a time axis 

% Limits for time axis
dmin = datetime(1997,09,01);
dmax = targetDateTimeEnd;
   
% Ticks for time axis
xtickTimes = datetime(1998,01,01):calyears(2):targetDateTimeEnd;
xtickTimes.Format = 'yyyy'; % use 'yyyy' to display only the year
xtickLabels = string(xtickTimes);

% Chla break point
yChlaBreakValue = 3; % mg m-3

% Colours
colourChlProducts = colormap(brewermap(nChlaProducts,'*Spectral')); 

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.29],'Color','w') 

for iProduct = 1:nChlaProducts
 
    % Extract variables
    thisProductName = prodName{iProduct};
    y = chlaStruct.(thisProductName).CHL_dailyMean;
    x = chlaStruct.(thisProductName).CHL_time;
    
    % Create axes for the bottom part of the plot
    if (iProduct == 1)
        ax1 = axes('Position', [0.1 0.12 0.85 0.40]); % position in within the figure axes (DO NOT CHANGE)
        ax1.XAxis.FontSize = 12;
        ax1.YAxis.FontSize = 12;
    end
    hold on
    
    %plot(ax1,x,y,'Color',colourChlDataset(i,:),'LineWidth',widthLineChlDataset(i));
    scatter(ax1,x,y,4,colourChlProducts(iProduct,:),'filled');

    xlim([dmin dmax])
    ylim([0 yChlaBreakValue])
    if (iProduct == 1)
        xticks(xtickTimes)
        xticklabels(xtickLabels)
        xtickangle(45)
    end

    % Create axes for the upper part of the plot
    if (iProduct == 1)
        ax2 = axes('Position', [0.1 0.52 0.85 0.42]); % position in within the figure axes (DO NOT CHANGE)
        ax2.XAxis.Visible = 'off';
        ax2.YScale = 'log';
        ax2.XAxis.FontSize = 12;
        ax2.YAxis.FontSize = 12;
    end
    hold on
    
    %plot(ax2,x,y,'Color',colourChlDataset(i,:),'LineWidth',widthLineChlDataset(i));
    scatter(ax2,x,y,4,colourChlProducts(iProduct,:),'filled');

    xlim([dmin dmax])
    ylim([yChlaBreakValue 100])
    if (iProduct == 1)
        yticks([yChlaBreakValue 10 100])
        yticklabels({num2str(yChlaBreakValue),'10','100'})
    end

    hold on 
    
end
hold off

% xlabel('Time')
yl = ylabel('Chlorophyll a (mg m^{-3})','FontSize',12);
yl.Position(2) = 3.5; % change vertical position of ylabel

lg = legend(prodLabel);
% lg.Orientation = 'vertical';
lg.NumColumns = 2;
lg.Position(1) = 0.10; lg.Position(2) = 0.75;
lg.ItemTokenSize = [15,1];
lg.FontSize = 9;
set(lg,'Box','on')

exportgraphics(gcf,fullfile('.','figures','chla_timeseries.png'),'Resolution',600)
clear ax1 ax2

%% Plot separately – show outliers removal method

for iProduct = 1:nChlaProducts

    % Extract variables
    thisProductName = prodName{iProduct};
    y = chlaStruct.(thisProductName).CHL_dailyMean;
    x = chlaStruct.(thisProductName).CHL_time;

    % Detect and replace outliers in data
    % Define outliers as points more than three local scaled MAD away from 
    % the local median within a sliding window. Find the location of the 
    % outlier in y relative to the points in x with a window size of 7 days. 
    % Fill the outlier with the computed threshold value using the method 
    % 'clip', and plot the original and filled data.
    yFill = filloutliers(y,'pchip','movmedian',days(7),'SamplePoints',x);
    
    % Prepare axes environment
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.41 0.27],'Color','w') 
    ax = axes('Position',[0.09 0.14 0.70 0.74]); % position in within the figure axes (DO NOT CHANGE)
    hold on
    
    plot(ax,x,yFill,'-k','LineWidth',1); hold on;
    scatter(ax,x,y,4,colourChlProducts(iProduct,:),'filled'); hold on;
    
    xlim([x(1) x(end)])
    ylim([0 max(y)])
    xtickangle(45)

    yl = ylabel('Chlorophyll a (mg m^{-3})','FontSize',12);
    yl.Position(1) = yl.Position(1) - 100;
    title(prodLabel(iProduct),'FontSize',12)
    box on

    lg = legend('With filled outliers','With outliers','FontSize',12);
    lg.Orientation = 'vertical';
    lg.Position(1) = 0.78; lg.Position(2) = 0.74;
    lg.ItemTokenSize = [15,1];
    set(lg,'Box','off')
    
    set(ax,'FontSize',12)
    
    exportgraphics(gcf,fullfile('.','figures',...
        strcat('chla_timeseries_filled',prodName{iProduct},'.png')),'Resolution',600)
    clear ax

end

%% Plot separately – show trends on data with outliers removed

for iProduct = 1:nChlaProducts

    % Extract variables
    thisProductName = prodName{iProduct};
    y = chlaStruct.(thisProductName).CHL_dailyMean;
    x = chlaStruct.(thisProductName).CHL_time;

    % Detect and replace outliers in data
    yFill = filloutliers(y,'clip','movmedian',days(7),'SamplePoints',x);

    % Fit a linear trend to the data once all NaNs have been removed
    xd = datenum(x); % transform datetime into double
    [p,S,mu] = polyfit(xd(~isnan(yFill)),yFill(~isnan(yFill)),1);
    yTrend = polyval(p,xd,S,mu);
    
    % Calculate seasonal cycle (use the 'climatology' function in the 
    % Climate Data Toolbox for Matlab) 
    ynonnan = yFill;
    ynonnan(isnan(ynonnan)) = 0;
    ySeasonalCycle = climatology(ynonnan,x,'full');

    % Prepare axes environment
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.39 0.27],'Color','w')
    ax = axes('Position', [0.09 0.14 0.70 0.74]); % position in within the figure axes (DO NOT CHANGE)
    hold on
    
    p1 = plot(ax,x,ySeasonalCycle,'Color','#5aae61','LineWidth',1.0); hold on;
    p2 = scatter(ax,x,yFill,2,'k','filled'); hold on;
    p3 = plot(ax,x,yTrend,'Color',[0 0 1],'LineWidth',1.5); hold off; % same as polyplot

    xlim([x(1) x(end)])
    ylim([0 max(ySeasonalCycle)+max(ySeasonalCycle)*0.50]) % ylim([0 max(ySeasonalCycle)+max(ySeasonalCycle)*0.20])
    xtickangle(45)
    ytickformat('%.1f')

    yl = ylabel('Chlorophyll a (mg m^{-3})','FontSize',12);
    yl.Position(1) = yl.Position(1) - 70;

    title(prodLabel(iProduct),'FontSize',12)
    box on

    lg = legend([p2 p1 p3], {'Data','Seasonal cycle','Trend'},'FontSize',12);
    lg.Orientation = 'vertical';
    lg.Position(1) = 0.79; lg.Position(2) = 0.69;
    lg.ItemTokenSize = [15,1];
    set(lg,'Box','off')
    
    set(ax,'FontSize',12)
   
    exportgraphics(gcf,fullfile('.','figures',...
        strcat('chla_timeseries_trends',prodName{iProduct},'.png')),'Resolution',600)
    clear ax

end

%% Comparison of two products with different resolution

% Initialise structures
prodA = struct();
prodB = struct();

for iProduct = 1:nChlaProducts
    thisProductName = prodName{iProduct};
    switch thisProductName
        case 'global_occcimulti_r1km'
            prodA.chl = chlaStruct.(thisProductName).CHL_dailyMean;
            prodA.time = chlaStruct.(thisProductName).CHL_time;
            prodA.label = prodLabel{iProduct};
            % Detect and replace outliers in data
            prodA.chlFilled = filloutliers(prodA.chl,'clip','movmedian',days(7),'SamplePoints',prodA.time); 
        case 'global_occcimulti_r4km'
            prodB.chl = chlaStruct.(thisProductName).CHL_dailyMean;
            prodB.time = chlaStruct.(thisProductName).CHL_time;
            prodB.label = prodLabel{iProduct};
            % Detect and replace outliers in data
            prodB.chlFilled = filloutliers(prodB.chl,'clip','movmedian',days(7),'SamplePoints',prodB.time); 
    end
end
            
% Find concurrent dates: which dates in prodB.time are also in prodA.time
[isMatchBinA, ~] = ismember(prodB.time, prodA.time);
% and viceversa
[isMatchAinB, ~] = ismember(prodA.time, prodB.time);

% Filter prodA and prodB to keep only the matching dates
prodB.timeMatched = prodB.time(isMatchBinA);
prodB.chlFilledMatched = prodB.chlFilled(isMatchBinA);
prodA.timeMatched = prodA.time(isMatchAinB);
prodA.chlFilledMatched = prodA.chlFilled(isMatchAinB);

% Remove concurrent entries with NaN values
X = prodA.chlFilledMatched;
Y = prodB.chlFilledMatched;
highRes = X(~isnan(X) & ~isnan(Y));
lowRes = Y(~isnan(X) & ~isnan(Y));

% Calculate statistics
columnHeaders = {'N','slope','intercept','R2','r','RMSE','log_RMSE',...
    'ME','log_ME','MAE','log_MAE','MRPE','MAPE'};
variableTypes = {'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double'};
T = table('Size',[1,numel(columnHeaders)],...
    'VariableNames',columnHeaders,'VariableTypes',variableTypes);
[T{1,:}] = calculateMatchupRegressionStatistics(highRes,lowRes);

% Prepare figure environment
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.25 0.35],'Color','w') 
ax = axes('Position', [0.18 0.20 0.77 0.70]); % position in within the figure axes (DO NOT CHANGE)

c = getScatterHeatPlotColormapIndices(highRes,lowRes,80);
scatter(highRes,lowRes,50,c,'.');
hold on
colormap(brewermap(1000,'*YlGnBu'))
cb = colorbar;
cb.Label.String = 'Density of points';
set(cb,'FontSize',12);
caxis([0 80]);
set(gca,'Xscale','log','Yscale','log')
box on

xlim([0.1 15]); ylim([0.1 15])   
xticks([1e-1 0.5 1 2 3 5 10]); yticks([1e-1 0.5 1 2 3 5 10])
xticklabels({'0.1','0.5','1','2','3','5','10'}); yticklabels({'0.1','0.5','1','2','3','5','10'})

yl = ylabel({['Chl a (mg m^{-3})'],... 
             [prodB.label]},'FontSize',12);
xl = xlabel({['Chl a (mg m^{-3})'],... 
             [prodA.label]},'FontSize',12);
yl.Position(1) = yl.Position(1) - 0.01;
xl.Position(2) = xl.Position(2) - 0.005;
        
% Plot 1:1 reference line
xref = 10.^(-1.8:.1:1.8); yref = xref;
plot(xref,yref,'k-')
hold off;
        
xt = min(xlim)+0.15*min(xlim); 
yt = max(ylim)-0.15*max(ylim);
    
% Add statistics
text(xt, yt,... % text position relative to axis
     {[cell2mat(strcat({'R^{2} = '},num2str(T.R2,'%.2f')))],...
     [cell2mat(strcat({'RMSE = '},num2str(T.RMSE,'%.2f')))],...
     [cell2mat(strcat({'\delta = '},num2str(T.ME,'%.2f')))],...
     [cell2mat(strcat({'{\it N} = '},num2str(T.N,'%.0f')))]},...
    'FontSize',12,'Horiz','left','Vert','top');
    
set(ax,'FontSize',12)
box on

exportgraphics(gcf,fullfile('.','figures','chla_comparison_HR_vs_LR.png'),'Resolution',600)
clear prodA prodB

%% Comparison of two products with same resolution

% Initialise structures
prodA = struct();
prodB = struct();

% Notice prodA.lat = prodB.lat
for iProduct = 1:nChlaProducts
    thisProductName = prodName{iProduct};
    switch thisProductName
        case 'global_olci_r4km'
            prodA.chl = chlaStruct.(thisProductName).dataset;
            prodA.time = chlaStruct.(thisProductName).time;
            prodA.lat = chlaStruct.(thisProductName).lat;
            prodA.lon = chlaStruct.(thisProductName).lon;
            prodA.label = prodLabel{iProduct};
        case 'global_aquamodis_r4km'
            prodB.chl = chlaStruct.(thisProductName).dataset;
            prodB.time = chlaStruct.(thisProductName).time;
            prodB.lat = chlaStruct.(thisProductName).lat;
            prodB.lon = chlaStruct.(thisProductName).lon;
            prodB.label = prodLabel{iProduct};
    end
end

% Use ismember to find which dates in prodA.time are also present in prodB.time
[isMatchAinB, ~] = ismember(prodA.time, prodB.time);
% and viceversa
[isMatchBinA, ~] = ismember(prodB.time, prodA.time);

% Filter prodA and prodB to keep only the matching dates
prodA.timeMatched = prodA.time(isMatchAinB);
prodA.chlMatched = prodA.chl(:,:,isMatchAinB);
prodB.timeMatched = prodB.time(isMatchBinA);
prodB.chlMatched = prodB.chl(:,:,isMatchBinA);

% Average data over time
prodB.chlMatchedMean = mean(prodB.chlMatched(:,:,:),3,'omitnan');
prodA.chlMatchedMean = mean(prodA.chlMatched(:,:,:),3,'omitnan');

% Subtract
newSat = prodA.chlMatchedMean;
oldSat = prodB.chlMatchedMean;
chlaDiff = newSat(:,:) - oldSat(:,:);

areaStudy = m_shaperead(pathAreaStudyShapefile); % lat/lon coordinates
enduranceSite = m_shaperead(pathEnduranceShapefile); % lat/lon coordinates

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.27 0.35],'Color','w') 
ax = axes('Position', [0.10 0.10 0.77 0.70]); % position in within the figure axes (DO NOT CHANGE)
    
m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy);
m_pcolor(prodA.lon,prodA.lat,chlaDiff) 
shading flat
colormap(brewermap(1000,'*RdBu'))
caxis([-1 1])
m_gshhs_i('color','k'); % intermediate resolution coastline
m_grid('linewi',2,'tickdir','out','FontSize',12);
m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','r');

cb = colorbar;
cb.Label.String = 'Chlorophyll a difference (mg m^{-3})';
set(cb,'ytick',[-1:0.2:1],'yticklabel',[-1:0.2:1],...
    'tickdir','out','FontSize',12);
 
tl = title({['[',prodA.label,']',' –'],... 
            ['[',prodB.label,']']},'FontSize',12);
tl.Position(2) = tl.Position(2) + 0.001;

set(ax,'FontSize',12)

exportgraphics(gcf,fullfile('.','figures','chla_comparison_OLCI_vs_MODIS.png'),'Resolution',600)

% end % plotChlorophyllFromOceanColourProducts