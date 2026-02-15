function plotAnomalyPlotsForDrivers(AScmems,ASbicep,ASnasa,ASoccci,datasetConfig)

% PLOTANOMALYPLOTSFORDRIVERS Explore driver data using cycle plots as well 
% bar anomaly plots.
%
%   INPUT:
%       AScmems       - table with time-series data from CMEMS for our area of study
%       ASnasa        - table with time-series data from NASA for our area of study
%       ASoccci       - table with time-series data from OC-CCI for our area of study
%       datasetConfig - shapefile with our area of study
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

% Get driver data
[S,~,~,~] = extractDriverDataFromProducts(AScmems,ASbicep,ASnasa,ASoccci);

% Dataset names in S
[datasetNames] = getTimeSeriesDatasetNames(S);

% Dataset indexes
idxsCmemsDatasetsFirstRound = 1:8; % change as necessary
idxsCmemsDatasetsSecondRound = 9:16; % change as necessary
idxsBicepDatasets = 17:22; % change as necessary
idxsChlaDatasets = 23:30; % change as necessary   
         
%% Compress the data into cycles

% CMEMS
nSubplots = 8;
for iRound = 1:2
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.70],'Color','w')
    haxis = zeros(nSubplots,1);
    if (iRound == 1)
        idxsDatasets = idxsCmemsDatasetsFirstRound; 
    else
        idxsDatasets = idxsCmemsDatasetsSecondRound;
    end 
    for iDataset = 1:nSubplots
        thisDatasetIdx = idxsDatasets(iDataset); 
        haxis(iDataset) = subaxis(4,2,iDataset,'Spacing',0.015,'Padding',0.015,'Margin',0.06);
        setCyclePlotFeatures(haxis,iDataset,S,thisDatasetIdx,datasetConfig,nSubplots,datasetNames);
    end
    exportgraphics(gcf,fullfile('.','figures',strcat('cycles_cmems_',num2str(iRound),'.png')),'Resolution',600)
    clear ax
end

% BICEP
nSubplots = 6;
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.60],'Color','w')
haxis = zeros(nSubplots,1);
for iDataset = 1:nSubplots
    thisDatasetIdx = idxsBicepDatasets(iDataset); 
    haxis(iDataset) = subaxis(3,2,iDataset,'Spacing',0.015,'Padding',0.015,'Margin',0.06);
    setCyclePlotFeatures(haxis,iDataset,S,thisDatasetIdx,datasetConfig,nSubplots,datasetNames);
end
exportgraphics(gcf,fullfile('.','figures','cycles_bicep.png'),'Resolution',600)
clear ax

% Chla
nSubplots = 8;
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.70],'Color','w')
haxis = zeros(nSubplots,1);
for iDataset = 1:nSubplots
    thisDatasetIdx = idxsChlaDatasets(iDataset); 
    haxis(iDataset) = subaxis(4,2,iDataset,'Spacing',0.015,'Padding',0.015,'Margin',0.06);
    setCyclePlotFeatures(haxis,iDataset,S,thisDatasetIdx,datasetConfig,nSubplots,datasetNames);
end
exportgraphics(gcf,fullfile('.','figures','cycles_chla.png'),'Resolution',600)
clear ax

%% Anomaly bar plots showing deviations from average    

% CMEMS
nSubplots = 8;
for iRound = 1:2
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.70],'Color','w')
    haxis = zeros(nSubplots,1);
    if (iRound == 1)
        idxsDatasets = idxsCmemsDatasetsFirstRound; 
    else
        idxsDatasets = idxsCmemsDatasetsSecondRound;
    end 
    for iDataset = 1:nSubplots
        thisDatasetIdx = idxsDatasets(iDataset); 
        haxis(iDataset) = subaxis(4,2,iDataset,'Spacing',0.015,'Padding',0.015,'Margin',0.06);
        setAnomalyPlotFeatures(haxis,iDataset,S,thisDatasetIdx,datasetConfig,nSubplots,datasetNames)
    end % iDataset
    exportgraphics(gcf,fullfile('.','figures',strcat('anomalies_cmems_',num2str(iRound),'.png')),'Resolution',600)
    clear ax
end

% BICEP
nSubplots = 6;
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.60],'Color','w')
haxis = zeros(nSubplots,1);
for iDataset = 1:nSubplots
    thisDatasetIdx = idxsBicepDatasets(iDataset); 
    haxis(iDataset) = subaxis(3,2,iDataset,'Spacing',0.015,'Padding',0.015,'Margin',0.06);
    setAnomalyPlotFeatures(haxis,iDataset,S,thisDatasetIdx,datasetConfig,nSubplots,datasetNames)
end
exportgraphics(gcf,fullfile('.','figures','anomalies_bicep.png'),'Resolution',600)
clear ax

% Chla
nSubplots = 8;
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.70],'Color','w')
haxis = zeros(nSubplots,1);
for iDataset = 1:nSubplots
    thisDatasetIdx = idxsChlaDatasets(iDataset); 
    haxis(iDataset) = subaxis(4,2,iDataset,'Spacing',0.015,'Padding',0.015,'Margin',0.06);
    setAnomalyPlotFeatures(haxis,iDataset,S,thisDatasetIdx,datasetConfig,nSubplots,datasetNames)
end
exportgraphics(gcf,fullfile('.','figures','anomalies_chla.png'),'Resolution',600)
clear ax

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function setCyclePlotFeatures(haxis,iDataset,S,thisDatasetIdx,datasetConfig,...
    nSubplots,datasetNames)

monthNames = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

% Move subplots
ax(iDataset).pos = get(haxis(iDataset),'Position');
if (iDataset == 1 || iDataset == 2)
    ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.045;
elseif (iDataset == 3 || iDataset == 4)
    ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.015;
elseif (iDataset == 5 || iDataset == 6)
    ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.015;
elseif (iDataset == 7 || iDataset == 8)
    ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.045;  
end  
if (iDataset == 2 || iDataset == 4 || iDataset == 6 || iDataset == 8)
    ax(iDataset).pos(1) = ax(iDataset).pos(1) + 0.06;
else
    ax(iDataset).pos(1) = ax(iDataset).pos(1) + 0.005;
end
set(haxis(iDataset),'Position',ax(iDataset).pos)
               
% Extract variables
data = S(thisDatasetIdx).sceneAverage;
time = S(thisDatasetIdx).time;

% Convert the datetime vector to just dates (with time set to midnight)
dateOnlyVector = dateshift(time, 'start', 'day');

% Evaluate if the data is daily or not
dateDifferences = days(diff(dateOnlyVector));
uniqueDifferences = unique(dateDifferences); % check the mode of the date differences
if all(uniqueDifferences == 1)
    frequency = 'daily regular';
elseif all(ismember(uniqueDifferences, [28,29,30,31]))
    frequency = 'monthly';
elseif any(uniqueDifferences == 1)
    frequency = 'daily irregular';
else
    frequency = 'other';
end

% Plot       
if (strcmp(frequency,'daily regular') || strcmp(frequency,'daily irregular')) 
    scatter(haxis(iDataset),doy(datenum(time)),data,10,datenum(time),'filled')
    xlim([1 365])
    xlabelStr = 'Day of the year';
elseif (strcmp(frequency,'monthly'))
    scatter(haxis(iDataset),month(time),data,20,datenum(time),'filled')
    xlim([1 12]) 
    xlabelStr = 'Month of the year';
    xticks(1:1:12); 
    xticklabels(monthNames); 
        
end
%         minVar = cmemsDatasetConfig(thisDatasetIdx).minVar;
%         maxVar = cmemsDatasetConfig(thisDatasetIdx).maxVar;
%         ylim([minVar maxVar])

yearsInRange = unique(year(time));
nYearsInRange = numel(yearsInRange);
cmap = brewermap(nYearsInRange,'*GnBu');
colormap(cmap)

tickAndLabelFreq = 3;
cb = colorbar;
cbdate('yyyy')
dTk = diff(cb.Limits)/(tickAndLabelFreq*length(cmap));
cbTicks = cb.Limits(1)+dTk:tickAndLabelFreq*dTk:cb.Limits(2)-dTk;
% Only show tick labels every 3 years when there are more than 12
% years on display, so that the colour bar does not look too crammed
cbTickLabels = cell(size(yearsInRange));
if (nYearsInRange > 12)
    % Set labels for every three years and leave the rest empty
    for i = 1:length(yearsInRange)
        if mod(i - 1, tickAndLabelFreq) == 0 % Check if the index is in our freq (0, 3, 6, ...)
            cbTickLabels{i} = num2str(yearsInRange(i));
        else
            cbTickLabels{i} = ''; % Leave the label empty
        end
    end
else
    cbTickLabels = string(yearsInRange);
end
set(cb,'Ticks',cbTicks,'TickLabels',cbTickLabels)

yl = ylabel(datasetConfig(thisDatasetIdx).unitsVar,'FontSize',10);
% yl.Position(1) = yl.Position(1) - 10; % change vertical position of ylabel
ytickformat(datasetConfig(thisDatasetIdx).tickLabelFmt);

if (iDataset == nSubplots)
    xl = xlabel(xlabelStr,'FontSize',12);
    xl.Position(1) = xl.Position(1) - 1.3*xl.Position(1); 
    xl.Position(2) = xl.Position(2) - 0.05; 
end

box on

title(datasetNames{thisDatasetIdx},'FontSize',12);

end % setCyclePlotFeatures

% *************************************************************************

function setAnomalyPlotFeatures(haxis,iDataset,S,thisDatasetIdx,datasetConfig,nSubplots,datasetNames)

% Colours for bars
positiveColor = '#31a354';
negativeColor = '#08519c';

% ytick configuration
yTickMax = datasetConfig(thisDatasetIdx).maxPercent;
yTickMin = -1*datasetConfig(thisDatasetIdx).maxPercent;
if (yTickMax > 30)
    yTickStep = 20;
    yTickFmt = '%.0f%%';
elseif (yTickMax > 20 && yTickMax <= 30)
    yTickStep = 10;
    yTickFmt = '%.0f%%';
elseif (yTickMax > 10 && yTickMax <= 20)
    yTickStep = 5;
    yTickFmt = '%.0f%%';
elseif (yTickMax <= 10 && yTickMax > 2)
    yTickStep = 2;
    yTickFmt = '%.0f%%';
elseif (yTickMax <= 2)
     yTickStep = 0.5;
     yTickFmt = '%.1f%%';
end

% Move subplots
ax(iDataset).pos = get(haxis(iDataset),'Position');
if (iDataset == 1 || iDataset == 2)
    ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.045;
elseif (iDataset == 3 || iDataset == 4)
    ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.015;
elseif (iDataset == 5 || iDataset == 6)
    ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.015;
elseif (iDataset == 7 || iDataset == 8)
    ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.045;  
end  
if (iDataset == 2 || iDataset == 4 || iDataset == 6 || iDataset == 8)
    ax(iDataset).pos(1) = ax(iDataset).pos(1) + 0.06;
else
    ax(iDataset).pos(1) = ax(iDataset).pos(1) + 0.005;
end
set(haxis(iDataset),'Position',ax(iDataset).pos)

% Extract variables
y = S(thisDatasetIdx).yearlyMean.sceneAverage;
x = unique(year(S(thisDatasetIdx).yearlyMean.time)); 

% Calculate annual mean
longTermMean = mean(y,'omitnan');
diffFromLongTermMean = 100*((y - longTermMean)/longTermMean);
idxPositiveYears = diffFromLongTermMean >= 0;
idxNegativeYears = diffFromLongTermMean < 0;

% Plot bars: annual deviation from the long-term mean
h = bar(haxis(iDataset),x,diag(diffFromLongTermMean),'stacked','BaseValue',0);
hold on
set(h(idxPositiveYears),'FaceColor',positiveColor) 
set(h(idxNegativeYears),'FaceColor',negativeColor)

% Plot fit line on top
time = datenum(x); % transform datetime into double
[p,str,mu] = polyfit(time(~isnan(diffFromLongTermMean)),diffFromLongTermMean(~isnan(diffFromLongTermMean)),1);
yDiffTrend = polyval(p,time,str,mu);
plot(haxis(iDataset),x,yDiffTrend,'Color','k','LineWidth',2);
hold off

% Axis arrangement
yLim = [yTickMin yTickMax];
yTicks = [yTickMin:yTickStep:yTickMax];
yTickLabels = arrayfun(@(y) sprintf(yTickFmt, y), yTicks, 'UniformOutput', false);

ylim(yLim)
yticks(yTicks)
yticklabels(yTickLabels)

grid on

firstYear = x(1);
lastYear = x(end);

%         if (firstYear == 1997)
%             xticks(datetime(1998,01,01):calyears(2):datetime(2022,12,31))
%             xticklabels({'1998','2000','2002','2004','2006','2008','2010','2012','2014','2016','2018','2020','2022'})
%         elseif (firstYear == 2016)
%             xticks(datetime(2016,01,01):calyears(1):datetime(2023,12,31))
%             xticklabels({'2016','2017','2018','2019','2020','2021','2022','2023'})
%         end
%         xtickangle(45)

% Add text annotation for the average
xt = min(xlim);
yt = max(ylim);

text(xt, yt,... % text position relative to axis
    {[strcat(num2str(firstYear,'%.0f'),'–',num2str(lastYear,'%.0f'),' avg:',num2str(longTermMean,'%.2f'),datasetConfig(thisDatasetIdx).unitsVar)]},...
    'FontSize',10,'Horiz','right','Vert','top');

yl = ylabel({[cell2mat(strcat({'Difference from the '}))],...
             [cell2mat(strcat(num2str(firstYear,'%.0f'),{'–'},num2str(lastYear,'%.0f'),'avg.'))]},'FontSize',10);
yl.Position(1) = yl.Position(1);

box on

title(datasetNames{thisDatasetIdx},'FontSize',12);

end % setAnomalyPlotFeatures

% *************************************************************************

end % plotAnomalyPlotsForDrivers