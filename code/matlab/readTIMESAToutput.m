clc
clear all

fullPathRootDir         = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/';
fullPathMainDir         = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/OCCCI_L3_data/';
fullPathTimeSeriesDir   = strcat(fullPathMainDir,'data_timeseries_nc/');
fullPathDataDir         = strcat(fullPathMainDir,'timesat_analysis/');
fullPathBinaryDataDir   = strcat(fullPathDataDir,'timesat_chla_bin/');
fullPathPlotsDir        = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/plots/phenology/';
filenameInfo            = 'CHLA_1998_2023.txt';
filenameTimeSeriesFit   = 'CHLA_1998_2023_fit.tts';
filenameTimeSeriesRaw   = 'CHLA_1998_2023_raw.tts';
filenameSeasonality     = 'CHLA_1998_2023_TS.tpa';
fullPathAreaStudyCoords = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/coords_boundarybox/bbox_50km';
fullPathEnduranceCoords = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/coords_endurance/endurance_2023';
addpath(genpath(fullPathMainDir))

%% Read original data and Endurance coords

% This was created by the script processSatelliteChla.m
load(strcat(fullPathRootDir,'continuousSatelliteChlaFiveDay4km.mat'),...
    'continuousSatelliteChlaFiveDay4km','timeVectorCompleteFiveDay4km',...
    'latSatellite','lonSatellite')

areaStudy = m_shaperead(fullPathAreaStudyCoords); % lat/lon coordinates
enduranceSite = m_shaperead(fullPathEnduranceCoords); % lat/lon coordinates

%% Read the .txt file with run information

% Open the text file for reading
txtFilePath = strcat(fullPathDataDir,filenameInfo);
infoID = fopen(txtFilePath,'r');
if infoID == -1
    error(['Could not open file: ', txtFilePath]);
end

% Initialize variables to store extracted information
startDate = '';
endDate = '';

% Iterate through the file
while ~feof(infoID)
    line = fgetl(infoID);
    if startsWith(line, 'Start date is')
        startDate = sscanf(line, 'Start date is %s');
    elseif startsWith(line, 'End date is')
        endDate = sscanf(line, 'End date is %s');
    end
end
fclose(infoID);

%% Read the file with fit time-series (*.tts)

binaryFilePath = strcat(fullPathDataDir,filenameTimeSeriesFit);
timeseriesFitID = fopen(binaryFilePath,'rb'); % 'rb' stands for read binary
if timeseriesFitID == -1
    error(['Could not open file: ', binaryFilePath]);
end

% Read header information
nyears = fread(timeseriesFitID, 1, 'int32');
nptperyear = fread(timeseriesFitID, 1, 'int32');
rowstart = fread(timeseriesFitID, 1, 'int32');
rowstop = fread(timeseriesFitID, 1, 'int32');
colstart = fread(timeseriesFitID, 1, 'int32');
colstop = fread(timeseriesFitID, 1, 'int32');

% Calculate the number of pixels
nRows = rowstop - rowstart + 1;
nCols = colstop - colstart + 1;

% Initialize the data and indices matrices
timeseriesFitIndxs = zeros(nRows * nCols, 2, 'int32');
timeseriesFitData  = zeros(nRows * nCols, nptperyear * nyears, 'single');

% Loop through each pixel and read the time-series data and indeces
for iRow = 1:nRows
    for iCol = 1:nCols
        timeseriesFitIndxs((iRow - 1) * nCols + iCol, :) = fread(timeseriesFitID, 2, 'int32');
        timeseriesFitData((iRow - 1) * nCols + iCol, :) = fread(timeseriesFitID, nptperyear * nyears, 'single');
    end
end
fclose(timeseriesFitID);

% Transform into double-precision values
timeseriesFitIndxs = double(timeseriesFitIndxs);
timeseriesFitData = double(timeseriesFitData);

%% Read the file with raw time-series (*.tts)

binaryFilePath = strcat(fullPathDataDir,filenameTimeSeriesRaw);
timeseriesRawID = fopen(binaryFilePath,'rb'); % 'rb' stands for read binary
if timeseriesRawID == -1
    error(['Could not open file: ', binaryFilePath]);
end

% Read header information
nyears = fread(timeseriesRawID, 1, 'int32');
nptperyear = fread(timeseriesRawID, 1, 'int32');
rowstart = fread(timeseriesRawID, 1, 'int32');
rowstop = fread(timeseriesRawID, 1, 'int32');
colstart = fread(timeseriesRawID, 1, 'int32');
colstop = fread(timeseriesRawID, 1, 'int32');

% Calculate the number of pixels
nRows = rowstop - rowstart + 1;
nCols = colstop - colstart + 1;

% Initialize the data and indices matrices
timeseriesRawIndxs = zeros(nRows * nCols, 2, 'int32');
timeseriesRawData  = zeros(nRows * nCols, nptperyear * nyears, 'single');

% Loop through each pixel and read the time-series data and indeces
for iRow = 1:nRows
    for iCol = 1:nCols
        timeseriesRawIndxs((iRow - 1) * nCols + iCol, :) = fread(timeseriesRawID, 2, 'int32');
        timeseriesRawData((iRow - 1) * nCols + iCol, :) = fread(timeseriesRawID, nptperyear * nyears, 'single');
    end
end
fclose(timeseriesRawID);

% Transform into double-precision values
timeseriesRawIndxs = double(timeseriesRawIndxs);
timeseriesRawData = double(timeseriesRawData);

%% Read the file with seasonality parameters (*.tpa)

binaryFilePath = strcat(fullPathDataDir,filenameSeasonality);
seasonalityID = fopen(binaryFilePath,'rb'); % 'rb' stands for read binary
if seasonalityID == -1
    error(['Could not open file: ', binaryFilePath]);
end

% Read header information
nyears = fread(seasonalityID, 1, 'int32');
nptperyear = fread(seasonalityID, 1, 'int32');
rowstart = fread(seasonalityID, 1, 'int32');
rowstop = fread(seasonalityID, 1, 'int32');
colstart = fread(seasonalityID, 1, 'int32');
colstop = fread(seasonalityID, 1, 'int32');

% Calculate the number of pixels
nRows = rowstop - rowstart + 1;
nCols = colstop - colstart + 1;

% Initialize the data cell array
seasonalityIndxs = zeros(nRows * nCols, 3, 'int32');
seasonalityData = cell(nRows, nCols);

% Loop through each pixel
for iRow = 1:nRows
    for iCol = 1:nCols
        
        % Read row, col, and n for the current pixel
        pixelInfo = fread(seasonalityID, 3, 'int32');
        seasonalityIndxs((iRow - 1) * nCols + iCol, :) = pixelInfo;
        n = pixelInfo(3);
        
        % Read seasonality parameters for each season
        seasonalityParams = fread(seasonalityID, [13, n], 'single'); % there are 13 phenology metrics
        seasonalityData{iRow, iCol} = struct('row', pixelInfo(1),...
                                             'col', pixelInfo(2),...
                                             'n', n,...
                                             'seasonality', seasonalityParams);
    end
end
fclose(seasonalityID);

% Extract data for first bloom period of every year and then for the second
% bloom period of every year
seasonalityDataFirstBloom = NaN(nRows,nCols,13,ceil(n/2)); 
seasonalityDataSecondBloom = NaN(nRows,nCols,13,ceil(n/2));

for iRow = 1:nRows
    for iCol = 1:nCols
        endSecondSeason = seasonalityData{iRow,iCol}.seasonality(2,2);
        if (endSecondSeason < 365) 
            seasonalityDataFirstBloom(iRow,iCol,:,:) =... 
                seasonalityData{iRow,iCol}.seasonality(:,1:2:n);
            seasonalityDataSecondBloom(iRow,iCol,:,1:ceil(n/2)-1) =... 
                seasonalityData{iRow,iCol}.seasonality(:,2:2:n-1);
        elseif (endSecondSeason >= 365)
            seasonalityDataFirstBloom(iRow,iCol,:,1:ceil(n/2)-1) =... 
                seasonalityData{iRow,iCol}.seasonality(:,2:2:n-1);
            seasonalityDataSecondBloom(iRow,iCol,:,:) =... 
                seasonalityData{iRow,iCol}.seasonality(:,1:2:n);
        end
        
    end
end

% Day metrics (day of bloom start, day of bloom end and day of mid season)
% need to be converted into day number from 1 to 365
dateVector = datetime(startDate):datetime(endDate);
for iMetric = [1,2,5]
    for iRow = 1:nRows
        for iCol = 1:nCols
            % First season
            data = squeeze(seasonalityDataFirstBloom(iRow,iCol,iMetric,:));
            indices = round(data);
            dayOfYear = day(dateVector(indices(indices > 0)),'dayofyear');
            data(data > 0) = dayOfYear;
            seasonalityDataFirstBloom(iRow,iCol,iMetric,:) = data;
            % Second season
            data = squeeze(seasonalityDataSecondBloom(iRow,iCol,iMetric,:));
            indices = round(data);
            dayOfYear = day(dateVector(indices(indices > 0)),'dayofyear');
            data(data > 0) = dayOfYear;
            seasonalityDataSecondBloom(iRow,iCol,iMetric,:) = data;
        end
    end    
end

%% Save

phenologyLon = lonSatellite;
phenologyLat = latSatellite;
phenologyTime = dateVector';

save(strcat(fullPathRootDir,'phenologyMetrics.mat'),...
    'seasonalityDataFirstBloom','seasonalityDataSecondBloom',...
    'phenologyTime','phenologyLat','phenologyLon')

%% Plot phenology metrics

startYear = year(datetime(startDate));
endYear = year(datetime(endDate));
yearLabels = num2cell(startYear:endYear);

metricLabels = {'Day of bloom start',... % 1
                'Day of bloom end',... % 2
                'Bloom duration (days)',... % 3
                'Base level (mg chla m^{-3})',... % 4
                'Day of mid season',... % 5
                'Peak intensity (mg chla m^{-3})',... % 6
                'Amplitude of the bloom (mg chla m^{-3})',... % 7
                'rate of increase at the beginning of the season',... % 8
                'rate of decrease at the end of the season',... % 9
                'large seasonal integral',... % 10
                'small seasonal integral',... % 11
                'value for the start of the season',... % 12
                'value for the end of the season'}; % 13

% Play with this to inspect variable limits by looking at the histogram of
% values
for iBloom = 1:2
    switch iBloom
        case 1
            thisBloomDataArray = seasonalityDataFirstBloom;
        case 2
            thisBloomDataArray = seasonalityDataSecondBloom;
    end
    for iMetric = 1:length(metricLabels)
        thisMetricAndYearDataArray = squeeze(thisBloomDataArray(:,:,iMetric,:));
        thisMetricAndYearDataArray(isnan(thisMetricAndYearDataArray)) = 0;
        [sortedArray,sortedIndices] = sort(reshape(thisMetricAndYearDataArray,[],1));
%         figure()
%         histogram(sortedArray,100)  %'BinLimits',[1e-3,2]
%         set(gca,'YScale','log')
    end
end
% AA = squeeze(thisMetricAndYearDataArray(:,:,1)); 
        
for iBloom = 1:2
    
    switch iBloom
        case 1
            nYears = ceil(n/2);
            thisBloomDataArray = seasonalityDataFirstBloom;
            axisLimits = [1,    50;...  % day bloom start
                          80,   200;... % day bloom end
                          30,   170;... % bloom duration
                          0.30, 1.0;... % base level
                          20,   120;... % day mid bloom
                          1.0,  4.0];   % bloom peak
        case 2
            nYears = ceil(n/2);
            thisBloomDataArray = seasonalityDataSecondBloom;
            axisLimits = [170,  280;...
                          290,  330;...
                          30,   170;...
                          0.30, 1.0;...
                          200,  290;...
                          1.0,  4.0];
    end
   
    for iMetric = 1:6 %length(metricLabels)
        
        figure()               
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.75],'Color','w') 
        haxis = zeros(5,5);
    
        for iYear = 1:(nYears-1) % DO NOT SHOW LAST YEAR (2023)
            
            haxis(iYear) = subaxis(5,5,iYear,'Spacing',0.015,'Padding',0.010,'Margin',0.04);
            ax(iYear).pos = get(haxis(iYear),'Position');
            
            % Move subplots
            if (iYear == 1 || iYear == 6 || iYear == 11 || iYear == 16 || iYear == 21 || iYear == 26)
                ax(iYear).pos(1) = ax(iYear).pos(1) + 0.025;
            elseif (iYear == 2 || iYear == 7 || iYear == 12 || iYear == 17 || iYear == 22 || iYear == 27)
                ax(iYear).pos(1) = ax(iYear).pos(1) + 0;
            elseif (iYear == 3 || iYear == 8 || iYear == 13 || iYear == 18 || iYear == 23 || iYear == 28)
                ax(iYear).pos(1) = ax(iYear).pos(1) - 0.025;
            elseif (iYear == 4 || iYear == 9 || iYear == 14 || iYear == 19 || iYear == 24 || iYear == 29)
                ax(iYear).pos(1) = ax(iYear).pos(1) - 0.050;  
            elseif (iYear == 5 || iYear == 10 || iYear == 15 || iYear == 20 || iYear == 25 || iYear == 30)
                ax(iYear).pos(1) = ax(iYear).pos(1) - 0.075;
            end  
            set(haxis(iYear),'Position',ax(iYear).pos) 
    
            thisMetricAndYearDataArray = squeeze(thisBloomDataArray(:,:,iMetric,iYear));
            thisMetricAndYearDataArray(thisMetricAndYearDataArray == 0) = NaN;
            m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
            m_pcolor(lonSatellite,latSatellite,thisMetricAndYearDataArray) 
            shading flat
            cmap = brewermap(1000,'*RdYlBu');
            colormap(cmap)
            caxis(axisLimits(iMetric,:))
            
%             if (iMetric == 1 || iMetric == 2 || iMetric == 5)
%                                 
%             % Define the tick positions representing days of the year (1 to 365)
%             tickPositions = get(colorbar, 'Ticks');
%             
%             % Calculate the corresponding month indices (1 to 12)
%             monthIndices = ceil(tickPositions / (365 / 12));
% 
%             % Calculate the corresponding day in each month
%             dayInMonth = mod(tickPositions - 1, 30) + 1; % Assuming each month has 30 days
% 
%             % Create labels representing the 1st day of each month
%             tickLabels = strcat(num2str(dayInMonth'), {' '}, datestr(datenum(0, monthIndices, 1), 'mmm'));
% 
%             % Set the colorbar tick labels to represent the 1st day of each month
%             set(colorbar, 'Ticks', tickPositions, 'TickLabels', tickLabels);
%             
%             end
            
            set(gca, 'color', 'w'); % or whatever color you want for the background
            
            % Show axis labels for border subplots
            if (iYear == 1 || iYear == 6 || iYear == 11 || iYear == 16)
                m_grid('linewi',1,'tickdir','out','FontSize',8,...
                    'xticklabel',[]);
            elseif (iYear == 22 || iYear == 23 || iYear == 24 || iYear == 25)
                m_grid('linewi',1,'tickdir','out','FontSize',8,...
                    'yticklabel',[]);
            elseif (iYear == 21)
                m_grid('linewi',1,'tickdir','out','FontSize',8)
            else
                m_grid('linewi',1,'tickdir','out','FontSize',8,...
                    'xticklabel',[],'yticklabel',[]);
            end

            m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
%             xlabel('Longitude');
%             ylabel('Latitude');
            title(yearLabels{iYear},'FontSize',9)
            box on
            hold on
        end
        
%         % Add a common title for the entire group of subplots
%         annotation('textbox', [0.1 0.9 0.8 0.1], 'String', metricLabels{iMetric}, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14);
%         
        % Colour bar settings
        cb = colorbar('Location','eastoutside');

        % Set tick marks to be middle of each range
%         dTk = diff(cb.Limits)/(2*length(cmap));
%         set(cb,'Ticks',[cb.Limits(1)+dTk:2*dTk:cb.Limits(2)-dTk],...
%             'TickLabels',({'1','2','3','4','5','6','7','8','9','10','11','12','13','14'}))
        cb.Position(1) = cb.Position(1) + 0.06;
        cb.Position(2) = cb.Position(2) + 0.07;
        cb.Position(3) = 0.025; % WIDTH
        cb.Position(4) = 0.70; % LENGTH
        cb.Label.String = metricLabels{iMetric}; 
        cb.FontSize = 8.5;

        set(gcf,'PaperPositionMode','auto')
        print(gcf,fullfile(fullPathPlotsDir,strcat('phenology_season_',num2str(iBloom),'_metric_',num2str(iMetric))),'-dpdf','-r0')
    
    end % iMetric
    hold on  
end % iBloom
hold off

%% Plot fits for this run

% Randomly sample the vector timeseriesFitIndxs four times
nSamples = 4; 
randomIndxs = randperm(height(timeseriesFitIndxs), nSamples);

% To focus the x axis
startDate = datetime('2015-01-01');
endDate = datetime('2018-12-31'); % datetime('2023-12-31');

figure()               
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.65 0.45],'Color','w') 
haxis = zeros(nSamples,1);

for iSample = 1:nSamples
    
    haxis(iSample) = subaxis(2,2,iSample,'Spacing',0.018,'Padding',0.020,'Margin',0.1);
    
    ax(iSample).pos = get(haxis(iSample),'Position');
    if (iSample == 1 || iSample == 2)
        ax(iSample).pos(2) = ax(iSample).pos(2) + 0.05;
    end
    ax(iSample).pos(1) = ax(iSample).pos(1) - 0.04;
    set(haxis(iSample),'Position',ax(iSample).pos) 
    
    plot(haxis(iSample),dateVector,timeseriesFitData(randomIndxs(iSample),:),...
        'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); hold on
    scatter(haxis(iSample),dateVector,timeseriesRawData(randomIndxs(iSample),:),...
        10,'b','filled'); hold off
    xlim([startDate endDate])
    currentYLim = ylim;
    ylim([0 currentYLim(2)]);
    grid on;
    ylabel('Chla (mg m^{-3})','FontSize',11)
    if (iSample == 3 || iSample == 4)
        xlabel('Time (date)','FontSize',11);
    end
    title(sprintf('Row=%u and Col=%u', timeseriesFitIndxs(randomIndxs(iSample),1), timeseriesFitIndxs(randomIndxs(iSample),2)))
    box on
end

lg = legend('Fit data (Savitzky-Golay)','Raw data (5-day)');  
lg.Orientation = 'vertical';
lg.Position(1) = 0.84; lg.Position(2) = 0.86;
lg.ItemTokenSize = [15,1];
lg.FontSize = 11;
set(lg,'Box','off')

set(gcf,'PaperPositionMode','manual')
orient(gcf,'landscape')
print(gcf,fullfile(fullPathPlotsDir,'plot_checkfit'),'-dpdf','-r0') 

%% Compare fits of the thre runs

% randomIndx = randperm(89*150, 1);

fullPathTimeSeriesFit1 = strcat(fullPathMainDir,'timesat_analysis_1/CHLA_1998_2022_fit.tts');
fullPathTimeSeriesRaw1 = strcat(fullPathMainDir,'timesat_analysis_1/CHLA_1998_2022_raw.tts');
fullPathTimeSeriesInfo1 = strcat(fullPathMainDir,'timesat_analysis_1/CHLA_1998_2022.txt');

fullPathTimeSeriesFit2 = strcat(fullPathMainDir,'timesat_analysis_2/CHLA_1998_2023_fit.tts');
fullPathTimeSeriesRaw2 = strcat(fullPathMainDir,'timesat_analysis_2/CHLA_1998_2023_raw.tts');
fullPathTimeSeriesInfo2 = strcat(fullPathMainDir,'timesat_analysis_2/CHLA_1998_2023.txt');

nFits = 2;

% Initialize empty arrays
timeseriesFitIndxs = [];
timeseriesFitData = [];
timeseriesRawData = [];
            
figure()               
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.75 0.45],'Color','w') 
haxis = zeros(nFits,1);

for iFit = 1:nFits
    
    switch iFit
        case 1
            binaryFilePathFit = fullPathTimeSeriesFit1;
            binaryFilePathRaw = fullPathTimeSeriesRaw1;
            txtFilePath = fullPathTimeSeriesInfo1;
        case 2
            binaryFilePathFit = fullPathTimeSeriesFit2;
            binaryFilePathRaw = fullPathTimeSeriesRaw2;
            txtFilePath = fullPathTimeSeriesInfo2;
    end
    
    % Read the info file to form the date vector
    infoID = fopen(txtFilePath,'r');

    startDate = '';
    endDate = '';

    while ~feof(infoID)
        line = fgetl(infoID);
        if startsWith(line, 'Start date is')
            startDate = sscanf(line, 'Start date is %s');
        elseif startsWith(line, 'End date is')
            endDate = sscanf(line, 'End date is %s');
        end
    end
    fclose(infoID);

    dateVector = datetime(startDate):datetime(endDate);
    
    % Read the .tts files to extract the time-series data
    for i = 1:2
        
        if (i == 1)
            timeseriesID = fopen(binaryFilePathFit,'rb'); 
        else
            timeseriesID = fopen(binaryFilePathRaw,'rb'); 
        end

        nyears = fread(timeseriesID, 1, 'int32');
        nptperyear = fread(timeseriesID, 1, 'int32');
        rowstart = fread(timeseriesID, 1, 'int32');
        rowstop = fread(timeseriesID, 1, 'int32');
        colstart = fread(timeseriesID, 1, 'int32');
        colstop = fread(timeseriesID, 1, 'int32');

        nRows = rowstop - rowstart + 1;
        nCols = colstop - colstart + 1;

        timeseriesIndxs = zeros(nRows * nCols, 2, 'int32');
        timeseriesData  = zeros(nRows * nCols, nptperyear * nyears, 'single');

        for iRow = 1:nRows
            for iCol = 1:nCols
                timeseriesIndxs((iRow - 1) * nCols + iCol, :) = fread(timeseriesID, 2, 'int32');
                timeseriesData((iRow - 1) * nCols + iCol, :) = fread(timeseriesID, nptperyear * nyears, 'single');
            end
        end
        fclose(timeseriesID);
        
        if (i == 1)
            timeseriesFitIndxs = double(timeseriesIndxs);
            timeseriesFitData = double(timeseriesData);
        else
            timeseriesRawData = double(timeseriesData);
        end
        
    end

    % Figure configuration
    
    haxis(iFit) = subaxis(nFits,1,iFit,'Spacing',0.018,'Padding',0.020,'Margin',0.1);
    
    ax(iFit).pos = get(haxis(iFit),'Position');
    if (iFit == 1)
        ax(iFit).pos(2) = ax(iFit).pos(2) + 0.08;
    elseif (iFit == 2)    
        ax(iFit).pos(2) = ax(iFit).pos(2) + 0.04;
    end
    set(haxis(iFit),'Position',ax(iFit).pos) 

    plot(haxis(iFit),dateVector,timeseriesFitData(randomIndx,:)); hold on
    scatter(haxis(iFit),dateVector,timeseriesRawData(randomIndx,:),10,'filled'); hold on
    xlim([datetime('1998-01-01') datetime('2000-12-31')])
    % ylim([0 1])
    grid on;
    ylabel('Chla (mg m^{-3})','FontSize',11)
    if (iFit == 3)
        xlabel('Time (date)','FontSize',11);
    end
    switch iFit
        case 1
            title(sprintf('Run with lower limit of 0.4 mg chla m^{-3}, pixel corresponding to row=%u and col=%u', timeseriesFitIndxs(randomIndx,1), timeseriesFitIndxs(randomIndx,2)))
        case 2
            title(sprintf('Run with lower limit of 0.2 mg chla m^{-3}, pixel corresponding to row=%u and col=%u', timeseriesFitIndxs(randomIndx,1), timeseriesFitIndxs(randomIndx,2)))
    end
    
    box on
    
end

lg = legend('Fit data','Raw data');  
lg.Orientation = 'vertical';
lg.Position(1) = 0.89; lg.Position(2) = 0.88;
lg.ItemTokenSize = [15,1];
lg.FontSize = 11;
set(lg,'Box','off')

set(gcf,'PaperPositionMode','manual')
orient(gcf,'landscape')
print(gcf,fullfile(fullPathPlotsDir,'plot_checkfit'),'-dpdf','-r0') 

