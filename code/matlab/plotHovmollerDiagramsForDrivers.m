function plotHovmollerDiagramsForDrivers(AScmems,ASbicep,ASnasa,ASoccci,datasetConfig)

% PLOTHOVMOLLERDIAGRAMSFORDRIVERS Plot Hovmoller diagrams for driver data
%
%   INPUT: 
%       AScmems           - table with time-series data from CMEMS for our area of study
%       ASbicep
%
%   OUTPUT:
%       S                 - data structure with data extracted from selected datasets
%       datasetLabel      - string containing selected product labels for plotting
%       cmemsDatasetNames - string containing selected product names
% 
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 6 May 2024
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

%% Get data into the right format

for iDataset = 1:numel(datasetNames)
    
    time = S(iDataset).time;
    
    % Convert the datetime vector to just dates (with time set to midnight)
    dateOnlyVector = dateshift(time, 'start', 'day');
    
    % Evaluate if the data is daily or not. If it is daily, we have to
    % handle leap years by removing Feb 29th, so that the year has exactly
    % 365 days and not 366
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
    
    if (strcmp(frequency,'daily regular'))
        % Make datetimes corresponding to data
        t = dateOnlyVector(1):dateOnlyVector(end);
        % Find all the leap days 
        idl = month(t)==2 & day(t)==29; 
        % Remove leap days
        timeSeriesData = S(iDataset).sceneAverage;
        timeSeriesTime = dateOnlyVector;
        timeSeriesData = timeSeriesData(~idl);
        timeSeriesTime = timeSeriesTime(~idl);
    else
        timeSeriesData = S(iDataset).sceneAverage;
        timeSeriesTime = dateOnlyVector;
    end
    
    % Prepare data to build Hovmoller array
    TS = table(timeSeriesData,timeSeriesTime,'VariableNames',{'data','date'});
    TS.year = year(TS.date);
    TS.month = month(TS.date);
    yearVector = min(TS.year):1:max(TS.year);
    
    % Initialise Hovmoller array
    if (strcmp(frequency,'daily regular'))
        dataHovmollerArray = NaN(365,numel(yearVector));
        for iYear = 1:numel(yearVector)
            thisYearData = TS.data(TS.year == yearVector(iYear));
            % If this year is not complete (< 365 days of data)
            if (numel(thisYearData) < 365)
                dataHovmollerArray(1:numel(thisYearData),iYear) = thisYearData;
            else
                dataHovmollerArray(:,iYear) = thisYearData;
            end
        end
    elseif (strcmp(frequency,'monthly'))
        dataHovmollerArray = NaN(12,numel(yearVector));
        for iYear = 1:numel(yearVector)
            thisYearData = TS.data(TS.year == yearVector(iYear));
            dataHovmollerArray(:,iYear) = thisYearData;
        end
    elseif (strcmp(frequency,'daily irregular'))
        dataHovmollerArray = NaN(365,numel(yearVector));
        for iYear = 1:numel(yearVector)
            thisYearData = TS.data(TS.year == yearVector(iYear));
            thisYearDates = TS.date(TS.year == yearVector(iYear));
            % Calculate the day of the year for each date in the vector
            idxDayOfYear = day(thisYearDates, 'dayofyear');
            dataHovmollerArray(idxDayOfYear,iYear) = thisYearData;
        end
    end

    S(iDataset).hovmoller = dataHovmollerArray;
    
end

%% Plot features

colourHovmoller = [[1, 1, 1]; brewermap(100,'*YlGnBu')]; % append a row of white to colour values outside lower limit
monthNames = {'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'};

% To add month labels in the correct day position 
daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
cumulativeDays = cumsum(daysInMonth); % calculate the cumulative sum of days in each month
monthStartPositions = [1, cumulativeDays + 1]; % add a starting point (0) for the beginning of the year

%% Plot CMEMS Hovmoller diagrams

% We will plot in two rounds, with 8 subplots in each

nSubplots = 8;

for iRound = 1:2
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.70],'Color','w')
    haxis = zeros(nSubplots,1);
    
    if (iRound == 1)
        idxsDatasets = idxsCmemsDatasetsFirstRound; 
    else
        idxsDatasets = idxsCmemsDatasetsSecondRound;
    end 

    for iDataset = 1:nSubplots
        
        thisDatasetIdx = idxsDatasets(iDataset); 

        haxis(iDataset) = subaxis(4,2,iDataset,'Spacing',0.015,'Padding',0.010,'Margin',0.06);
        ax(iDataset).pos = get(haxis(iDataset),'Position');
        
        % Move subplots
        if (iDataset == 1 || iDataset == 2)
            ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.040;
        elseif (iDataset == 3 || iDataset == 4)
            ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.015;
        elseif (iDataset == 5 || iDataset == 6)
            ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.010;
        elseif (iDataset == 7 || iDataset == 8)
            ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.035;  
        end  
        if (iDataset == 2 || iDataset == 4 || iDataset == 6 || iDataset == 8)
            ax(iDataset).pos(1) = ax(iDataset).pos(1) + 0.030;
        else
            ax(iDataset).pos(1) = ax(iDataset).pos(1) - 0.025;
        end
        set(haxis(iDataset),'Position',ax(iDataset).pos) 
        
        hovmollerArray = S(thisDatasetIdx).hovmoller;
        yearVector = unique(year(S(thisDatasetIdx).time));
    
        imagesc(yearVector,1:365,hovmollerArray); % X-axis is years, Y-axis is days in year
        colormap(colourHovmoller)
        minVar = datasetConfig(thisDatasetIdx).minVar;
        maxVar = datasetConfig(thisDatasetIdx).maxVar;
        caxis([minVar maxVar])

        cb = colorbar;
        cb.Label.String = datasetConfig(thisDatasetIdx).unitsVar;
        cb.FontSize = 9;
        cb.Ruler.TickLabelFormat = datasetConfig(thisDatasetIdx).tickLabelFmt;
        
        % Customize the y-axis
        yticks(monthStartPositions); % Set the y-axis tick positions to the start of each month
        yticklabels(monthNames); % Set the y-axis tick labels to month names
        
        title(datasetNames{thisDatasetIdx},'FontSize',12);

    end % iDataset
    
    exportgraphics(gcf,fullfile('.','figures',strcat('hovmoller_cmems_',num2str(iRound),'.png')),'Resolution',600)

end % iRound

%% Plot BICEP Hovmoller diagrams

% We will plot 6 subplots

nSubplots = 6;

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.60],'Color','w')
haxis = zeros(nSubplots,1);

for iDataset = 1:nSubplots

    thisDatasetIdx = idxsBicepDatasets(iDataset); 

    haxis(iDataset) = subaxis(3,2,iDataset,'Spacing',0.015,'Padding',0.010,'Margin',0.06);
    ax(iDataset).pos = get(haxis(iDataset),'Position');

    % Move subplots
    if (iDataset == 1 || iDataset == 2)
        ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.040;
    elseif (iDataset == 3 || iDataset == 4)
        ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.015;
    elseif (iDataset == 5 || iDataset == 6)
        ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.010;
    end  
    if (iDataset == 2 || iDataset == 4 || iDataset == 6)
        ax(iDataset).pos(1) = ax(iDataset).pos(1) + 0.030;
    else
        ax(iDataset).pos(1) = ax(iDataset).pos(1) - 0.025;
    end
    set(haxis(iDataset),'Position',ax(iDataset).pos) 
        
    hovmollerArray = S(thisDatasetIdx).hovmoller;
    yearVector = unique(year(S(thisDatasetIdx).time));

    imagesc(yearVector,1:12,hovmollerArray);
    colormap(colourHovmoller)
    minVar = datasetConfig(thisDatasetIdx).minVar;
    maxVar = datasetConfig(thisDatasetIdx).maxVar;
    caxis([minVar maxVar])

    cb = colorbar;
    cb.Label.String = datasetConfig(thisDatasetIdx).unitsVar;
    cb.FontSize = 9;
    cb.Ruler.TickLabelFormat = datasetConfig(thisDatasetIdx).tickLabelFmt;

    % Customize the y-axis
    yticks(1:12); 
    yticklabels(monthNames);

    title(datasetNames{thisDatasetIdx},'FontSize',12);

end % iDataset

exportgraphics(gcf,fullfile('.','figures','hovmoller_bicep.png'),'Resolution',600)

%% Plot chla Hovmoller diagrams

% We will plot 8 subplots

nSubplots = 8;

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.70],'Color','w')
haxis = zeros(nSubplots,1);

for iDataset = 1:nSubplots

    thisDatasetIdx = idxsChlaDatasets(iDataset); 

    haxis(iDataset) = subaxis(4,2,iDataset,'Spacing',0.015,'Padding',0.010,'Margin',0.06);
    ax(iDataset).pos = get(haxis(iDataset),'Position');

    % Move subplots
    if (iDataset == 1 || iDataset == 2)
        ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.040;
    elseif (iDataset == 3 || iDataset == 4)
        ax(iDataset).pos(2) = ax(iDataset).pos(2) + 0.015;
    elseif (iDataset == 5 || iDataset == 6)
        ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.010;
    elseif (iDataset == 7 || iDataset == 8)
        ax(iDataset).pos(2) = ax(iDataset).pos(2) - 0.035;  
    end  
    if (iDataset == 2 || iDataset == 4 || iDataset == 6 || iDataset == 8)
        ax(iDataset).pos(1) = ax(iDataset).pos(1) + 0.030;
    else
        ax(iDataset).pos(1) = ax(iDataset).pos(1) - 0.025;
    end
    set(haxis(iDataset),'Position',ax(iDataset).pos) 
        
    hovmollerArray = S(thisDatasetIdx).hovmoller;
    yearVector = unique(year(S(thisDatasetIdx).time));

    imagesc(yearVector,1:365,log10(hovmollerArray));
    colormap(colourHovmoller)
    minVar = log10(0.2); %log10(chlaDatasetConfig(iDataset).minVar)
    maxVar = log10(10); %log10(chlaDatasetConfig(iDataset).maxVar)
    caxis([minVar maxVar])

    cb = colorbar;
    cb.Label.String = datasetConfig(thisDatasetIdx).unitsVar;
    cb.FontSize = 9;
    set(cb,'ytick',log10([0.3 0.5 1 2 3 5 10]),'yticklabel',[0.3 0.5 1 2 3 5 10]);
    cb.Ruler.TickLabelFormat = datasetConfig(thisDatasetIdx).tickLabelFmt;
  
    % Customize the y-axis
    yticks(monthStartPositions);
    yticklabels(monthNames);

    title(datasetNames{thisDatasetIdx},'FontSize',12);

end % iDataset

exportgraphics(gcf,fullfile('.','figures','hovmoller_chla.png'),'Resolution',600)

end % plotHovmollerDiagramsForDrivers