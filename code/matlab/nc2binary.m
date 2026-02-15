clc
clear all

fullPathMainDir         = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/OCCCI_L3_data/';
%fullPathMainDir         = '/home/anna/Documents/OCCCI_L3_data/';
fullPathDataDir         = strcat(fullPathMainDir,'timesat_analysis/');
fullPathBinaryDataDir   = strcat(fullPathDataDir,'timesat_chla_bin/');
fullPathBinaryWeightDir = strcat(fullPathDataDir,'timesat_weights_bin/');
fullPathToolsDir        = '/Users/Anna/LocalDocuments/Academic/Tools/';
fullPathPlotsDir        = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/plots/phenology/';
addpath(genpath(fullPathToolsDir))
addpath(genpath(fullPathMainDir))

% Box of 100 x 100 km around the Endurance CCS site
minLat = 53.758;
maxLat = 54.656;
minLon = 0.265;
maxLon = 1.801;

occciVarName.Latitude = 'lat';
occciVarName.Longitude = 'lon';
occciVarName.Time = 'time';

%% Read nc files

for iDataset = 1:2
    
    switch iDataset
        case 1
            outputFileName = 'data_timeseries_nc/occci_1km_chl_9710.nc';
        case 2
            outputFileName = 'data_timeseries_nc/occci_1km_chl_1023.nc';
%         case 3    
%             outputFileName = 'data_timeseries_nc/occci_1km_reflectance_9710.nc';
%         case 4
%             outputFileName = 'data_timeseries_nc/occci_1km_reflectance_1023.nc';
    end
    
    [D,lat,lon,time,colourVarNames,colourVarUnits] =... 
        extractVariables(strcat(fullPathMainDir,outputFileName),...
        occciVarName,minLat,maxLat,minLon,maxLon);
    
    % Concatenate arrays along the time dimension
    switch iDataset
        case 1
            D_chl_firsthalf = D;
            time_firsthalf = time;
        case 2
            D_chl_all = cat(3, D_chl_firsthalf, D);
            time_all = cat(1, time_firsthalf, time);
%         case 3    
%             D_reflectance_firsthalf = D;
%         case 4
%             D_reflectance_all = cat(3, D_reflectance_firsthalf, D);
    end
    
    if (iDataset == 1)
        D_varnames = colourVarNames;
        D_varunits = colourVarUnits;
    elseif (iDataset > 1 && rem(iDataset,2)) % if it's an odd number
        D_varnames = [D_varnames; colourVarNames];
        D_varunits = [D_varunits; colourVarUnits];
    end
 
end

% Inspect the chlorophyll array. For instance, how does the scene
% containing the maximum value looks like? The image shows that the pixel
% containing the maximum value is an outlier as the surrounding pixels have
% a mux lower value and are surrounded by clouds (artifact of chla
% algorithm)
%
% chla = D_chl_all(:,:,:,1);
% maxValue = max(chla,[],'all');
% maxLinearIndex = find(chla == maxValue, 1);
% 
% % Convert the linear index to subscripts (positions) in each dimension
% [maxIdx1, maxIdx2, maxIdx3] = ind2sub(size(chla), maxLinearIndex);
% 
% thisDataArray = chla(:,:,maxIdx3);
% thisDataArray(isnan(thisDataArray)) = -1; 
% figure()
% pcolor(lon,lat,thisDataArray)
% caxis([0, 100]);
% shading flat
% box on
% colorbar
% colormap('jet');
% cmap = colormap;
% cmap(1, :) = [1, 1, 1];  % set the FIRST color (-1 for NaNs) to white
% colormap(cmap);
% xlabel('Longitude');
% ylabel('Latitude');
  
%% Save into binary files

% Choose start and end date
startDate = datetime('1998-01-01');
endDate = datetime('2023-12-25');

% Find dates
iFirstDate = find(timeVectorAllDays == startDate); % '2012-07-15'
iLastDate = find(timeVectorAllDays == endDate); % '2023-02-18'
nImages = iLastDate - iFirstDate + 1; % the no. images per year must be an integer number -change dates until that condition is met
nYears = years(timeVectorAllDays(iLastDate) - timeVectorAllDays(iFirstDate));
nImagesPerYear = nImages/round(nYears,0);

% Check dimensions
[m,n] = size(squeeze(chla(:,:,1))); % we need to use "flipud(rot90(dataArray))" because that's what gets plotted in TIMESAT
nRowsPerFile = m; % 89 (were no. rows = no. latitude points)
nColsPerFile = n; % 150 (were no. columns = no. longitude points)

% Open .txt file for writing run information
strStartYear = datestr(timeVectorAllDays(iFirstDate),'yyyy');
strEndYear = datestr(timeVectorAllDays(iLastDate),'yyyy');
infoFileID = fopen(strcat(fullPathDataDir,'CHLA_',strStartYear,'_',strEndYear,'.txt'),'w');
fprintf(infoFileID,'Start date is %s\n',datestr(timeVectorAllDays(iFirstDate),'yyyy-mm-dd'));
fprintf(infoFileID,'End date is %s\n',datestr(timeVectorAllDays(iLastDate),'yyyy-mm-dd'));
fprintf(infoFileID,'The number of figures analysed is %d\n',nImages);
fprintf(infoFileID,'The number of years analysed is %f\n',nYears);
fprintf(infoFileID,'Each image has dimensions %d rows x %d columns\n',nRowsPerFile,nColsPerFile);
fclose(infoFileID);

% Check if the folder exists and, if so, delete it and create a new one
if exist(fullPathBinaryDataDir, 'dir') == 7
    rmdir(fullPathBinaryDataDir, 's');
end
mkdir(fullPathBinaryDataDir);

if exist(fullPathBinaryWeightDir, 'dir') == 7
    rmdir(fullPathBinaryWeightDir, 's');
end
mkdir(fullPathBinaryWeightDir);

for iDay = iFirstDate:iLastDate
    
    strDate = string(timeVectorAllDays(iDay));
    strDateWithoutDashes = replace(strDate, '-', '');
    
    % We need to apply the following transformation (flipud and rot90) for
    % TIMESAT to show us the data in the right spatial axis. This is
    % because MATLAB works column-wise and TIMESAT works row-wise
    dataArray = flipud(rot90(squeeze(chlaAllDays(:,:,iDay)))); 
%     dataArray(isnan(dataArray)) = minValue;
    
    % Weights array: set to '0.1' NaN values and set to '1.0' the rest
    weightsArray = dataArray;
    weightsArray(isnan(dataArray)) = 0.0;
    weightsArray(~isnan(dataArray)) = 1.0;
    weightsArray(dataArray > 0 & (dataArray <= 0.10 | dataArray >= 40)) = 0.10;
 
    % Binary files for data
    binaryFileName = strcat('chla_',strDateWithoutDashes,'.bin');
    dataID = fopen(strcat(fullPathBinaryDataDir,binaryFileName),'wb');
    fwrite(dataID,dataArray,'float32','ieee-le');
    fclose(dataID);
    
    % Binary files for weights (must be of the same type as the data)
    binaryFileName = strcat('flag_',strDateWithoutDashes,'.bin');
    weightID = fopen(strcat(fullPathBinaryWeightDir,binaryFileName),'wb');
    fwrite(weightID,weightsArray,'float32','ieee-le');
    fclose(weightID);
    
end

% Prepare file lists
for iFolder = 1:2
    switch iFolder
        case 1
            binaryDir = fullPathBinaryDataDir;
        case 2
            binaryDir = fullPathBinaryWeightDir;
    end
    fileList = dir(fullfile(binaryDir,'*.bin'));
    filePathsAndNames = fullfile(binaryDir,{fileList.name}');
    listContent = [num2str(nImages);filePathsAndNames]; % append on top of the list the total number of files
    fid = fopen(strcat(binaryDir,'listbinfiles.txt'),'w');
    fprintf(fid,'%s\n',listContent{:});
    fclose(fid);
end

% % Check that the binary files that we have created look alright
% binaryFilePath = strcat(fullPathBinaryDataDir,'chla_20200531.bin'); % 31st May 2020
% fid = fopen(binaryFilePath,'rb');  % 'rb' stands for read binary
% if fid == -1
%     error(['Could not open file: ', binaryFilePath]);
% end
% thisDataArray = fread(fid,[150,89],'float32'); 
% fclose(fid)
% 
% figure()
% pcolor(lon,lat,flipud(rot90(thisDataArray)))
% caxis([0.1, 4]);
% shading flat; box on; colorbar
% colormap('jet');
% cmap = colormap;
% cmap(1, :) = [1, 1, 1];  % set the FIRST color (-1 for NaNs) to white
% colormap(cmap);
% xlabel('Longitude');
% ylabel('Latitude');
% title('Chla 31st May 2020')
% % set(gcf,'PaperPositionMode','auto')
% % print(gcf,fullfile(fullPathDataDir,'plot_20000531_chl'),'-dpdf','-r0')  

%% Clasify data by season

chla = chlaAllDays(:,:,iFirstDate:iLastDate);
chla(isnan(chla)) = 0;

% Classify data by season
dates = timeVectorAllDays(iFirstDate:iLastDate);

% Define the start and end dates for each season
seasons = {'Winter','Spring','Summer','Autumn'};
seasonStartEnd = {
    datetime('1997-12-21'), datetime('1998-03-20');
    datetime('1998-03-21'), datetime('1998-06-20');
    datetime('1998-06-21'), datetime('1998-09-20');
    datetime('1998-09-21'), datetime('1998-12-20');
};

% Initialize cell arrays to store data for each season
winterData = cell(1, endDate.Year - startDate.Year + 1);
springData = cell(1, endDate.Year - startDate.Year + 1);
summerData = cell(1, endDate.Year - startDate.Year + 1);
autumnData = cell(1, endDate.Year - startDate.Year + 1);

for thisYear = startDate.Year:endDate.Year
    for iSeason = 1:length(seasons)
        
        if ~strcmp(seasons{iSeason},'Winter') 
            seasonStart = datetime(thisYear, seasonStartEnd{iSeason, 1}.Month, seasonStartEnd{iSeason, 1}.Day);
        elseif strcmp(seasons{iSeason},'Winter') % in the winter, start the season in the year before
            seasonStart = datetime(thisYear-1, seasonStartEnd{iSeason, 1}.Month, seasonStartEnd{iSeason, 1}.Day);
        end
        seasonEnd = datetime(thisYear, seasonStartEnd{iSeason, 2}.Month, seasonStartEnd{iSeason, 2}.Day);
        
        % Extract data for the current season and year
        currYearAndSeason = chla(:,:,dates >= seasonStart & dates <= seasonEnd);
        
        % Store the data in the respective cell arrays
        switch seasons{iSeason}
            case 'Winter'
                winterData{thisYear - startDate.Year + 1} = currYearAndSeason;
            case 'Spring'
                springData{thisYear - startDate.Year + 1} = currYearAndSeason;
            case 'Summer'
                summerData{thisYear - startDate.Year + 1} = currYearAndSeason;
            case 'Autumn'
                autumnData{thisYear - startDate.Year + 1} = currYearAndSeason;
        end
    end
end

%% Plot histograms

% Histogram with all the data

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.29],'Color','w') 

% Main plot
histogram(chla,100,'BinLimits',[1e-2,100],'FaceColor', [0.5, 0.5, 0.5])
set(gca,'YScale','log');
grid on;
mainAxes = gca;
mainAxes.XGrid = 'on';  % 'on' for main grid, 'off' for no grid
mainAxes.YGrid = 'on';  % 'on' for main grid, 'off' for no grid
mainAxes.MinorGridAlpha = 0;  % Set alpha (transparency) to 0 for secondary grid lines
yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]);
yticklabels({'1', '10', '10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}'});
xlabel('Chlorophyll a (mg m^{-3})');
ylabel('Frequency');
title('Distribution of satellite observations from the OC-CCI 1 km product');

% Inset plot
insetAxes = axes('Position', [0.4, 0.43, 0.47, 0.45]); % define the position of the inset plot
histogram(insetAxes,chla,100,'BinLimits',[1e-2,1],'FaceColor', [0.5, 0.5, 0.5])
set(insetAxes,'YScale','log');
grid on;
insetAxes.XGrid = 'on';  % 'on' for main grid, 'off' for no grid
insetAxes.YGrid = 'on';  % 'on' for main grid, 'off' for no grid
insetAxes.MinorGridAlpha = 0;  % Set alpha (transparency) to 0 for secondary grid lines
insetAxes.XAxis.FontSize = 8;
insetAxes.YAxis.FontSize = 8;
yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6]);
yticklabels({'1','10','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'});
xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
ylabel('Frequency','FontSize',8)

set(gcf, 'PaperPositionMode', 'auto')
print(gcf,fullfile(fullPathPlotsDir,'histogram_satellite_all'),'-dpdf','-r0')

% Histogram by season

colourScheme = brewermap(length(seasons),'*Spectral'); % flip upside down to have yellow colour for spring and blue for winter

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
        case 'Spring'
            seasonData = springData;
        case 'Summer'
            seasonData = summerData;
        case 'Autumn'
            seasonData = autumnData;
        case 'Winter'
            seasonData = winterData;
    end
    
    % Reshape cell array data into one column
    concatenatedData = cat(3, seasonData{:});
    columnData = reshape(concatenatedData, [], 1);
    
    % Main plot
    
    % Plot all data
    histogram(haxis(iSeason),chla,100,'BinLimits',[1e-2,100],'FaceColor',[0.7,0.7,0.7],'EdgeColor','none')
    hold on
    % Plot seasonal data
    histogram(haxis(iSeason),columnData,100,'BinLimits',[1e-2,100],'FaceColor',colourScheme(iSeason,:))
    hold on
    ylim([1 2e7])
    set(gca,'YScale','log');
    grid on;
    mainSeasonalAxes = gca;
    mainSeasonalAxes.XGrid = 'on';  % 'on' for main grid, 'off' for no grid
    mainSeasonalAxes.YGrid = 'on';  % 'on' for main grid, 'off' for no grid
    mainSeasonalAxes.MinorGridAlpha = 0;  % Set alpha (transparency) to 0 for secondary grid lines
    yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]);
    yticklabels({'1','10','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}'});
    
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
                                          0.49*ax(iSeason).pos(4)]); % define the position of the inset plot
    
    histogram(insetSeasonalAxes,columnData,100,'BinLimits',[1e-2,1],'FaceColor',colourScheme(iSeason,:))
    hold on
    ylim([1 2e5])
    set(insetSeasonalAxes,'YScale','log');
    grid on;
    insetSeasonalAxes.XGrid = 'on';  % 'on' for main grid, 'off' for no grid
    insetSeasonalAxes.YGrid = 'on';  % 'on' for main grid, 'off' for no grid
    insetSeasonalAxes.MinorGridAlpha = 0;  % Set alpha (transparency) to 0 for secondary grid lines
    insetSeasonalAxes.XAxis.FontSize = 8;
    insetSeasonalAxes.YAxis.FontSize = 8;
    yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6]);
    yticklabels({'1','10','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'});
    xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
    ylabel('Frequency','FontSize',8)   
end
hold off

set(gcf, 'PaperPositionMode', 'auto')
print(gcf,fullfile(fullPathPlotsDir,'histogram_satellite_byseason'),'-dpdf','-r0')
