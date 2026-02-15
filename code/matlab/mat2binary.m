clc
clear all

fullPathMainDir         = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/OCCCI_L3_data/';
%fullPathMainDir         = '/home/anna/Documents/OCCCI_L3_data/';
fullPathTimeSeriesDir   = strcat(fullPathMainDir,'data_timeseries_nc/');
fullPathDataDir         = strcat(fullPathMainDir,'timesat_analysis/');
fullPathBinaryDataDir   = strcat(fullPathDataDir,'timesat_chla_bin/');
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

% Check if the folder exists and, if so, delete it and create a new one
if exist(fullPathDataDir, 'dir') == 7
    rmdir(fullPathDataDir, 's');
end
mkdir(fullPathDataDir)

%% Load time-series data product

% This was created by the script processSatelliteChla.m
load(strcat(fullPathTimeSeriesDir,'continuousSatelliteChlaFiveDay4km.mat'),...
    'continuousSatelliteChlaFiveDay4km','timeVectorCompleteFiveDay4km',...
    'latSatellite','lonSatellite')

% Inspect the chlorophyll array. For instance, how does the scene
% containing the maximum value looks like? 
chla = continuousSatelliteChlaFiveDay4km(:,:,:);
maxValue = max(chla,[],'all');
maxLinearIndex = find(chla == maxValue, 1);

% Convert the linear index to subscripts (positions) in each dimension
[maxIdx1, maxIdx2, maxIdx3] = ind2sub(size(chla), maxLinearIndex);
thisDataArray = chla(:,:,maxIdx3);
thisDataArray(isnan(thisDataArray)) = -1; 
figure()
pcolor(lonSatellite,latSatellite,thisDataArray)
caxis([0, 10]);
shading flat
box on
colorbar
colormap('jet');
cmap = colormap;
cmap(1, :) = [1, 1, 1];  % set the FIRST color (-1 for NaNs) to white
colormap(cmap);
xlabel('Longitude');
ylabel('Latitude');
  
%% Save into binary files

% Choose start and end date
startDate = datetime('1998-01-01');
endDate = datetime('2023-12-25');
timeVectorDaily = (startDate:days(1):endDate)';
timeVectorDaily.Format = 'yyyy-MM-dd';

% Find dates
iFirstDate = find(timeVectorDaily == startDate); % '2012-07-15'
iLastDate = find(timeVectorDaily == endDate); % '2023-02-18'
nImages = iLastDate - iFirstDate + 1; % the no. images per year must be an integer number -change dates until that condition is met
nYears = years(timeVectorDaily(iLastDate) - timeVectorDaily(iFirstDate));
nImagesPerYear = nImages/round(nYears,0);

% Check dimensions
[m,n] = size(squeeze(chla(:,:,1))); % we need to use "flipud(rot90(dataArray))" because that's what gets plotted in TIMESAT
nRowsPerFile = m; % 89 (where no. rows = no. latitude points)
nColsPerFile = n; % 150 (where no. columns = no. longitude points)

% Open .txt file for writing run information
strStartYear = datestr(timeVectorDaily(iFirstDate),'yyyy');
strEndYear = datestr(timeVectorDaily(iLastDate),'yyyy');
infoFileID = fopen(strcat(fullPathDataDir,'CHLA_',strStartYear,'_',strEndYear,'.txt'),'w');
fprintf(infoFileID,'Start date is %s\n',datestr(timeVectorDaily(iFirstDate),'yyyy-mm-dd'));
fprintf(infoFileID,'End date is %s\n',datestr(timeVectorDaily(iLastDate),'yyyy-mm-dd'));
fprintf(infoFileID,'The number of figures analysed is %d\n',nImages);
fprintf(infoFileID,'The number of years analysed is %f\n',nYears);
fprintf(infoFileID,'Each image has dimensions %d rows x %d columns\n',nRowsPerFile,nColsPerFile);
fclose(infoFileID);

% Check if the folder exists and, if so, delete it and create a new one
if exist(fullPathBinaryDataDir, 'dir') == 7
    rmdir(fullPathBinaryDataDir, 's');
end
mkdir(fullPathBinaryDataDir);

iDayFiveDay = 1;
for iDayDaily = iFirstDate:iLastDate
    
    strDate = string(timeVectorDaily(iDayDaily));
    strDateWithoutDashes = replace(strDate, '-', '');
    
    % We need to apply the following transformation (flipud and rot90) for
    % TIMESAT to show us the data in the right spatial axis. This is
    % because MATLAB works column-wise and TIMESAT works row-wise
    
    if (~ismember(timeVectorDaily(iDayDaily), timeVectorCompleteFiveDay4km))
        dataArray = NaN(nRowsPerFile,nColsPerFile);
    else
        dataArray = flipud(rot90(squeeze(continuousSatelliteChlaFiveDay4km(:,:,iDayFiveDay)))); 
        iDayFiveDay = iDayFiveDay + 1;
    end
    
    % Weights array: set to '0.1' NaN values and set to '1.0' the rest
%     weightsArray = dataArray;
%     weightsArray(isnan(dataArray)) = 0.0;
%     weightsArray(~isnan(dataArray)) = 1.0;
%     weightsArray(dataArray > 0 & (dataArray <= 0.10 | dataArray >= 40)) = 0.10;
 
    % Binary files for data
    binaryFileName = strcat('chla_',strDateWithoutDashes,'.bin');
    dataID = fopen(strcat(fullPathBinaryDataDir,binaryFileName),'wb');
    fwrite(dataID,dataArray,'float32','ieee-le');
    fclose(dataID);
    
    % Binary files for weights (must be of the same type as the data)
%     binaryFileName = strcat('flag_',strDateWithoutDashes,'.bin');
%     weightID = fopen(strcat(fullPathBinaryWeightDir,binaryFileName),'wb');
%     fwrite(weightID,weightsArray,'float32','ieee-le');
%     fclose(weightID);
    
end

% Prepare file lists
for iFolder = 1 %:2
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

% Check that the binary files that we have created look alright
binaryFilePath = strcat(fullPathBinaryDataDir,'chla_20200530.bin'); % 30th May 2020
fid = fopen(binaryFilePath,'rb');  % 'rb' stands for read binary
if fid == -1
    error(['Could not open file: ', binaryFilePath]);
end
thisDataArray = fread(fid,[nColsPerFile,nRowsPerFile],'float32'); 
fclose(fid)

figure()
pcolor(lonSatellite,latSatellite,flipud(rot90(thisDataArray)))
caxis([0.1, 4]);
shading flat; box on; colorbar
colormap('jet');
cmap = colormap;
cmap(1, :) = [1, 1, 1];  % set the FIRST color (-1 for NaNs) to white
colormap(cmap);
xlabel('Longitude');
ylabel('Latitude');
title('Chla 30th May 2020')
set(gcf,'PaperPositionMode','auto')
print(gcf,fullfile(fullPathDataDir,'plot_20000530_chl'),'-dpdf','-r0')  
