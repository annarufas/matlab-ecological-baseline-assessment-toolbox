
% This script generates an interpolated and filled chlorophyll a product
% for the Endurance site. Run it and the output will be saved in fullPathMainDir
% as continuousSatelliteChlaFiveDay4km.mat. As it runs, you have the
% option to turn on/off the possibility of plotting

clc
clear all

fullPathMainDir         = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/';
fullPathDataOCCCIdir    = strcat(fullPathMainDir,'OCCCI_L3_data/data_timeseries_nc/');
fullPathDataCMEMSdir    = strcat(fullPathMainDir,'CMEMS_data/data_timeseries_CMC_nc/');
fullPathPlotsDir        = strcat(fullPathMainDir,'plots/analysisOC/processChlaProduct/');
fullPathAreaStudyCoords = strcat(fullPathMainDir,'coords_boundarybox/bbox_50km');
fullPathEnduranceCoords = strcat(fullPathMainDir,'coords_endurance/endurance_2023');

addpath(genpath('/Users/Anna/LocalDocuments/Academic/Tools/')) 
addpath(genpath(strcat(fullPathMainDir,'CMEMS_data/')))
addpath(genpath(strcat(fullPathMainDir,'OCCCI_L3_data/'))) 

areaStudy = m_shaperead(fullPathAreaStudyCoords); % lat/lon coordinates
enduranceSite = m_shaperead(fullPathEnduranceCoords); % lat/lon coordinates

%% Choices

isCMEMSdataReady = 1; % 0=no, 1=yes
isOCCCIdataReady = 1;

isPlotSatelliteProductInterComparison = 1;
isPlotScenePercentageMissingData = 1;
isPlotSatelliteHistogram = 1;
isPlotHighValuesMask = 1;
isPlotInterpolation = 1;
isPlotDistributionPersistentMissingData = 1;
isPlotEuphoticLayer = 1;
isPlotCMEMSregridded = 1;
isPlotFilledProduct = 1;
isPlotComparisonWithInsituData = 1;

%% Load ocean colour data

if (isCMEMSdataReady == 0)
    disp("Reading CMEMS products...")
    ncreadCMEMSdata(fullPathDataCMEMSdir)
end
if (isOCCCIdataReady == 0)
    disp("Reading OC-CCI products (1-day, 1 km & 1-day, 4 km & 8-day, 4 km & 5-day, 1 km...")
    ncreadTimeseriesOCCCIdata(areaStudy,fullPathDataOCCCIdir)
end

load(strcat(fullPathDataCMEMSdir,'OC_cmems_bbox50km.mat'))
load(strcat(fullPathDataOCCCIdir,'OC_occci_bbox50km.mat'))

%% Process satellite data products

% Get rid of the year 1997 (lots of gaps), calculate season and calculate 
% the distribution of missing data

occciDatasets = {OC_occci.dataset_name};

startDateSatellite = datetime('1998-01-01');
fracMissingPixelToConsiderInvalidImage = 98; % percentage (%) figure

for iDataset = 1:length(occciDatasets)
    
    switch occciDatasets{iDataset} 

        case 'occci_1km_1day_chl'
            posIdProd = find(strcmp({OC_occci.dataset_name}, occciDatasets{iDataset}));
            iStartDate = find(OC_occci(posIdProd).time == startDateSatellite);
            chlaSatellite = OC_occci(posIdProd).dataset(:,:,iStartDate:end,1); % chla is in position 1
            timeVector = OC_occci(posIdProd).time(iStartDate:end);
            timeVectorComplete = (timeVector(1):timeVector(end))';
            
            % Find the missing days
            missingDays = setdiff(timeVectorComplete, timeVector);
            missingDayIdxs = find(ismember(timeVectorComplete, missingDays));
            % For those missing days, generate an array of NaN values
            [nRows,nCols,~] = size(chlaSatellite);
            chlaSatelliteComplete = NaN(nRows,nCols,length(timeVectorComplete));
            iImageWithData = 1;
            for iImage = 1:length(timeVectorComplete)
                if (ismember(iImage,missingDayIdxs))
                    chlaSatelliteComplete(:,:,iImage) = NaN(nRows,nCols);
                else
                    chlaSatelliteComplete(:,:,iImage) = chlaSatellite(:,:,iImageWithData);
                    iImageWithData = iImageWithData + 1;
                end
            end
            
        case 'occci_4km_5day_chl'
            posIdProd = find(strcmp({OC_occci.dataset_name}, occciDatasets{iDataset}));
            iStartDate = find(OC_occci(posIdProd).time == startDateSatellite);
            chlaSatelliteComplete = OC_occci(posIdProd).dataset(:,:,iStartDate:end,1); % chla is in position 1
            timeVectorComplete = OC_occci(posIdProd).time(iStartDate:end);

            % There are no missing days in the 5-day dataset, no need of any
            % further processing
            
%             % Initialise the complete time vector
%             timeVector = OC_occci(posIdProd).time(iStartDate:end);
%             timeVectorComplete = timeVector(1);
%             % Loop to generate the time vector
%             while timeVectorComplete(end) + 5 <= timeVector(end)
%                 newDate = timeVectorComplete(end) + 5;
%                 % Check if the new date is in the same year
%                 if year(newDate) ~= year(timeVectorComplete(end))
%                     % Adjust the date to the 1st of January of the next year
%                     newDate = datetime(['01/01/' num2str(year(newDate))], 'InputFormat', 'dd/MM/yyyy');
%                 end
%                 timeVectorComplete = [timeVectorComplete; newDate];
%             end

       case 'occci_4km_8day_chl'
            posIdProd = find(strcmp({OC_occci.dataset_name}, occciDatasets{iDataset}));
            iStartDate = find(OC_occci(posIdProd).time == startDateSatellite);
            chlaSatelliteComplete = OC_occci(posIdProd).dataset(:,:,iStartDate:end,1); % chla is in position 1
            timeVectorComplete = OC_occci(posIdProd).time(iStartDate:end);
            
            % There are no missing days in the 8-day dataset, no need of any
            % further processing
 
    end
    
    timeVectorComplete.Format = 'yyyy-MM-dd';
    [nRows,nCols,nExpectedImages] = size(chlaSatelliteComplete);
    nPixelsInImage = nRows*nCols;
    
    % Function to calculate missing pixel statistics
    [fracMissingPixels,isValidImage,fracMissingPixelsByMonth,...
        imageDistribByMonthAndValidity,M] = calculateMissingPixelsStatistics(...
        chlaSatelliteComplete,timeVectorComplete,...
        nExpectedImages,nPixelsInImage,fracMissingPixelToConsiderInvalidImage);

    % Add some return arrays
    switch occciDatasets{iDataset} 
        case 'occci_1km_1day_chl'
            chlaSatelliteCompleteOneDay1km = chlaSatelliteComplete;
            timeVectorCompleteOneDay1km = timeVectorComplete;
            seasonOneDay1km = M.season';
            fracMissingPixelsOneDay1km = fracMissingPixels;
            isValidImageOneDay1km = logical(isValidImage);
            imageDistribByMonthAndValidityOneDay1km = imageDistribByMonthAndValidity;
            fracMissingPixelsByMonthOneDay1km = fracMissingPixelsByMonth;
        case 'occci_4km_5day_chl'
            chlaSatelliteCompleteFiveDay4km = chlaSatelliteComplete;
            timeVectorCompleteFiveDay4km = timeVectorComplete;
            seasonFiveDay4km = M.season';
            fracMissingPixelsFiveDay4km = fracMissingPixels;
            isValidImageFiveDay4km = logical(isValidImage);
            imageDistribByMonthAndValidityFiveDay4km = imageDistribByMonthAndValidity;
            fracMissingPixelsByMonthFiveDay4km = fracMissingPixelsByMonth;
        case 'occci_4km_8day_chl'
            chlaSatelliteCompleteEightDay4km = chlaSatelliteComplete;
            timeVectorCompleteEightDay4km = timeVectorComplete;
            seasonEightDay4km = M.season';
            fracMissingPixelsEightDay4km = fracMissingPixels;
            isValidImageEightDay4km = logical(isValidImage);
            imageDistribByMonthAndValidityEightDay4km = imageDistribByMonthAndValidity;
            fracMissingPixelsByMonthEightDay4km = fracMissingPixelsByMonth;
    end
    
end

% Some plots

if (isPlotSatelliteProductInterComparison)

    plotMonthlyDistributionOfValidVsInvalidImages(occciDatasets,...
        fracMissingPixelToConsiderInvalidImage,...
        imageDistribByMonthAndValidityOneDay1km,...
        imageDistribByMonthAndValidityFiveDay4km,...
        imageDistribByMonthAndValidityEightDay4km,...
        fullPathPlotsDir)

    plotMonthlyDistributionOfValidVsInvalidPixels(occciDatasets,...
        fracMissingPixelsByMonthOneDay1km,...
        fracMissingPixelsByMonthFiveDay4km,...
        fracMissingPixelsByMonthEightDay4km,...
        fullPathPlotsDir)

    plotOverallDistributionOfInvalidPixels(occciDatasets,...
        fracMissingPixelsOneDay1km,isValidImageOneDay1km,...
        fracMissingPixelsFiveDay4km,isValidImageFiveDay4km,...
        fracMissingPixelsEightDay4km,isValidImageEightDay4km,...
        fullPathPlotsDir)

    plotSeasonalDistributionOfInvalidPixels(occciDatasets,...
        fracMissingPixelsOneDay1km,isValidImageOneDay1km,seasonOneDay1km,...
        fracMissingPixelsFiveDay4km,isValidImageFiveDay4km,seasonFiveDay4km,...
        fracMissingPixelsEightDay4km,isValidImageEightDay4km,seasonEightDay4km,...
        fullPathPlotsDir)

    plotPercentageMissingDataOverTime(occciDatasets,...
        fracMissingPixelsOneDay1km,timeVectorCompleteOneDay1km,...
        fracMissingPixelsFiveDay4km,timeVectorCompleteFiveDay4km,...
        fracMissingPixelsEightDay4km,timeVectorCompleteEightDay4km,...
        fullPathPlotsDir)

    plotPercentageMissingDataForSpecificYears(occciDatasets,OC_occci,...
        chlaSatelliteCompleteOneDay1km,fracMissingPixelsOneDay1km,timeVectorCompleteOneDay1km,...
        chlaSatelliteCompleteFiveDay4km,fracMissingPixelsFiveDay4km,timeVectorCompleteFiveDay4km,...
        chlaSatelliteCompleteEightDay4km,fracMissingPixelsEightDay4km,timeVectorCompleteEightDay4km,...
        areaStudy,enduranceSite,fullPathPlotsDir)

end

% Our chosen product is the OC-CCI 5-day, 4 km resolution product. The
% following analyses are based on that product only.
chlaProductChosen = 'occci_4km_5day_chl';
posIdChlaProdChosen = find(strcmp({OC_occci.dataset_name}, chlaProductChosen));

%% Remove outliers

% 1) Explore anomalously high values.
% 2) For those pixels that have a high value (> 10 mg chla m-3), let's 
% explore how the surrounding pixels look like.
% 3) Filter out those anomaly high values.

[nSatelliteLat,nSatelliteLon,nSatelliteTimeSteps] = size(chlaSatelliteCompleteFiveDay4km);
nPixelsInImage = nSatelliteLat*nSatelliteLon;
latSatellite = OC_occci(posIdChlaProdChosen).lat;
lonSatellite = OC_occci(posIdChlaProdChosen).lon;

if (isPlotScenePercentageMissingData)
    plotScenePercentageMissingDataByYear(...
        chlaSatelliteCompleteFiveDay4km,timeVectorCompleteFiveDay4km,...
        latSatellite,lonSatellite,areaStudy,enduranceSite,fullPathPlotsDir)
    plotScenePercentageMissingDataByMonth(...
        chlaSatelliteCompleteFiveDay4km,timeVectorCompleteFiveDay4km,...
        latSatellite,lonSatellite,areaStudy,enduranceSite,fullPathPlotsDir)
end

if (isPlotSatelliteHistogram)
%     chlaSatelliteCompleteFiveDay4km(isnan(chlaSatelliteCompleteFiveDay4km)) = 0;
    plotHistogramOfSatelliteValues(chlaSatelliteCompleteFiveDay4km,fullPathPlotsDir)
end

% Apply a mask to those very high values (defined by a threshold)
chlaHighThreshold = 10; % mg chla m-3
indices = find(chlaSatelliteCompleteFiveDay4km > chlaHighThreshold);
[subscripts(:, 1), subscripts(:, 2), subscripts(:, 3)] = ind2sub(size(chlaSatelliteCompleteFiveDay4km), indices);
uniqueSceneIndexes = unique(subscripts(:, 3)); 
chlaSatelliteCompleteFiveDay4kmMasked = chlaSatelliteCompleteFiveDay4km;
neighborhoodSize = 3; % 3x3 neighborhood
for iScene = 1:length(uniqueSceneIndexes)
    thisScene = squeeze(chlaSatelliteCompleteFiveDay4km(:,:,uniqueSceneIndexes(iScene)));
    thisScene(thisScene == 0) = NaN;
    % Perform a local percentile-based thresholding using a custom filter
    anomalyMask = customAnomalyFilter(thisScene, neighborhoodSize, chlaHighThreshold);
    % Replace values in data corresponding to mask equal to 1 with NaN
    thisScene(anomalyMask == 1) = NaN;
    % Update the satellite product
    chlaSatelliteCompleteFiveDay4kmMasked(:,:,uniqueSceneIndexes(iScene)) = thisScene;
end

% Plot scenes with high chlorophyll values before and after applying the mask
if (isPlotHighValuesMask)
    plotComparisonScenesBeforeAndAfterHighValueMask(uniqueSceneIndexes,...
        chlaSatelliteCompleteFiveDay4km,chlaSatelliteCompleteFiveDay4kmMasked,...
        timeVectorCompleteFiveDay4km,latSatellite,lonSatellite,...
        areaStudy,enduranceSite,fullPathPlotsDir)
end

%% Fill in data gaps

% 1) Apply linear interpolation, as explained in Racault et al. (2014), to
% fill in missing data.
% 2) Fill in persistent missing data (winter months) by inserting the data 
% from the CMEMS regional biogeochemical model reanalysis product of chl.

% Apply linear interpolation, as explained in Racault et al. (2014), to fill in persistent data gaps
chlaSatelliteCompleteFiveDay4kmInt = Racault2014interpolation(chlaSatelliteCompleteFiveDay4kmMasked);

% Compute the percentage of missing data before and after interpolation
nanCount = sum(isnan(reshape(chlaSatelliteCompleteFiveDay4kmMasked,[],1)));
nonNanCount = sum(~isnan(reshape(chlaSatelliteCompleteFiveDay4kmMasked,[],1)));
fracNanBefore = 100 * (nanCount / (nanCount + nonNanCount));
nanCount = sum(isnan(reshape(chlaSatelliteCompleteFiveDay4kmInt,[],1)));
nonNanCount = sum(~isnan(reshape(chlaSatelliteCompleteFiveDay4kmInt,[],1)));
fracNanAfter = 100 * (nanCount / (nanCount + nonNanCount));
disp(['Fraction of NaN pixels before interpolation: ',num2str(fracNanBefore)]);
disp(['Fraction of NaN pixels after interpolation: ',num2str(fracNanAfter)]);

% Randomly choose a scene and plot it before and after interpolation
if (isPlotInterpolation)
    % nDates = numel(timeVectorCompleteFiveDay4km);
    % randomTimeIndex = randi(nDates);
    randomTimeIndex = 1469;
    plotComparisonScenesBeforeAndAfterInterpolation(randomTimeIndex,...
        chlaSatelliteCompleteFiveDay4kmMasked,...
        chlaSatelliteCompleteFiveDay4kmInt,...
        latSatellite,lonSatellite,timeVectorCompleteFiveDay4km,...
        areaStudy,enduranceSite,fullPathPlotsDir)
end

% Explore where the persistent missing data are after interpolation
fracMissingPixelToConsiderInvalidImage = 98;
[fracMissingPixels,isValidImage,fracMissingPixelsByMonth,...
    imageDistribByMonthAndValidity,M] = calculateMissingPixelsStatistics(...
    chlaSatelliteCompleteFiveDay4kmInt,timeVectorCompleteFiveDay4km,...
    nSatelliteTimeSteps,nPixelsInImage,fracMissingPixelToConsiderInvalidImage);

if (isPlotDistributionPersistentMissingData)
    plotMonthlyDistributionOfValidVsInvalidImagesOneDataset(fracMissingPixelToConsiderInvalidImage,...
        imageDistribByMonthAndValidity,fullPathPlotsDir)
end

% Fill in the persistent data gaps (winter months) inserting the chla from 
% the CMEMS regional biogeochemical model reanalysis product

% Extract CMEMS data
cmemsChlaProduct = 'mod_bgc_reg_chl';
posIdCmemsChlaProd = find(strcmp({OC_cmems.dataset_name}, cmemsChlaProduct));

latCmemsModel = OC_cmems(posIdCmemsChlaProd).lat;
lonCmemsModel = OC_cmems(posIdCmemsChlaProd).lon;
timeCmemsModelAll = OC_cmems(posIdCmemsChlaProd).time;
iStartDate = find(timeCmemsModelAll == startDateSatellite);
timeCmemsModel = OC_cmems(posIdCmemsChlaProd).time(iStartDate:end);
chlaCmemsModel = OC_cmems(posIdCmemsChlaProd).dataset(:,:,iStartDate:end,:); % 'mod_bgc_reg_chl'
depthsCmemsModel = OC_cmems(posIdCmemsChlaProd).depth; % there is only data in the first 6 depth levels (0-30 m)
         
nCmemsLat = length(latCmemsModel);
nCmemsLon = length(lonCmemsModel);
nCmemsTime = length(timeCmemsModel);
nCmemsTimeAll = length(timeCmemsModelAll);

% Plot MLD and euphotic layer depth from CMEMS model to get an idea of the
% water column depth that needs to be averaged
if (isPlotEuphoticLayer)
    cmemsMldProduct = 'mod_phy_reg_mld'; % m
    posIdCmemsMldProd = find(strcmp({OC_cmems.dataset_name}, cmemsMldProduct));
    cmemsKdProduct = 'mod_bgc_reg_kd'; % d-1
    posIdCmemsKdProd = find(strcmp({OC_cmems.dataset_name}, cmemsKdProduct));

    plotCMEMSmixedLayerDepthAndEuphoticZoneDepth(...
        OC_cmems(posIdCmemsMldProd).dataset,...
        OC_cmems(posIdCmemsKdProd).dataset,...
        nCmemsLat,nCmemsLon,timeCmemsModelAll,nCmemsTimeAll,fullPathPlotsDir)
end

% Based on the plot above, we'll take data from the first 30 m depth
targetZeu = 30; % m 
[~,iZeu] = min(abs(depthsCmemsModel - targetZeu));
chlaCmemsModelEuphotic = zeros(nCmemsLat,nCmemsLon,nCmemsTime);
for iTime = 1:nCmemsTime
    for iLon = 1:nCmemsLon
        for iLat = 1:nCmemsLat
            chlaCmemsModelEuphotic(iLat,iLon,iTime) = mean(chlaCmemsModel(iLat,iLon,iTime,1:iZeu),4,'omitnan');
        end
    end
end

% Regrid CMEMS model output from daily, 7 km resolution to 5-day, 4 km
% resolution

% For that, first create the new time grid, where each year must start on 
% Jan 1st.
newTimeCmemsModel = timeCmemsModel(1);
while newTimeCmemsModel(end) + 5 <= timeCmemsModel(end)
    newDate = newTimeCmemsModel(end) + 5;
    % Check if the new date is in the same year
    if year(newDate) ~= year(newTimeCmemsModel(end))
        % Adjust the date to the 1st of January of the next year
        newDate = datetime(['01/01/' num2str(year(newDate))], 'InputFormat', 'dd/MM/yyyy');
    end
    newTimeCmemsModel = [newTimeCmemsModel; newDate];
end

% New grid with query points for interpolation/regridding
[qX,qY,qT] = ndgrid(latSatellite,lonSatellite,datenum(newTimeCmemsModel));

% Old grid
[X,Y,T] = ndgrid(latCmemsModel,lonCmemsModel,datenum(timeCmemsModel));

% Create interpolant
F = griddedInterpolant(X,Y,T,chlaCmemsModelEuphotic);

% Evaluate interpolant on new grid
chlaCmemsModelEuphoticFiveDay4km = F(qX,qY,qT);

% Check that the regridding went as expected
if (isPlotCMEMSregridded)
    plotComparisonCMEMSchlaBeforeAndAfterRegridding(...
        chlaCmemsModelEuphotic,latCmemsModel,lonCmemsModel,timeCmemsModel,...
        chlaCmemsModelEuphoticFiveDay4km,latSatellite,lonSatellite,newTimeCmemsModel,...
        areaStudy,enduranceSite,fullPathPlotsDir)
end

% Co-locate missing data in the satellite product and fill it with the CMEMS data
chlaSatelliteCompleteFiveDay4kmIntMissingFilled = NaN(size(chlaSatelliteCompleteFiveDay4kmInt));
isPermanentlyMissing = zeros(size(chlaSatelliteCompleteFiveDay4kmInt));
nSatelliteTimeSteps = numel(timeVectorCompleteFiveDay4km);
for iLon = 1:nSatelliteLon
    for iLat = 1:nSatelliteLat
        for iTimeStep = 1:nSatelliteTimeSteps
            if (iTimeStep > length(newTimeCmemsModel))
                break
            else
                if (isnan(chlaSatelliteCompleteFiveDay4kmInt(iLat,iLon,iTimeStep)))
                    chlaSatelliteCompleteFiveDay4kmIntMissingFilled(iLat,iLon,iTimeStep) =...
                        chlaCmemsModelEuphoticFiveDay4km(iLat,iLon,iTimeStep);
                    isPermanentlyMissing(iLat,iLon,iTimeStep) = 1;
                elseif (~isnan(chlaSatelliteCompleteFiveDay4kmInt(iLat,iLon,iTimeStep)))
                    chlaSatelliteCompleteFiveDay4kmIntMissingFilled(iLat,iLon,iTimeStep) =...
                        chlaSatelliteCompleteFiveDay4kmInt(iLat,iLon,iTimeStep);
                end
            end
        end
    end
end
isPermanentlyMissing = logical(isPermanentlyMissing);

% For the last year of the satellite product, for which we do not have
% CMEMS data to fill in persistent missing data, use the CMEMS data of the
% year before.
startDate = '2023-01-01';
iStartDate = find(timeVectorCompleteFiveDay4km == startDate);
for iLon = 1:nSatelliteLon
    for iLat = 1:nSatelliteLat
        for iTimeStep = iStartDate:nSatelliteTimeSteps
            if (isnan(chlaSatelliteCompleteFiveDay4kmInt(iLat,iLon,iTimeStep)))
                thisYearDate = timeVectorCompleteFiveDay4km(iTimeStep);
                yearAgoDate = thisYearDate - days(365);
                iTimeStepYearBefore = find(timeVectorCompleteFiveDay4km == yearAgoDate);
                chlaSatelliteCompleteFiveDay4kmIntMissingFilled(iLat,iLon,iTimeStep) =...
                    chlaCmemsModelEuphoticFiveDay4km(iLat,iLon,iTimeStepYearBefore);  
                isPermanentlyMissing(iLat,iLon,iTimeStep) = 1;
            elseif (~isnan(chlaSatelliteCompleteFiveDay4kmInt(iLat,iLon,iTimeStep)))
                 chlaSatelliteCompleteFiveDay4kmIntMissingFilled(iLat,iLon,iTimeStep) =...
                     chlaSatelliteCompleteFiveDay4kmInt(iLat,iLon,iTimeStep);
            end
        end
    end
end

continuousSatelliteChlaFiveDay4km = chlaSatelliteCompleteFiveDay4kmIntMissingFilled;

% Check that the filling process went well
if (isPlotFilledProduct)
    plotSatelliteTimeSeriesWithGapsFilledWithCMEMSdata(...
        chlaSatelliteCompleteFiveDay4kmIntMissingFilled,...
        isPermanentlyMissing,...
        latSatellite,lonSatellite,nSatelliteLat,nSatelliteLon,...
        timeVectorCompleteFiveDay4km,fullPathPlotsDir)
end

% Compare satellite data and CMEMS model data with in situ (SmartBuoy) data
if (isPlotComparisonWithInsituData)
    plotComparisonSmartBuoyDataWithSatelliteDataAndModelReanalysis(...
        chlaSatelliteCompleteFiveDay4kmMasked,timeVectorCompleteFiveDay4km,...
        chlaCmemsModelEuphoticFiveDay4km,newTimeCmemsModel,...
        fullPathMainDir,fullPathPlotsDir)
end


%% Save

save(strcat(fullPathMainDir,'continuousSatelliteChlaFiveDay4km.mat'),...
    'continuousSatelliteChlaFiveDay4km','timeVectorCompleteFiveDay4km',...
    'latSatellite','lonSatellite')

%% Local plotting functions to this script

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMonthlyDistributionOfValidVsInvalidImages(occciDatasets,...
    fracToTest,imageDistribByMonthAndValidityOneDay1km,...
    imageDistribByMonthAndValidityFiveDay4km,...
    imageDistribByMonthAndValidityEightDay4km,fullPathPlotsDir)

% Plot the frequency distribution by month of invalid vs valid images

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.67],'Color','w')
iSubplot = 0;

for iDataset = 1:length(occciDatasets)

    switch occciDatasets{iDataset} 
        case 'occci_1km_1day_chl'
            barData = imageDistribByMonthAndValidityOneDay1km;
            titleStr = {'Daily, 1 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        case 'occci_4km_5day_chl'
            barData = imageDistribByMonthAndValidityFiveDay4km;
            titleStr = {'5-day, 4 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        case 'occci_4km_8day_chl'
            barData = imageDistribByMonthAndValidityEightDay4km;
            titleStr = {'8-day, 4 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        otherwise
            continue
    end
    
    haxis(iSubplot) = subaxis(3,1,iSubplot,'Spacing',0.01,'Padding',0.025,'Margin', 0.10);
    
    ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.03;
    if (iSubplot == 1)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.025;
    elseif (iSubplot == 2)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2);
    elseif (iSubplot == 3)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.025;
    end
    set(haxis(iSubplot),'Position',ax(iSubplot).pos)
    
    hbar = bar(haxis(iSubplot),barData,'BarWidth',1.0);
    % Chance the colors of each bar segment
    colors = parula(size(hbar,2)); 
    colors = repelem(colors,size(hbar,1),1); 
    colors = mat2cell(colors,ones(size(colors,1),1),3);
    set(hbar,{'FaceColor'},colors) 
    
    if (iSubplot > 1)
        ylim([0 200])
    end
    
    title(titleStr,'FontSize',12);
    if (iSubplot == 3)
        xlabel('Month');
    end
    xticks(1:12);
    xticklabels(months);
    ylabel('Frequency');
    if (iSubplot == 1)
        lg = legend('Valid','Invalid');
        lg.Position(1) = 0.85; lg.Position(2) = 0.85; 
        lg.Orientation = 'vertical';
        lg.FontSize = 11; 
        lg.ItemTokenSize = [10,20];
        set(lg,'Box','on','Color','w')
    end
    
    % Create a wider title spanning both subplots
    annotation('textbox',[0 0.5 1 0.5],... 
        'String',sprintf('Distribution of invalid images (> %2.0f%% of invalid pixels)',fracToTest),...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');

end

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,strcat('plot_comparison_valid_images_',num2str(fracToTest))),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    strcat('plot_comparison_valid_images_',num2str(fracToTest),'.png')),'Resolution',600)

end % plotMonthlyDistributionOfValidVsInvalidImages

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMonthlyDistributionOfValidVsInvalidImagesOneDataset(fracToTest,...
    imageDistribByMonthAndValidity,fullPathPlotsDir)

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.42 0.25],'Color','w')
axh = axes('Position', [0.10, 0.15, 0.75, 0.70]);

hbar = bar(axh,imageDistribByMonthAndValidity,'BarWidth',1.0);

% Chance the colors of each bar segment
colors = parula(size(hbar,2)); 
colors = repelem(colors,size(hbar,1),1); 
colors = mat2cell(colors,ones(size(colors,1),1),3);
set(hbar,{'FaceColor'},colors) 

ylim([0 200])
title(sprintf('Distribution of invalid images (> %2.0f%% of invalid pixels) after interpolation',fracToTest),'FontSize',12);
xlabel('Month');
xticks(1:12);
xticklabels(months);
ylabel('Frequency');

lg = legend('Valid','Invalid');
lg.Position(1) = 0.85; lg.Position(2) = 0.73; 
lg.Orientation = 'vertical';
lg.FontSize = 11; 
lg.ItemTokenSize = [10,20];
set(lg,'Box','on','Color','w')

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,strcat('plot_distribution_valid_images_after_interpolation_',num2str(fracToTest))),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    strcat('plot_distribution_valid_images_after_interpolation_',num2str(fracToTest),'.png')),'Resolution',600)

end % plotMonthlyDistributionOfValidVsInvalidImagesOneDataset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMonthlyDistributionOfValidVsInvalidPixels(occciDatasets,...
    fracMissingPixelsByMonthOneDay1km,fracMissingPixelsByMonthFiveDay4km,...
    fracMissingPixelsByMonthEightDay4km,fullPathPlotsDir)

% For valid images, plot the frequency distribution of invalid pixels by month

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.35 0.67],'Color','w')
iSubplot = 0;

for iDataset = 1:length(occciDatasets)
    
    switch occciDatasets{iDataset} 
        case 'occci_1km_1day_chl'
            barData = fracMissingPixelsByMonthOneDay1km;
            titleStr = {'Distribution of the mean percentage of missing pixels in valid images',... 
                'in the daily, 1 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        case 'occci_4km_5day_chl'
            barData = fracMissingPixelsByMonthFiveDay4km;
            titleStr = {'Distribution of the mean percentage of missing pixels in valid images',... 
                'in the 5-day, 4 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        case 'occci_4km_8day_chl'
            barData = fracMissingPixelsByMonthEightDay4km;
            titleStr = {'Distribution of the mean percentage of missing pixels in valid images',... 
                'in the 8-day, 4 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        otherwise
            continue
    end
    
    haxis(iSubplot) = subaxis(3,1,iSubplot,'Spacing',0.01,'Padding',0.025,'Margin', 0.10);
    
    ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    if (iSubplot == 1)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.025;
    elseif (iSubplot == 2)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2);
    elseif (iSubplot == 3)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.025;
    end
    set(haxis(iSubplot),'Position',ax(iSubplot).pos)

    bar(haxis(iSubplot),barData,'BarWidth',1.0,'FaceColor',[0.7,0.7,0.7]);
    xticks(1:12);
    xticklabels(months);
    xlim([0.5, length(months) + 0.5]);
    yticks(0:0.2:1);
    yticklabels({'0','20','40','60','80'});
    ylim([0 0.8])
    histaxes = gca;
    histaxes.XGrid = 'off';  
    histaxes.YGrid = 'on';      
    title(titleStr,'FontSize',12);
    if (iSubplot == 3)
        xlabel('Month');
    end
    ylabel('% of missing pixels per image');
end

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,'plot_comparison_missing_pixels_bymonth'),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'plot_comparison_missing_pixels_bymonth.png'),'Resolution',600)

end % plotMonthlyDistributionOfValidVsInvalidPixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotOverallDistributionOfInvalidPixels(occciDatasets,...
    fracMissingPixelsOneDay1km,isValidImageOneDay1km,...
    fracMissingPixelsFiveDay4km,isValidImageFiveDay4km,...
    fracMissingPixelsEightDay4km,isValidImageEightDay4km,...
    fullPathPlotsDir)

% For valid images, plot the histogram of invalid pixels overall 

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.35 0.67],'Color','w')
iSubplot = 0;

for iDataset = 1:length(occciDatasets)
    
    switch occciDatasets{iDataset} 
        case 'occci_1km_1day_chl'
            histData = fracMissingPixelsOneDay1km(isValidImageOneDay1km == 1);
            titleStr = {'Distribution of the percentage of missing pixels in valid images',... 
                'in the daily, 1 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        case 'occci_4km_5day_chl'
            histData = fracMissingPixelsFiveDay4km(isValidImageFiveDay4km == 1);
            titleStr = {'Distribution of the percentage of missing pixels in valid images',... 
                'in the 5-day, 4 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        case 'occci_4km_8day_chl'
            histData = fracMissingPixelsEightDay4km(isValidImageEightDay4km == 1);
            titleStr = {'Distribution of the percentage of missing pixels in valid images',... 
                'in the 8-day, 4 km resolution OC-CCI dataset (1998-2023)'};
            iSubplot = iSubplot + 1;
        otherwise
            continue
    end
    
    haxis(iSubplot) = subaxis(3,1,iSubplot,'Spacing',0.010,'Padding',0.025,'Margin',0.10);
    
    ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    if (iSubplot == 1)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.025;
    elseif (iSubplot == 2)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2);
    elseif (iSubplot == 3)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.025;
    end
    set(haxis(iSubplot),'Position',ax(iSubplot).pos)

    histogram(haxis(iSubplot),histData,50,'BinLimits',[0,1],'FaceColor',[0.5, 0.5, 0.5])
    set(gca,'YScale','log');
    ylim([1 1e3])
    yticks([1, 1e1, 1e2, 1e3]);
    yticklabels({'1','10','10^{2}','10^{3}'});
    grid on;
    mainAxes = gca;
    mainAxes.XGrid = 'on'; 
    mainAxes.YGrid = 'on'; 
    mainAxes.MinorGridAlpha = 0;
    xticks(0:0.1:1);
    xticklabels({'0','10','20','30','40','50','60','70','80','90','100'});
    if (iSubplot == 2)
        xlabel('% of missing pixels per image');
    end
    ylabel('Frequency');
    title(titleStr);
end

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,'plot_comparison_missing_pixels_overall'),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'plot_comparison_missing_pixels_overall.png'),'Resolution',600)

end % plotOverallDistributionOfInvalidPixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSeasonalDistributionOfInvalidPixels(occciDatasets,...
    fracMissingPixelsOneDay1km,isValidImageOneDay1km,seasonOneDay1km,...
    fracMissingPixelsFiveDay4km,isValidImageFiveDay4km,seasonFiveDay4km,...
    fracMissingPixelsEightDay4km,isValidImageEightDay4km,seasonEightDay4km,...
    fullPathPlotsDir)

% For valid images, plot the histogram of invalid pixels by season

seasons = {'Winter','Spring','Summer','Autumn'};
colourScheme = brewermap(length(seasons),'*Spectral');

for iDataset = 1:length(occciDatasets)
    
    switch occciDatasets{iDataset} 
        case 'occci_1km_1day_chl'
            overallData = fracMissingPixelsOneDay1km(isValidImageOneDay1km == 1);
        case 'occci_4km_5day_chl'
            overallData = fracMissingPixelsFiveDay4km(isValidImageFiveDay4km == 1);
        case 'occci_4km_8day_chl'
            overallData = fracMissingPixelsEightDay4km(isValidImageEightDay4km == 1);
        otherwise
            continue
    end
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.40],'Color','w') 
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
                if (strcmp(occciDatasets{iDataset},'occci_1km_1day_chl'))
                    seasonData = fracMissingPixelsOneDay1km(isValidImageOneDay1km == 1 &...
                        seasonOneDay1km == "Spring");
                elseif (strcmp(occciDatasets{iDataset},'occci_4km_5day_chl'))
                    seasonData = fracMissingPixelsFiveDay4km(isValidImageFiveDay4km == 1 &...
                        seasonFiveDay4km == "Spring");
                elseif (strcmp(occciDatasets{iDataset},'occci_4km_8day_chl'))
                    seasonData = fracMissingPixelsEightDay4km(isValidImageEightDay4km == 1 &...
                        seasonEightDay4km == "Spring");
                end
            case 'Summer'
                if (strcmp(occciDatasets{iDataset},'occci_1km_1day_chl'))
                    seasonData = fracMissingPixelsOneDay1km(isValidImageOneDay1km == 1 &...
                        seasonOneDay1km == "Summer");
                elseif (strcmp(occciDatasets{iDataset},'occci_4km_5day_chl'))
                    seasonData = fracMissingPixelsFiveDay4km(isValidImageFiveDay4km == 1 &...
                        seasonFiveDay4km == "Summer");
                elseif (strcmp(occciDatasets{iDataset},'occci_4km_8day_chl'))
                    seasonData = fracMissingPixelsEightDay4km(isValidImageEightDay4km == 1 &...
                        seasonEightDay4km == "Summer");
                end
            case 'Autumn'
                if (strcmp(occciDatasets{iDataset},'occci_1km_1day_chl'))
                    seasonData = fracMissingPixelsOneDay1km(isValidImageOneDay1km == 1 &...
                        seasonOneDay1km == "Autumn");
                elseif (strcmp(occciDatasets{iDataset},'occci_4km_5day_chl'))
                    seasonData = fracMissingPixelsFiveDay4km(isValidImageFiveDay4km == 1 &...
                        seasonFiveDay4km == "Autumn");
                elseif (strcmp(occciDatasets{iDataset},'occci_4km_8day_chl'))
                    seasonData = fracMissingPixelsEightDay4km(isValidImageEightDay4km == 1 &...
                        seasonEightDay4km == "Autumn");
                end
            case 'Winter'
                if (strcmp(occciDatasets{iDataset},'occci_1km_1day_chl'))
                    seasonData = fracMissingPixelsOneDay1km(isValidImageOneDay1km == 1 &...
                        seasonOneDay1km == "Winter");
                elseif (strcmp(occciDatasets{iDataset},'occci_4km_5day_chl'))
                    seasonData = fracMissingPixelsFiveDay4km(isValidImageFiveDay4km == 1 &...
                        seasonFiveDay4km == "Winter");
                elseif (strcmp(occciDatasets{iDataset},'occci_4km_8day_chl'))
                    seasonData = fracMissingPixelsEightDay4km(isValidImageEightDay4km == 1 &...
                        seasonEightDay4km == "Winter");
                end
        end

        % Plot all data
        histogram(haxis(iSeason),overallData,50,'BinLimits',[0,1],'FaceColor',[0.7,0.7,0.7],'EdgeColor','none')
        hold on
        % Plot seasonal data
        histogram(haxis(iSeason),seasonData,50,'BinLimits',[0,1],'FaceColor',colourScheme(iSeason,:))
        hold on
        set(gca,'YScale','log');
        ylim([1 1e3])
        grid on;
        mainSeasonalAxes = gca;
        mainSeasonalAxes.XGrid = 'on'; 
        mainSeasonalAxes.YGrid = 'on'; 
        mainSeasonalAxes.MinorGridAlpha = 0; 
%         yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]);
%         yticklabels({'1','10','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}'});

        if (iSeason == 1 || iSeason == 3)
            ylabel('Frequency');
        end
        if (iSeason == 3 || iSeason == 4)
            xlabel('% missing pixels per image');
        end
        title(seasons{iSeason});

    end % iSeason
    hold off

%     set(gcf,'PaperPositionMode','auto')
%     print(gcf,fullfile(fullPathPlotsDir,strcat('comparison_missing_pixels_byseason_',occciDatasets{iDataset})),'-dpdf','-r0')
    exportgraphics(gcf,fullfile(fullPathPlotsDir,...
        strcat('comparison_missing_pixels_byseason_',occciDatasets{iDataset},'.png')),'Resolution',600)

end

end % plotSeasonalDistributionOfInvalidPixels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotPercentageMissingDataOverTime(occciDatasets,...
    fracMissingPixelsOneDay1km,timeVectorOneDay1km,...
    fracMissingPixelsFiveDay4km,timeVectorFiveDay4km,...
    fracMissingPixelsEightDay4km,timeVectorEightDay4km,...
    fullPathPlotsDir)

% Plot the percentage of missing data over time

for iDataset = 1:length(occciDatasets)
    
    switch occciDatasets{iDataset} 
        case 'occci_1km_1day_chl'
            timeVector = timeVectorOneDay1km;
            fracMissingPixels = fracMissingPixelsOneDay1km;
        case 'occci_4km_5day_chl'
            timeVector = timeVectorFiveDay4km;
            fracMissingPixels = fracMissingPixelsFiveDay4km;
        case 'occci_4km_8day_chl'
            timeVector = timeVectorEightDay4km;
            fracMissingPixels = fracMissingPixelsEightDay4km;
        otherwise
            continue
    end
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.48 0.45],'Color','w')
    haxis = zeros(2,1);
    
    for iSubplot = 1:2
        
        haxis(iSubplot) = subaxis(2,1,iSubplot,'Spacing',0.01,'Padding',0.01,'Margin', 0.10);

        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        if (iSubplot == 1)
            ax(1).pos(2) = ax(1).pos(2) + 0.07;
        elseif (iSubplot == 2)
            ax(2).pos(2) = ax(2).pos(2) - 0.03;
        end
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1) + 0.03;
        set(haxis(iSubplot),'Position',ax(iSubplot).pos)
    
        if (iSubplot == 1)
            scatter(haxis(iSubplot),timeVector,100.*fracMissingPixels,6,'b','filled'); 
            xlim([datetime('1998-01-01') datetime('2022-12-31')])
            xlabel('Time');
            ylabel('% missing pixels per image');
            title('Time distribution of missing data')
        elseif (iSubplot == 2)
            histogram(haxis(iSubplot),100.*fracMissingPixels,50,'BinLimits',[0,100],'FaceColor',[0.7,0.7,0.7],'EdgeColor','black')
            set(gca,'YScale','log');
            ylim([1 5500])
            histaxes = gca;
            histaxes.XGrid = 'on';  
            histaxes.YGrid = 'on';  
            histaxes.MinorGridAlpha = 0;
            xlabel('% missing pixels per image');
            ylabel('Frequency');
            title('Histogram of missing data')
        end
        grid on
        box on
    end
    
%     set(gcf,'PaperPositionMode','auto')
%     print(gcf,fullfile(fullPathPlotsDir,strcat('plot_percentage_missing_data_',occciDatasets{iDataset})),'-dpdf','-r0')  
    exportgraphics(gcf,fullfile(fullPathPlotsDir,...
        strcat('plot_percentage_missing_data_',occciDatasets{iDataset},'.png')),'Resolution',600)
    
end

end % plotPercentageMissingDataOverTime
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotPercentageMissingDataForSpecificYears(occciDatasets,OC_occci,...
    chlaSatelliteOneDay1km,fracMissingPixelsOneDay1km,timeVectorOneDay1km,...
    chlaSatelliteFiveDay4km,fracMissingPixelsFiveDay4km,timeVectorFiveDay4km,...
    chlaSatelliteEightDay4km,fracMissingPixelsEightDay4km,timeVectorEightDay4km,...
    areaStudy,enduranceSite,fullPathPlotsDir)

% Plot percentage of missing data for 1998 and 2020

for iDataset = 1:length(occciDatasets)
    
    switch occciDatasets{iDataset} 
        case 'occci_1km_1day_chl'
            
            iFirstDate1998 = find(timeVectorOneDay1km == '1998-01-01');
            iEndDate1998 = find(timeVectorOneDay1km == '1998-12-31');
            iFirstDate2020 = find(timeVectorOneDay1km == '2020-01-01');
            iEndDate2020 = find(timeVectorOneDay1km == '2020-12-31');
            nImages1998 = iEndDate1998 - iFirstDate1998 + 1;
            nImages2020 = iEndDate2020 - iFirstDate2020 + 1;
            [nRows,nCols,~] = size(chlaSatelliteOneDay1km);
            meanFracMisingData1998 = zeros(nRows,nCols);
            meanFracMisingData2020 = zeros(nRows,nCols);
            for iRow = 1:nRows
                for iCol = 1:nCols
                    meanFracMisingData1998(iRow,iCol) = 100 * (sum(isnan(squeeze(chlaSatelliteOneDay1km(iRow,iCol,iFirstDate1998:iEndDate1998))))/nImages1998); 
                    meanFracMisingData2020(iRow,iCol) = 100 * (sum(isnan(squeeze(chlaSatelliteOneDay1km(iRow,iCol,iFirstDate2020:iEndDate2020))))/nImages2020); 
                end
            end
            
            timeVector = timeVectorOneDay1km;
            fracMissingPixels = fracMissingPixelsOneDay1km;
            posIdProd = find(strcmp({OC_occci.dataset_name}, occciDatasets{iDataset}));
            lat = OC_occci(posIdProd).lat; 
            lon = OC_occci(posIdProd).lon; 

        case 'occci_4km_5day_chl'
            
            iFirstDate1998 = find(timeVectorFiveDay4km == '1998-01-01');
            [~,iEndDate1998] = min(abs(timeVectorFiveDay4km - '1998-12-31'));
            iFirstDate2020 = find(timeVectorFiveDay4km == '2020-01-01');
            [~,iEndDate2020] = min(abs(timeVectorFiveDay4km - '2020-12-31'));
            nImages1998 = iEndDate1998 - iFirstDate1998 + 1;
            nImages2020 = iEndDate2020 - iFirstDate2020 + 1;
            [nRows,nCols,~] = size(chlaSatelliteFiveDay4km);
            meanFracMisingData1998 = zeros(nRows,nCols);
            meanFracMisingData2020 = zeros(nRows,nCols);
            for iRow = 1:nRows
                for iCol = 1:nCols
                    meanFracMisingData1998(iRow,iCol) = 100 * (sum(isnan(squeeze(chlaSatelliteFiveDay4km(iRow,iCol,iFirstDate1998:iEndDate1998))))/nImages1998); 
                    meanFracMisingData2020(iRow,iCol) = 100 * (sum(isnan(squeeze(chlaSatelliteFiveDay4km(iRow,iCol,iFirstDate2020:iEndDate2020))))/nImages2020); 
                end
            end
            
            timeVector = timeVectorFiveDay4km;
            fracMissingPixels = fracMissingPixelsFiveDay4km;
            posIdProd = find(strcmp({OC_occci.dataset_name}, occciDatasets{iDataset}));
            lat = OC_occci(posIdProd).lat;
            lon = OC_occci(posIdProd).lon;
            
        case 'occci_4km_8day_chl'
            
            iFirstDate1998 = find(timeVectorEightDay4km == '1998-01-01');
            [~,iEndDate1998] = min(abs(timeVectorEightDay4km - '1998-12-31'));
            iFirstDate2020 = find(timeVectorEightDay4km == '2020-01-01');
            [~,iEndDate2020] = min(abs(timeVectorEightDay4km - '2020-12-31'));
            nImages1998 = iEndDate1998 - iFirstDate1998 + 1;
            nImages2020 = iEndDate2020 - iFirstDate2020 + 1;
            [nRows,nCols,~] = size(chlaSatelliteEightDay4km);
            meanFracMisingData1998 = zeros(nRows,nCols);
            meanFracMisingData2020 = zeros(nRows,nCols);
            for iRow = 1:nRows
                for iCol = 1:nCols
                    meanFracMisingData1998(iRow,iCol) = 100 * (sum(isnan(squeeze(chlaSatelliteEightDay4km(iRow,iCol,iFirstDate1998:iEndDate1998))))/nImages1998); 
                    meanFracMisingData2020(iRow,iCol) = 100 * (sum(isnan(squeeze(chlaSatelliteEightDay4km(iRow,iCol,iFirstDate2020:iEndDate2020))))/nImages2020); 
                end
            end
            
            timeVector = timeVectorEightDay4km;
            fracMissingPixels = fracMissingPixelsEightDay4km;
            posIdProd = find(strcmp({OC_occci.dataset_name}, occciDatasets{iDataset}));
            lat = OC_occci(posIdProd).lat;
            lon = OC_occci(posIdProd).lon;
            
        otherwise
            continue

    end
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.50],'Color','w') 
    haxis = zeros(2,2);

    for iSubplot = 1:4

        haxis(iSubplot) = subaxis(2,2,iSubplot,'Spacing',0.05,'Padding',0.05,'Margin', 0.05);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');

        if (iSubplot == 1)
            ax(iSubplot).pos(3) = 0.47;
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 
            plot(timeVector,100.*fracMissingPixels,...
                '-k','LineWidth',1,...
                'Marker','.','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b'); 
            xlim([datetime('1998-01-01') datetime('1998-12-31')])
        elseif (iSubplot == 2)
            ax(iSubplot).pos(1) = 0.68;
            ax(iSubplot).pos(3) = 0.28;
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 
            m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
            m_pcolor(lon,lat,meanFracMisingData1998) 
            title('% missing data per pixel 1998')
        elseif (iSubplot == 3)
            ax(iSubplot).pos(3) = 0.47;
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 
            plot(timeVector,100.*fracMissingPixels,...
                '-k','LineWidth',1,...
                'Marker','.','MarkerSize',10,'MarkerEdgeColor','b','MarkerFaceColor','b'); 
            xlim([datetime('2020-01-01') datetime('2020-12-31')])
        elseif (iSubplot == 4)
            ax(iSubplot).pos(1) = 0.68;
            ax(iSubplot).pos(3) = 0.28;
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 
            m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
            m_pcolor(lon,lat,meanFracMisingData2020)
            title('% missing data per pixel 2020')
        end

        grid on 
        box on

        if (iSubplot == 1 || iSubplot == 3)
            xlabel('Time');
            ylabel('% missing pixels per image');
        elseif (iSubplot == 2 || iSubplot == 4)
            m_grid('linewi',1,'tickdir','out','FontSize',8)
            m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
            shading flat
            colorbar
            colormap(brewermap(1000,'*Blues'))
            switch occciDatasets{iDataset} 
                case 'occci_1km_1day_chl'
                    caxis([60 100])
                case 'occci_4km_5day_chl'
                    caxis([0 50])
                case 'occci_4km_8day_chl'
                    caxis([0 50])
            end
            xlabel('Longitude');
            ylabel('Latitude'); 
        end

    %     set(gca, 'color', 'w'); % or whatever color you want for the background

    end

%     set(gcf,'PaperPositionMode','auto')
%     print(gcf,fullfile(fullPathPlotsDir,strcat('plot_percentage_missing_data_years_',occciDatasets{iDataset})),'-dpdf','-r0')  
    exportgraphics(gcf,fullfile(fullPathPlotsDir,...
        strcat('plot_percentage_missing_data_years_',occciDatasets{iDataset},'.png')),'Resolution',600)
end

end % plotPercentageMissingDataForSpecificYears

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotScenePercentageMissingDataByYear(chlaSatellite,timeVector,...
    latSatellite,lonSatellite,areaStudy,enduranceSite,fullPathPlotsDir)

yearVectorData = (min(year(timeVector)):1:max(year(timeVector)))';
yearFracMisingData = zeros(numel(latSatellite),numel(lonSatellite),numel(yearVectorData));
for iYear = 1:numel(yearVectorData)
    yearIdxs = find(year(timeVector) == yearVectorData(iYear));
    nImagesInYear = numel(yearIdxs);
    for iCol = 1:numel(lonSatellite)
        for iRow = 1:numel(latSatellite)
            yearFracMisingData(iRow,iCol,iYear) =...
                100*(sum(isnan(squeeze(chlaSatellite(iRow,iCol,yearIdxs))))/nImagesInYear);
        end
    end
end

figure()               
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.75],'Color','w') 
haxis = zeros(5,5);
        
for iYear = 1:(numel(yearVectorData)-1) % DO NOT SHOW LAST YEAR (2023)
            
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
    
    m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
    m_pcolor(lonSatellite,latSatellite,yearFracMisingData(:,:,iYear)) 
    shading flat
    colormap(brewermap(1000,'*Blues'))
    caxis([10 60])
    set(gca,'color','w'); % or whatever color you want for the background

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
    title(num2str(yearVectorData(iYear)),'FontSize',11)
    box on
    hold on
    
end

% Colour bar settings
cb = colorbar('Location','eastoutside');
cb.Position(1) = cb.Position(1) + 0.055;
cb.Position(2) = cb.Position(2) + 0.07;
cb.Position(3) = 0.020; % WIDTH
cb.Position(4) = 0.70; % LENGTH
cb.Label.String = '% missing data'; 
cb.FontSize = 10;

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,'plot_scene_percentage_missing_data_by_year'),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'plot_scene_percentage_missing_data_by_year.png'),'Resolution',600)

end % plotScenePercentageMissingDataByYear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotScenePercentageMissingDataByMonth(chlaSatellite,timeVector,...
    latSatellite,lonSatellite,areaStudy,enduranceSite,fullPathPlotsDir)

monthFracMisingData = zeros(numel(latSatellite),numel(lonSatellite),12);
for iMonth = 1:12
    monthIdxs = find(month(timeVector) == iMonth);
    nImagesInMonthForPeriod = numel(monthIdxs);
    for iCol = 1:numel(lonSatellite)
        for iRow = 1:numel(latSatellite)
            monthFracMisingData(iRow,iCol,iMonth) =...
                100*(sum(isnan(squeeze(chlaSatellite(iRow,iCol,monthIdxs))))/nImagesInMonthForPeriod);
        end
    end
end

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

figure()               
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.35 0.45],'Color','w') 
haxis = zeros(3,4);
        
for iMonth = 1:12
            
    haxis(iMonth) = subaxis(3,4,iMonth,'Spacing',0.015,'Padding',0.010,'Margin',0.04);
    ax(iMonth).pos = get(haxis(iMonth),'Position');

    % Move subplots
    if (iMonth == 1 || iMonth == 5 || iMonth == 9)
        ax(iMonth).pos(1) = ax(iMonth).pos(1);
    elseif (iMonth == 2 || iMonth == 6 || iMonth == 10)
        ax(iMonth).pos(1) = ax(iMonth).pos(1) - 0.025;
    elseif (iMonth == 3 || iMonth == 7 || iMonth == 11)
        ax(iMonth).pos(1) = ax(iMonth).pos(1) - 0.050;
    elseif (iMonth == 4 || iMonth == 8 || iMonth == 12)
        ax(iMonth).pos(1) = ax(iMonth).pos(1) - 0.075;
    end  
    set(haxis(iMonth),'Position',ax(iMonth).pos) 
    
    m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
    m_pcolor(lonSatellite,latSatellite,monthFracMisingData(:,:,iMonth)) 
    shading flat
    colormap(brewermap(1000,'*Blues'))
    caxis([10 100])
    set(gca,'color','w'); % or whatever color you want for the background

    % Show axis labels for border subplots
    if (iMonth == 1 || iMonth == 5)
        m_grid('linewi',1,'tickdir','out','FontSize',8,...
            'xticklabel',[]);
    elseif (iMonth == 10 || iMonth == 11 || iMonth == 12)
        m_grid('linewi',1,'tickdir','out','FontSize',8,...
            'yticklabel',[]);
    elseif (iMonth == 9)
        m_grid('linewi',1,'tickdir','out','FontSize',8)
    else
        m_grid('linewi',1,'tickdir','out','FontSize',8,...
            'xticklabel',[],'yticklabel',[]);
    end

    m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
    title(months{iMonth},'FontSize',11)
    box on
    hold on
    
end

% Colour bar settings
cb = colorbar('Location','eastoutside');
cb.Position(1) = cb.Position(1) + 0.075;
cb.Position(2) = cb.Position(2) + 0.055;
cb.Position(3) = 0.020; % WIDTH
cb.Position(4) = 0.70; % LENGTH
cb.Label.String = '% missing data in the period 1998-2022'; 
cb.FontSize = 10;

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,'plot_scene_percentage_missing_data_by_month'),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'plot_scene_percentage_missing_data_by_month.png'),'Resolution',600)

end % plotScenePercentageMissingDataByMonth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotHistogramOfSatelliteValues(chlaSatellite,fullPathPlotsDir)

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.29],'Color','w') 

% Main plot
histogram(chlaSatellite,100,'BinLimits',[1e-2,100],'FaceColor', [0.5, 0.5, 0.5])
set(gca,'YScale','log');
grid on;
mainAxes = gca;
mainAxes.XGrid = 'on';
mainAxes.YGrid = 'on';
mainAxes.MinorGridAlpha = 0;
yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]);
yticklabels({'1','10','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}','10^{7}','10^{8}'});
xlabel('Chlorophyll a (mg m^{-3})');
ylabel('Frequency');
title('Distribution of satellite observations from the OC-CCI, 5-day 4 km product');

% Inset plot
insetAxes = axes('Position', [0.4, 0.43, 0.47, 0.45]); % define the position of the inset plot
histogram(insetAxes,chlaSatellite,100,'BinLimits',[1e-2,2],'FaceColor', [0.5, 0.5, 0.5])
set(insetAxes,'YScale','log');
grid on;
insetAxes.XGrid = 'on';
insetAxes.YGrid = 'on';
insetAxes.MinorGridAlpha = 0;
insetAxes.XAxis.FontSize = 8;
insetAxes.YAxis.FontSize = 8;
yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6]);
yticklabels({'1','10','10^{2}','10^{3}','10^{4}','10^{5}','10^{6}'});
xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
ylabel('Frequency','FontSize',8)

% set(gcf,'PaperPositionMode', 'auto')
% print(gcf,fullfile(fullPathPlotsDir,'histogram_values_5day_product'),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'histogram_values_5day_product.png'),'Resolution',600)

end % plotHistogramOfSatelliteValues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotComparisonScenesBeforeAndAfterHighValueMask(uniqueSceneIndexes,...
    chlaSatelliteOriginal,chlaSatelliteMask,timeVectorSatellite,...
    latSatellite,lonSatellite,areaStudy,enduranceSite,fullPathPlotsDir)

% We have 24 cases where chla values > 10 mg chla m-3 

for iFigure = 1:2

    figure()               
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.75],'Color','w') 
    haxis = zeros(6,4);

    idx = 1;
    for i = 1:length(uniqueSceneIndexes)

        if (iFigure == 1)
            thisScene = squeeze(chlaSatelliteOriginal(:,:,uniqueSceneIndexes(idx)));
            thisScene(thisScene == 0) = NaN;
            figureSuffix = 'original';
        elseif (iFigure == 2)
            thisScene = squeeze(chlaSatelliteMask(:,:,uniqueSceneIndexes(idx)));
            thisScene(thisScene == 0) = NaN;
            figureSuffix = 'filtered';
        end

        haxis(i) = subaxis(6,4,i,'Spacing',0.015,'Padding',0.010,'Margin',0.04);
        ax(i).pos = get(haxis(i),'Position');

        % Move subplots
        if (i == 1 || i == 5 || i == 9 || i == 13 || i == 17 || i == 21)
            ax(i).pos(1) = ax(i).pos(1) + 0.025;
        elseif (i == 2 || i == 6 || i == 10 || i == 14 || i == 18 || i == 22)
            ax(i).pos(1) = ax(i).pos(1) + 0;
        elseif (i == 3 || i == 7 || i == 11 || i == 15 || i == 19 || i == 23)
            ax(i).pos(1) = ax(i).pos(1) - 0.025;
        elseif (i == 4 || i == 8 || i == 12 || i == 16 || i == 20 || i == 24)
            ax(i).pos(1) = ax(i).pos(1) - 0.050;  
        end  
        ax(i).pos(2) = ax(i).pos(2) + 0.035;
        set(haxis(i),'Position',ax(i).pos) 

        m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
        m_pcolor(lonSatellite,latSatellite,thisScene) 
%         imagesc(lon,lat,thisScene)
        shading flat
        colormap(brewermap(1000,'*RdYlBu'))
        caxis([0 15])
        set(gca,'color','w'); % or whatever color you want for the background
        m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
        m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
        title(datestr(timeVectorSatellite(uniqueSceneIndexes(idx)),'yyyy-mm-dd'),'FontSize',9)
        idx = idx + 1;
        box on
        hold on
    end

    % Colour bar settings
    cb = colorbar('Location','southoutside');
    cb.Position(1) = cb.Position(1) - 0.58;
    cb.Position(2) = cb.Position(2) - 0.07;
    cb.Position(3) = 0.60; % LENGTH
    cb.Position(4) = 0.015; % WIDTH
    cb.Label.String = 'Chla (mg m^{-3})'; 
    cb.FontSize = 8.5;

%     set(gcf,'PaperPositionMode','auto')
%     print(gcf,fullfile(fullPathPlotsDir,strcat('check_high_chla_5day_product_',figureSuffix)),'-dpdf','-r0')
    exportgraphics(gcf,fullfile(fullPathPlotsDir,...
        strcat('check_high_chla_5day_product_',figureSuffix,'.png')),'Resolution',600)
    
end

end % plotScenesWithMaskedValues

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotComparisonScenesBeforeAndAfterInterpolation(iTime,...
    chlaSatelliteMask,chlaSatelliteInt,latSatellite,lonSatellite,...
    timeVectorSatellite,areaStudy,enduranceSite,fullPathPlotsDir)

titleStr1 = datestr(timeVectorSatellite(iTime), 'yyyy-mm-dd');

figure()               
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.35],'Color','w') 
haxis = zeros(1,2);

for i = 1:2
    if (i == 1)
        data = chlaSatelliteMask(:,:,iTime);
        titleStr2 = 'Original data';
    elseif (i == 2)
        data = chlaSatelliteInt(:,:,iTime);
        titleStr2 = 'Interpolated data';
    end

    haxis(i) = subaxis(1,2,i,'Spacing',0.015,'Padding',0.010,'Margin',0.08);
    ax(i).pos = get(haxis(i),'Position');
    ax(i).pos(2) = ax(i).pos(2) + 0.035;
    set(haxis(i),'Position',ax(i).pos) 

    m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
    m_pcolor(lonSatellite,latSatellite,data)
    shading flat
    colormap(brewermap(1000,'*RdYlBu'))
    caxis([0 15])
    set(gca,'color','w'); % or whatever color you want for the background
    m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
    m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
    title([titleStr1, newline, titleStr2],'FontSize',11)
    box on
    hold on
end

% Colour bar settings
cb = colorbar('Location','southoutside');
cb.Position(1) = cb.Position(1) - 0.35;
cb.Position(2) = cb.Position(2) - 0.07;
cb.Position(3) = 0.60; % LENGTH
cb.Position(4) = 0.015; % WIDTH
cb.Label.String = 'Chla (mg m^{-3})'; 
cb.FontSize = 11;

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,strcat('interpolation_before_and_after_',titleStr1)),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    strcat('interpolation_before_and_after_',titleStr1,'.png')),'Resolution',600)
    
end % plotComparisonScenesBeforeAndAfterInterpolation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotCMEMSmixedLayerDepthAndEuphoticZoneDepth(mldCmemsModel,kdCmemsModel,...
    nCmemsLat,nCmemsLon,timeVectorCmemsModel,nTimeStepsCmemsModel,fullPathPlotsDir)

mldCmemsModelSceneAveraged = NaN(1,length(timeVectorCmemsModel)); % m
for iTimeStep = 1:nTimeStepsCmemsModel
    thisScene = reshape(mldCmemsModel(:,:,iTimeStep),[],1);
    mldCmemsModelSceneAveraged(iTimeStep) = mean(thisScene,'omitnan');
end

% Perform linear regression to estimate the trend
[pmld,~] = polyfit(datenum(timeVectorCmemsModel), mldCmemsModelSceneAveraged, 1);
% Generate the trend line
trendLineMld = polyval(pmld, datenum(timeVectorCmemsModel));

zeuCmemsModel = NaN(size(kdCmemsModel));
zeuCmemsModelSceneAveraged = NaN(1,length(timeVectorCmemsModel));
for iTimeStep = 1:nTimeStepsCmemsModel
    for iLon = 1:nCmemsLon
        for iLat = 1:nCmemsLat
            if (kdCmemsModel(iLat,iLon,iTimeStep) > 0 && ~isnan(kdCmemsModel(iLat,iLon,iTimeStep)))
                zeuCmemsModel(iLat,iLon,iTimeStep) = 4.6./kdCmemsModel(iLat,iLon,iTimeStep); % standard formula to calculate zeu
            end
        end
    end
    thisScene = reshape(zeuCmemsModel(:,:,iTimeStep),[],1);
    zeuCmemsModelSceneAveraged(iTimeStep) = mean(thisScene,'omitnan');
end

% Interpolate NaN value
iNan = find(isnan(zeuCmemsModelSceneAveraged));
zeuCmemsModelSceneAveraged(iNan) = interp1(...
    timeVectorCmemsModel(~isnan(zeuCmemsModelSceneAveraged)),...
    zeuCmemsModelSceneAveraged(~isnan(zeuCmemsModelSceneAveraged)),...
    timeVectorCmemsModel(iNan));

% Perform linear regression to estimate the trend
[pzeu,~] = polyfit(datenum(timeVectorCmemsModel), zeuCmemsModelSceneAveraged, 1);
% Generate the trend line
trendLineZeu = polyval(pzeu, datenum(timeVectorCmemsModel));

% Plot
startDate = min(timeVectorCmemsModel);
endDate = max(timeVectorCmemsModel);

% Identify winter months (e.g., December, January, February)
winterMonths = month(timeVectorCmemsModel) >= 12 | month(timeVectorCmemsModel) <= 2;

% Find the start and end indices of continuous winter months
winterStart = find(diff([0; winterMonths]) == 1);
winterEnd = find(diff([winterMonths; 0]) == -1);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.47],'Color','w') 

for iSubplot = 1:2
    
    haxis(iSubplot) = subaxis(2,1,iSubplot,'Spacing',0.010,'Padding',0.025,'Margin',0.10);
    ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    if (iSubplot == 1)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.025;
    elseif (iSubplot == 2)
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2);
    end
    set(haxis(iSubplot),'Position',ax(iSubplot).pos)

    if (iSubplot == 1)
        scatter(haxis(iSubplot),timeVectorCmemsModel,mldCmemsModelSceneAveraged,6,'b','filled'); hold on; 
        plot(haxis(iSubplot),timeVectorCmemsModel,trendLineMld,'r','LineWidth',1.5); hold on;
        ylabelStr = 'Mixed layer depth (m)';
    elseif (iSubplot == 2)
        scatter(haxis(iSubplot),timeVectorCmemsModel,zeuCmemsModelSceneAveraged,6,'b','filled'); hold on; 
        plot(haxis(iSubplot),timeVectorCmemsModel,trendLineZeu,'r','LineWidth',1.5); hold on;
        ylabelStr = 'Euphotic zone depth (m)';
    end

    set(gca,'Ydir','reverse')
    ylim([0 50])
    xlim([startDate endDate])
    
    % Plot a grey shaded area for winter months
    for i = 1:length(winterStart)
        fill(haxis(iSubplot),[timeVectorCmemsModel(winterStart(i)), timeVectorCmemsModel(winterStart(i)), timeVectorCmemsModel(winterEnd(i)), timeVectorCmemsModel(winterEnd(i))], ...
             [ylim, fliplr(ylim)], [0.7, 0.7, 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Winter Months');
    end

    ylabel(ylabelStr)
    if (iSubplot == 2)
        xlabel('Time')
    end
    grid on
    box on
    
end

% Create a wider title spanning both subplots
annotation('textbox',[0 0.5 1 0.5],... 
    'String','CMEMS reanalysis',...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');

% Add common legend
lg = legend('Data','Trend line');
lg.Position(1) = 0.87; lg.Position(2) = 0.835; 
lg.Orientation = 'vertical';
lg.ItemTokenSize = [10,20];
set(lg,'Box','on','Color','w')

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,'plot_cmems_mld_and_zeu'),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'plot_cmems_mld_and_zeu.png'),'Resolution',600)

end % plotCMEMSmixedLayerDepthAndEuphoticZoneDepth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotComparisonCMEMSchlaBeforeAndAfterRegridding(...
    chlaCmemsModelEuphoticOneDay7Km,latCmemsModel,lonCmemsModel,timeCmemsModel,...
    chlaCmemsModelEuphoticFiveDay4km,newLatCmemsModel,newLonCmemsModel,newTimeCmemsModel,...
    areaStudy,enduranceSite,fullPathPlotsDir)

testDate = '2006-08-24';
iTimeOriginalVector = find(timeCmemsModel == testDate);
iTimeNewVector = find(newTimeCmemsModel == testDate);

titleStr1 = datestr(timeCmemsModel(iTimeOriginalVector), 'yyyy-mm-dd');

figure()               
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.35],'Color','w') 
haxis = zeros(1,2);

for i = 1:2
    if (i == 1)
        data = squeeze(chlaCmemsModelEuphoticOneDay7Km(:,:,iTimeOriginalVector));
        lat = latCmemsModel;
        lon = lonCmemsModel;
        titleStr2 = 'Original CMEMS (daily, 7 km)';
    elseif (i == 2)
        data = squeeze(chlaCmemsModelEuphoticFiveDay4km(:,:,iTimeNewVector));
        lat = newLatCmemsModel;
        lon = newLonCmemsModel;
        titleStr2 = 'Regridded CMEMS (5-day, 4 km)';
    end

    haxis(i) = subaxis(1,2,i,'Spacing',0.015,'Padding',0.010,'Margin',0.08);
    ax(i).pos = get(haxis(i),'Position');
    ax(i).pos(2) = ax(i).pos(2) + 0.035;
    set(haxis(i),'Position',ax(i).pos) 

    m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
    m_pcolor(lon,lat,data)
    shading flat
    colormap(brewermap(1000,'*RdYlBu'))
    caxis([0 2])
    set(gca, 'color', 'w'); % or whatever color you want for the background
    m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
    m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
    title([titleStr1, newline, titleStr2],'FontSize',11)
    box on
    hold on
end

% Colour bar settings
cb = colorbar('Location','southoutside');
cb.Position(1) = cb.Position(1) - 0.35;
cb.Position(2) = cb.Position(2) - 0.07;
cb.Position(3) = 0.60; % LENGTH
cb.Position(4) = 0.015; % WIDTH
cb.Label.String = 'Chla (mg m^{-3})'; 
cb.FontSize = 11;

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,strcat('plot_regridding_before_and_after_',titleStr1)),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    strcat('plot_regridding_before_and_after_',titleStr1,'.png')),'Resolution',600)

end % plotComparisonCMEMSchlaBeforeAndAfterRegridding

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotSatelliteTimeSeriesWithGapsFilledWithCMEMSdata(...
    chlaSatelliteIntMissingFilled,isPermanentlyMissing,...
    latSatellite,lonSatellite,nSatelliteLat,nSatelliteLon,...
    timeVectorSatellite,fullPathPlotsDir)

% Randomly sample the vector chlaSatelliteCompleteFiveDay4kmIntFilled four times
nSamples = 4; 

% To focus the x axis
% startDate = datetime('1998-01-01');
% endDate = datetime('2002-01-01'); %datetime(endDateCmems)
startDate = datetime('2019-01-01');
endDate = datetime('2023-12-31'); 

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
    
    % Randomly sample
    randomIndxsLon = randi(nSatelliteLon, 1, 1);
    randomIndxsLat = randi(nSatelliteLat, 1, 1);
    dataArray = squeeze(chlaSatelliteIntMissingFilled(randomIndxsLat,randomIndxsLon,:));
    logicalArray = squeeze(isPermanentlyMissing(randomIndxsLat,randomIndxsLon,:));
    
    scatter(haxis(iSample),timeVectorSatellite(~logicalArray),dataArray(~logicalArray),...
        10,'b','filled'); hold on
    scatter(haxis(iSample),timeVectorSatellite(logicalArray),dataArray(logicalArray),...
        10,'r','filled'); hold on
    xlim([startDate endDate])
    % ylim([0 1])
    grid on;
    ylabel('Chla (mg m^{-3})','FontSize',11)
    if (iSample == 3 || iSample == 4)
        xlabel('Time (date)','FontSize',11);
    end
    title(sprintf('%2.2f°N, %2.2f°E',latSatellite(randomIndxsLat),lonSatellite(randomIndxsLon)))
    box on
end
hold off

lg = legend('Original (OC-CCI)','Fills (CMEMS)');  
lg.Orientation = 'vertical';
lg.Position(1) = 0.85; lg.Position(2) = 0.86;
lg.ItemTokenSize = [15,1];
lg.FontSize = 11;
set(lg,'Box','off')

% set(gcf,'PaperPositionMode','manual')
% orient(gcf,'landscape')
% print(gcf,fullfile(fullPathPlotsDir,'plot_data_gaps_filled'),'-dpdf','-r0') 
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'plot_data_gaps_filled.png'),'Resolution',600)

end % plotSatelliteTimeSeriesWithGapsFilledWithCMEMSdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotComparisonSmartBuoyDataWithSatelliteDataAndModelReanalysis(...
    chlaSatellite,timeVectorSatellite,chlaCmemsModelEuphotic,timeVectorCmemsModel,...
    fullPathMainDir,fullPathPlotsDir)

% Load the SmartBuoy dataset
fullPathSmartBuoyDir = strcat(fullPathMainDir,'CEFAS_SmartBuoy/');
load(strcat(fullPathSmartBuoyDir,'cefasSmartBuoy.mat')) % loads T

% Extract the chlorophyll and group it by buoy
buoyCategories = categories(T.buoy);
chlaSmartBuoy = [];
for iBuoy = 1:length(buoyCategories)
    thisBuoyChla = T(T.variable == "CHLOR" & T.buoy == buoyCategories(iBuoy),:);
    chlaSmartBuoy = [chlaSmartBuoy; thisBuoyChla];
end

% Satellite and model chlorophyll
[nSatRows,nSatCols,nSatTimeSteps] = size(chlaSatellite); 
[~,~,nModelTimeSteps] = size(chlaCmemsModelEuphotic); % num. rows and cols. coincides with satellite product
meanTimeStepSatelliteChla = zeros(1,nSatTimeSteps);
meanTimeStepModelChla = zeros(1,nModelTimeSteps);
for iTimeStep = 1:nSatTimeSteps
    thisTimeStepSatelliteChla = reshape(squeeze(chlaSatellite(:,:,iTimeStep)),[],1);
    meanTimeStepSatelliteChla(iTimeStep) = mean(thisTimeStepSatelliteChla,'omitnan');  
end
for iTimeStep = 1:nModelTimeSteps
    thisTimeStepModelChla = reshape(squeeze(chlaCmemsModelEuphotic(:,:,iTimeStep)),[],1); 
    meanTimeStepModelChla(iTimeStep) = mean(thisTimeStepModelChla,'omitnan');
end

% Plot settings
nBuoys = length(buoyCategories);
coloursBuoys = brewermap(nBuoys,'*Set1');
colourSatellite = [255 127 0]./255;
colourModel = [255 255 51]./255;
colourAllSatellitePixels = colormap(brewermap(nSatRows*nSatCols,'*paired'));

startDate = datetime('2006-01-01');
endDate = datetime('2013-09-15');

yBreakValue = 4; % mg chla m-3

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.55],'Color','w') 

for iSubplot = 1:2
        
    if (iSubplot == 1) % all satellite pixels
        
        i = 1;
        for iRow = 1:nSatRows
            for iCol = 1:nSatCols
                
                x = timeVectorSatellite;
                y = squeeze(chlaSatellite(iRow,iCol,:));
        
                % Create axes for the bottom part of the plot
                if (iRow == 1 && iCol == 1)    
                    ax1 = axes('Position', [0.1 0.55 0.85 0.25]); 
                    ax1.XGrid = 'on';
                    ax1.YGrid = 'on'; 
                end
                hold on

                scatter(ax1,x,y,4,colourAllSatellitePixels(i,:),'filled'); 
                ylim([0 yBreakValue])
                if (iRow == 1 && iCol == 1)
                    yticks([0.5 1 1.5 2 2.5 3 3.5 4])
                    yticklabels({'0.5','1.0','1.5','2.0','2.5','3.0','3.5','4.0'})
                end
                xlim([startDate endDate])
                grid on
                box on
                hold on
               
                % Create axes for the upper part of the plot
                if (iRow == 1 && iCol == 1)
                    ax2 = axes('Position', [0.10 0.80 0.85 0.15]); 
                    ax2.XAxis.Visible = 'off';
                    ax2.YScale = 'log';
                    ax2.XGrid = 'on';
                    ax2.YGrid = 'on'; 
                end
                hold on
                
                scatter(ax2,x,y,4,colourAllSatellitePixels(i,:),'filled'); 
                ylim([yBreakValue 10])
                if (iRow == 1 && iCol == 1)
                    yticks([10])
                    yticklabels({'10'})
                end
                xlim([startDate endDate])
                grid on
                box on
                hold on
                
                i = i + 1;
                
            end
        end

        xlabel('Time')
        yl = ylabel('Chla (mg m^{-3})');
        yl.Position(1) = yl.Position(1) - 50;
        yl.Position(2) = 1.5; 
        title('OC-CCI 5-day, 4 km resolution chla time-series for each pixel in an image')
        hold on
        
    elseif (iSubplot == 2) % comparison SmartBuoy, satellite and model reanalysis

        ax3 = axes('Position', [0.10 0.070 0.85 0.40]);

        for iBuoy = 1:nBuoys
            thisBuoyChla = find(chlaSmartBuoy.buoy == buoyCategories(iBuoy));
            x = chlaSmartBuoy.dateTime(thisBuoyChla);
            y = chlaSmartBuoy.value(thisBuoyChla);
            if (isempty(y))
%                 scatter(haxis(iSubplot),NaN,NaN,6,'Color',coloursBuoys(iBuoy,:),'filled'); hold on;
                scatter(ax3,NaN,NaN,6,coloursBuoys(iBuoy,:),'filled'); hold on;
            else
%                 scatter(haxis(iSubplot),x,y,6,'Color',coloursBuoys(iBuoy,:),'filled'); hold on;
                scatter(ax3,x,y,6,coloursBuoys(iBuoy,:),'filled'); hold on;
            end
        end
        scatter(ax3,timeVectorSatellite,meanTimeStepSatelliteChla,6,colourSatellite,'filled'); hold on
        scatter(ax3,timeVectorCmemsModel,meanTimeStepModelChla,6,colourModel,'filled'); hold on
        ylim([0 15])
        xlim([startDate endDate])
        grid on
        box on
        xlabel('Time')
        yl = ylabel('Chla (mg m^{-3})');
        yl.Position(1) = yl.Position(1) - 50;
        legend([buoyCategories; 'Satellite, 5-day, 4 km (image mean)'; 'CMEMS model (image mean)'],'Location','northwest');
        title('Comparison of SmartBuoy data with satellite and model reanalysis data')

    end
    
end

% set(gcf,'PaperPositionMode','auto')
% print(gcf,fullfile(fullPathPlotsDir,'plot_comparison_satellite_model_SmartBuoy'),'-dpdf','-r0')
exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'plot_comparison_satellite_model_SmartBuoy.png'),'Resolution',600)

end % plotComparisonSmartBuoyDataWithSatelliteDataAndModelReanalysis
