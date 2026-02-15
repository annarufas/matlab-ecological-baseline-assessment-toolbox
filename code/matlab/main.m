
% ======================================================================= %
%                                                                         %
% This script uses Matlab codes to analyse XXX data from XXX place. 
% The aim is to determine what proportion of phytoplankton variability can
% be attributed to XXX. Running this code will require MATLAB_R2021a or
% higher. The script has 7 sections:                                      %
%   Section 1 - Presets.                                                  %
%   Section 2 - Load the daatset and manipulate the data array.           %
%   Section 3 - Bin data monthly by depth horizon and propagate error.    %
%   Section 4 - Bin data monthly by unique depth and propagate error.     %
%   Section 5 - Bin data annually and propagate error.                    %
%   Section 6 - Calculate the number of data points based on various      % 
%               criteria.                                                 %
%   Section 7 - Save the data.                                            %              
%                                                                         %
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD                             %
%   Anna.RufasBlanco@earth.ox.ac.uk                                       %
%                                                                         %
%   Version 1.0 - Completed XX April 2024                                  %
%                                                                         %
% ======================================================================= %

close all; clear all; clc
addpath(genpath('./code/matlab/'));
addpath(genpath('./code/internal/'));
addpath(genpath('./data/raw/'));

% =========================================================================
%%
% -------------------------------------------------------------------------
% SECTION 1 - PRESETS
% -------------------------------------------------------------------------

% Filename and directory definitions

filenameCefasHplcRawData           = 'cefas_HPLC_2010_2019.csv';
pathBathymetryFile                 = fullfile('.','data','raw','ETOPO','ETOPO_2022_v1_30s_N90W180_surface.nc');
pathMimocMldDir                    = fullfile('.','data','raw','MIMOC');

filenamesCefasSmartBuoyData        = {'Dowsing_all.csv','NorthDogger_all.csv','OysterGrounds_all.csv'};
smartBuoyNames                     = {'Dowsing','NorthDogger','OysterGround'};
filenameSsbProgramPhData           = 'SSB_Cefas_updated_final_results.xlsx';
filenameUkoaProgramPhData          = 'UKOA_final_results_corrected.xlsx';
filenameSedimentCoreData           = 'data_cores_Malini_Grab_list_for_bks_42_25_and_43_21RRSamples.xlsx';
filenameCefasCruiseData            = 'C8611 Oxford Biogeochemistry Results v2.xlsx';

filenameCmemsMatchupsPh            = fullfile('data_matchups_pH_CMT_csv','cmems_ph_matchups.csv');
filenameCmemsMatchupsHplc          = fullfile('data_matchups_HPLC_CMT_csv','cmems_hplc_matchups.csv');
filenameNasaMatchupsHplc           = fullfile('data_matchups_HPLC_NASA_csv','nasa_hplc_matchups.csv');
filenameOccciMatchupsHplc          = fullfile('data_matchups_HPLC_OCCCI_csv','occci_hplc_matchups.csv');

filenameEnduranceCentralCoords     = 'endurance_centre.mat';
shapefileEndurance                 = 'endurance_2023';
shapefileAreaStudy                 = 'bbox_100km';

pathEnduranceShapefileDir          = fullfile('.','data','raw','coords_endurance');
pathEnduranceShapefile             = fullfile(pathEnduranceShapefileDir,shapefileEndurance);
pathAreaStudyShapefileDir          = fullfile('.','data','raw','coords_boundarybox');
pathAreaStudyShapefile             = fullfile(pathAreaStudyShapefileDir,shapefileAreaStudy);
pathOsparRegionsShapefile          = fullfile('.','data','raw','coords_ospar_regions_2017_01_002/ospar_regions_2017_01_002');
pathOsparSubregionsShapefile       = fullfile('.','data','raw','coords_ospar_subregions/OSPAR_subregions_20160418_3857');
pathOsparMpasShapefile             = fullfile('.','data','raw','coords_ospar_mpas/ospar_mpa_2021_07_001');
pathUkOffshoreGcsLicensesShapefile = fullfile('.','data','raw','coords_UKCS_Carbon_Storage_Licences_(ED50)/UKCS_Carbon_Storage_Licences_(ED50)');

filenameBicepAreaStudyDataProcessed = 'bicepBbox100km.mat';
filenameCmemsAreaStudyDataProcessed = 'cmemsBbox100km.mat';
filenameNasaAreaStudyDataProcessed  = 'nasaBbox100km.mat';
filenameOccciAreaStudyDataProcessed = 'occciBbox100km.mat';
filenameCarbonStocksAreaStudy       = 'bluecarbonBbox100km.mat';

% BICEP definitions
bicepTimeSeriesDatasetNames = {...
    'BICEP_POC_nc',    {'bicep_poc_4km'};
    'BICEP_NPP_nc',    {'bicep_npp_9km'};
    'BICEP_Cphyto_nc', {'bicep_cphyto_9km'};
};
bicepYearsVector = 1998:2020; % defines a common period for all datasets, will change with new product release versions

% CMEMS definitions
dirCmemsTimeseriesAreaStudyDataRaw = 'data_timeseries_areastudy_CMT_nc';
dirCmemsTimeseriesSmartBuoyDataRaw = 'data_timeseries_smartbuoy_CMT_nc';

cmemsDatasetNames = {...
    'obs_satell_glob_cmems_olci_4km_plk',...
    'obs_satell_glob_cmems_olci_4km_trns',...
    'obs_satell_reg_cmems_multi_1km_plk',...
    'obs_satell_reg_cmems_multi_1km_opt',...
    'obs_satell_reg_cmems_multi_1km_trns',...
    'obs_satell_reg_cmems_olci_300m_plk',...
    'mod_bgc_reg_chl',...
    'mod_bgc_reg_diat',...
    'mod_bgc_reg_dino',...
    'mod_bgc_reg_nano',...
    'mod_bgc_reg_pico',...
    'mod_bgc_reg_phy',...
    'mod_bgc_reg_npp',...
    'mod_bgc_reg_kd',...
    'mod_bgc_reg_no3',...
    'mod_bgc_reg_po4',...
    'mod_bgc_reg_o2',...
    'mod_bgc_reg_ph',...
    'mod_bgc_reg_pco2',...
    'mod_phy_reg_mld',...
    'mod_phy_reg_sal',...
    'mod_phy_reg_temp',...
    'mod_phy_reg_ssh',...
    'mod_phy_reg_velo'...
};

% NASA definitions
dirNasaTimeseriesAreaStudyDataRaw = {...
    'data_timeseries_areastudy_NASA_aquamodis_nc',  {'aquamodis_4km'};
    'data_timeseries_areastudy_NASA_meris_nc',      {'meris_4km'};
    'data_timeseries_areastudy_NASA_viirssnpp_nc',  {'viirssnpp_4km'};
    'data_timeseries_areastudy_NASA_viirsjpss1_nc', {'viirsjpss1_4km'};
};

% OC-CCI definitions
dirOccciTimeseriesAreaStudyDataRaw = 'data_timeseries_areastudy_OCCCI_nc';
occciAreaStudyDatasetNames = {...
    'occci_1km_1day', {'occci_1km_1day_chl_9710.nc',...
                       'occci_1km_1day_chl_1024.nc',...
                       'occci_1km_1day_waterclass_01to03_9710.nc',...
                       'occci_1km_1day_waterclass_01to03_1024.nc',...
                       'occci_1km_1day_waterclass_04to06_9710.nc',...
                       'occci_1km_1day_waterclass_04to06_1024.nc',...
                       'occci_1km_1day_waterclass_07to09_9710.nc',...
                       'occci_1km_1day_waterclass_07to09_1024.nc',...
                       'occci_1km_1day_waterclass_10to12_9710.nc',...
                       'occci_1km_1day_waterclass_10to12_1024.nc',...
                       'occci_1km_1day_waterclass_13to14_9710.nc',...
                       'occci_1km_1day_waterclass_13to14_1024.nc'};
    'occci_4km_1day', {'occci_4km_1day_chl_9710.nc',...
                       'occci_4km_1day_chl_1024.nc'};
    'occci_4km_5day', {'occci_4km_5day_chl_9710.nc',...
                       'occci_4km_5day_chl_1024.nc'};
    'occci_4km_8day', {'occci_4km_8day_chl_9710.nc',...
                       'occci_4km_8day_chl_1024.nc'};
};

% Parameters definitions
RADIUS_AREA_STUDY = 50; % km
AVG_EUPHOTIC_ZONE_DEPTH = 30; % m, estimated from CMEMS reanalysis  for our area of study
AVG_SHELF_DEPTH = 200; % m

%%
% -------------------------------------------------------------------------
% SECTION 2 - DEFINE THE FOOTPRINT OF OUR AREA OF STUDY
%
% Create a shapefile containing the coordinates of our study area, centered
% around the Endurance GCS site, with a radius extending 50 km in all 
% directions.
% -------------------------------------------------------------------------

createAreaStudyShapefile(RADIUS_AREA_STUDY,pathUkOffshoreGcsLicensesShapefile,...
    pathEnduranceShapefileDir,pathAreaStudyShapefileDir,filenameEnduranceCentralCoords,...
    shapefileEndurance,shapefileAreaStudy)

%%
% -------------------------------------------------------------------------
% SECTION 3 - DOWNLOAD TIME-SERIES DATA FOR OUR AREA OF STUDY
%
% Time-series data is from Copernicus Marine Service (CMEMS), OC-CC, NASA
% and BICEP.
% -------------------------------------------------------------------------

isCmemsAreaStudyDataReady = 1;
isBicepAreaStudyDataReady = 1;
isNasaAreaStudyDataReady  = 1;
isOccciAreaStudyDataReady = 1;

% Running the following code requires to have run 
% ./code/jupyter/downloadCMEMStimeseriesCMT.ipynb
if ~isCmemsAreaStudyDataReady
    [AScmems] = ncreadTimeseriesCmemsData(dirCmemsTimeseriesAreaStudyDataRaw,...
        cmemsDatasetNames,filenameCmemsAreaStudyDataProcessed);
else
    load(fullfile('.','data','processed',filenameCmemsAreaStudyDataProcessed),'cmemsData')
    AScmems = cmemsData;
    clear cmemsData
end

% Running the following code requires to have run 
% ./code/jupyter/downloadBICEPtimeseries.ipynb
if ~isBicepAreaStudyDataReady
    [ASbicep] = ncreadTimeseriesBicepData(bicepTimeSeriesDatasetNames,...
        bicepYearsVector,pathAreaStudyShapefile);
    [ASbicep] = addDepthIntegrationProductsToBicepArray(ASbicep,pathMimocMldDir,...
        pathAreaStudyShapefile,AScmems,filenameBicepAreaStudyDataProcessed);
else
    load(fullfile('.','data','processed',filenameBicepAreaStudyDataProcessed),'bicepData')
    ASbicep = bicepData;
    clear bicepData
end
    
% Running the following code requires to have run
% ./code/jupyter/downloadNASAtimeseries.ipynb
if ~isNasaAreaStudyDataReady
    [ASnasa] = ncreadTimeseriesNasaData(dirNasaTimeseriesAreaStudyDataRaw,...
        filenameNasaAreaStudyDataProcessed);
else
    load(fullfile('.','data','processed',filenameNasaAreaStudyDataProcessed),'nasaData')
    ASnasa = nasaData;
    clear nasaData
end

% Running the following code requires to have run 
% ./code/jupyter/downloadOCCCItimeseries.ipynb
if ~isOccciAreaStudyDataReady
    [ASoccci] = ncreadTimeseriesOccciData(dirOccciTimeseriesAreaStudyDataRaw,...
        occciAreaStudyDatasetNames,filenameOccciAreaStudyDataProcessed,pathAreaStudyShapefile);
else
    load(fullfile('.','data','processed',filenameOccciAreaStudyDataProcessed),'occciData')
    ASoccci = occciData;
    clear occciData
end

%%
% -------------------------------------------------------------------------
% SECTION 4 - CALCULATE AND PLOT WATER COLUMN CARBON STOCKS
%
% Precise regridding to a common grid and data gap filling is necessary for
% all products.
% -------------------------------------------------------------------------

isCalculateCarbonProducts = 1;

if isCalculateCarbonProducts
    isCheckRegridding = 1;
    isCheckSceneAreaCalc = 1;
    isCheckIntegrationCalc = 1;
    isPlotHistogramsAverageStocks = 1;
    isPlotScenesAverageStocks = 1;
    [AScarbon] = calculateAndPlotCarbonStocksFromOceanColourProducts(ASbicep,AScmems,...
        bicepTimeSeriesDatasetNames,pathAreaStudyShapefile,pathEnduranceShapefile,...
        filenameCarbonStocksAreaStudy,isCheckRegridding,isCheckSceneAreaCalc,...
        isCheckIntegrationCalc,isPlotHistogramsAverageStocks,isPlotScenesAverageStocks);
end

% Retrieve ID field names
populatedIDs = ~cellfun(@isempty, {AScarbon.ID});
nCarbonDatasets = sum(populatedIDs);
carbonDatasetNames = cell(nCarbonDatasets,1);
for i = 1:nCarbonDatasets
    carbonDatasetNames{i} = AScarbon(i).ID;
end

% Calculate NPP stock range
prodName = 'bicep_npp_9km';
idxProd = find(strcmp({AScarbon.ID}, prodName));
monthlyNppStocks = AScarbon(idxProd).sceneStockDistribByMonthAndYear_monthlyMean; % 12 months x N years
annualNppStocks = sum(monthlyNppStocks,1); % Gg C yr-1
annualNppStocks_Mt = annualNppStocks.*1e-3; % Mt C yr-1

annualNppStocks_Mt_cum = sum(annualNppStocks_Mt).*3.67; % Mt CO2 yr-1

fprintf('\nThe bounds for integrated NPP stocks are %4.1f to %4.1f Gg C/month',min(monthlyNppStocks),max(monthlyNppStocks))
fprintf('\nThe average for integrated NPP stocks is %4.1f Mt C yr-1',mean(annualNppStocks_Mt))
fprintf('\nThe bounds for integrated NPP stocks are %4.1f to %4.1f Mt C yr-1',min(annualNppStocks_Mt),max(annualNppStocks_Mt))
fprintf('\nIt is estimated that from 1998 to 2020, phytoplankton at the Endurance have (cumulatively) removed the equivalent of %4.1f million tonnes of carbon dioxide (MtCO2e) from the atmosphere',annualNppStocks_Mt_cum)

%%
% -------------------------------------------------------------------------
% SECTION 5 - EXPLORE THE PRODUCTS DOWNLOADED
%
% 
% -------------------------------------------------------------------------

isPlotHovmollerDiagrams   = 1;
isPlotAnomalyPlots        = 1;
isPlotOpticalWaterTypes   = 1;
isPlotPftProducts         = 1;
isPlotChlorophyllProducts = 1;

datasetConfig = struct(...
    'unitsVar', {'mg C m^{–3}', 'mg C m^{–3}', 'mg C m^{–3}', 'mg C m^{–3}',...
                 'm', 'mmol m^{-3}','mmol m^{-3}', 'mmol m^{-3}',... 
                 '', 'Pa', 'm', 'PSU',...
                 'ºC','m', 'm s^{-1}', 'm s^{-1}',...
                 'mg C m^{–2} d^{-1}', 'mg C m^{–3}', 'mg C m^{–3}', 'mg C m^{–3}',...
                 'mg C m^{–3}', 'mg C m^{–3}',...
                 'mg chla m^{–3}', 'mg chla m^{–3}', 'mg chla m^{–3}', 'mg chla m^{–3}',...
                 'mg chla m^{–3}', 'mg chla m^{–3}', 'mg chla m^{–3}', 'mg chla m^{–3}'},...    
    'tickLabelFmt', {'%.1f', '%.2f', '%.2f', '%.2f',...
                     '%.0f', '%.1f', '%.2f', '%.0f',...
                     '%.2f', '%.0f', '%.0f', '%.1f', ...
                     '%.0f', '%.1f', '%.2f', '%.2f',...
                     '%.0f', '%.0f', '%.0f', '%.0f',...
                     '%.0f', '%.0f',...
                     '%.1f', '%.1f', '%.1f', '%.1f',...
                     '%.1f', '%.1f', '%.1f', '%.1f'},...
    'minVar', {0, 0, 0, 0,... 
               20, 0.5, 0, 220,... 
               8.01, 30, 0, 33.9, ...
               2, -1, -0.20, -0.20, ...
               0, 0, 0, 0,... 
               0, 0,...
               0, 0, 0, 0,... 
               0, 0, 0, 0},...
    'maxVar', {0.5, 0.02, 0.2, 0.15,...
               50, 6, 0.4, 300,...
               8.17, 43, 50, 35, ...
               20, 0.1, 0.15, 0.15,...
               1000, 400, 80, 25,...
               20, 20,...
               5, 5, 5, 5,...
               5, 5, 5, 5},...
    'maxPercent', {50, 50, 50, 25,...
               2, 50, 50, 3,...
               1, 15, 15, 1, ...
               10, 20, 50, 50,...
               10, 20, 10, 20,...
               10, 10,...
               40, 40, 40, 20,...
               20, 40, 40, 40});

if isPlotHovmollerDiagrams           
    plotHovmollerDiagramsForDrivers(AScmems,ASbicep,ASnasa,ASoccci,datasetConfig)
end

if isPlotAnomalyPlots
    plotAnomalyPlotsForDrivers(AScmems,ASbicep,ASnasa,ASoccci,datasetConfig)    
end

if isPlotOpticalWaterTypes
    yearOfChoice = [2000,2010,2020];
    for i = 1:numel(yearOfChoice)
        plotOpticalWaterTypes(yearOfChoice(i),ASoccci,pathAreaStudyShapefile,...
            pathEnduranceShapefile)
    end
end

if isPlotPftProducts
    plotPftsFromOceanColourProducts(AScmems,ASbicep)
end

if isPlotChlorophyllProducts
    plotChlorophyllFromOceanColourProducts(AScmems,ASnasa,ASoccci,...
        pathAreaStudyShapefile,pathEnduranceShapefile)
end

%%
% -------------------------------------------------------------------------
% SECTION 6 - PREPARE, VISUALISE AND ANALYSE IN SITU HPLC DATA
%
% HPLC data were downloaded from CEFAS (https://data.cefas.co.uk/view/53).
% The dataset runs from 2010–2019. This section filters out some data 
% according to euphotic zone depth and shelf depth, plots it on a map
% and calculates phytoplankton size classes.
%
% -------------------------------------------------------------------------

% Define which actions get executed
isHplcDataReady                            = 0;
isPlotHplcData                             = 1;
isCalculateAndPlotPhytoplanktonSizeClasses = 1;
isPerformHplcMatchupAnalysis               = 1;

% Prepare/load the CEFAS HPLC data
if ~isHplcDataReady
    [HPLC] = prepareCefasHplcData(filenameCefasHplcRawData,AVG_EUPHOTIC_ZONE_DEPTH,...
        AVG_SHELF_DEPTH,pathBathymetryFile);
else
    load(fullfile('.','data','processed','cefasHPLCfiltered.mat'),'HPLC')
end

% Plot the CEFAS HPLC data
if isPlotHplcData
    plotCefasHplcData(pathEnduranceShapefile,pathAreaStudyShapefile,HPLC)
end

% Calculate fractional contribution of three main phytoplankton size
% classes to the HPLC data and make various plots
if isCalculateAndPlotPhytoplanktonSizeClasses 
    calculateAndPlotPhytoplanktonSizeClasses(HPLC)
end

% Running this function requires to have run:
% ./code/jupyter/downloadCMEMSmatchupsCMT.ipynb
% ./code/jupyter/downloadNASAmatchups.ipynb
% ./code/jupyter/downloadOCCCImatchups.ipynb
if isPerformHplcMatchupAnalysis
    matchupAnalysisHplcData(filenameCmemsMatchupsHplc,filenameNasaMatchupsHplc,...
        filenameOccciMatchupsHplc,HPLC,AScmems,ASnasa,ASoccci,pathEnduranceShapefile,...
        pathAreaStudyShapefile)
end

%%
% -------------------------------------------------------------------------
% SECTION 7 - BATHYMETRIC MAP OF THE NORTH SEA WITH SMARTBUOY, HPLC DATA
% AND SHIP TRANSECT DATA
%
% Plot a bathymetric map of the North Sea with the locations of in situ
% sampling
% -------------------------------------------------------------------------

isPlotNorthSeaBathymetricMap = 1;

if isPlotNorthSeaBathymetricMap 
    plotNorthSeaBathymetricMap(pathOsparRegionsShapefile,...
        pathOsparSubregionsShapefile,pathOsparMpasShapefile,...
        pathUkOffshoreGcsLicensesShapefile,pathAreaStudyShapefile,...
        HPLC,filenameSedimentCoreData,filenameCefasCruiseData)
end

%%
% -------------------------------------------------------------------------
% SECTION 8 - PREPARE IN SITU AND REANALYSIS DATASETS HANDED IN BY
% RESEARCHERS AT CEFAS
%
% -------------------------------------------------------------------------

isInsituPhDataReady            = 1;
isPlotInsituPhData             = 1;
isPerformMatchupAnalysisPhData = 1;
isSilvaSpmDataReady            = 0;

% Prepare the pH dataset
if ~isInsituPhDataReady
   [PH] = prepareInsituPhData(filenameSsbProgramPhData,filenameUkoaProgramPhData,...
        AVG_SHELF_DEPTH,pathBathymetryFile);
else
    load(fullfile('.','data','processed','greenwoodPhFiltered.mat'),'PH')
end

% Plot the pH data
if isPlotInsituPhData
    plotInsituPhData(pathEnduranceShapefile,pathAreaStudyShapefile,PH)
end

% Matchup analysis pH data with CMEMS biogeochemical reanalysis
if isPerformMatchupAnalysisPhData
    matchupAnalysisPhData(filenameCmemsMatchupsPh,PH)
end

% Prepare the SPM dataset
if ~isInsituPhDataReady
    [SPMclimatol] = prepareSilvaSpmData(pathAreaStudyShapefile);
else
    load(fullfile('.','data','processed','silvaSPM.mat'),'SPMclimatol','latSpm','lonSpm')
end

% Ship transect data...

%%
% -------------------------------------------------------------------------
% SECTION 5 - COMPARE IN SITU DATA WHERE THE SMARTBUOY SYSTEMS ARE WITH 
% PRODUCT DATA 
%
% bla bla bla
%
% -------------------------------------------------------------------------

% Define which actions get executed
isCefasSmartBuoyDataReady    = 1;
isCmemsSmartBuoyDataReady    = 0;
isPlotCefasSmartBuoyData     = 0;
isPlotComparisonCefasVsCmemsSmartBuoy = 0;

% Prepare/load the CEFAS SmartBuoy dataset
if ~isCefasSmartBuoyDataReady 
    [SB] = prepareCefasSmartBuoyData(filenamesCefasSmartBuoyData);
else
    load(fullfile('.','data','processed','cefasSmartBuoy.mat'),'SB')
end

% Plot the CEFAS SmartBuoy data
if isPlotCefasSmartBuoyData
    plotCefasSmartBuoyData(SB)
end

% Read .nc files downlaoded from the Copernicus Marine Data Store using the
% Copernicus Marine Toolbox (CMT). Last download: 20 April 2024.
if ~isCmemsSmartBuoyDataReady
    SBcmems = struct();
    for iSmartBuoy = 1:length(smartBuoyNames)
        thisBuoy = smartBuoyNames{iSmartBuoy};
        dirThisBuoy = fullfile(dirCmemsTimeseriesSmartBuoyDataRaw,thisBuoy);
        outputFilenameThisBuoy = strcat('cmemsSmartBuoy',thisBuoy,'.mat');
        [cmemsData] = ncreadTimeseriesCmemsData(dirThisBuoy,...
            cmemsDatasetNames,outputFilenameThisBuoy);
        SBcmems.(sprintf('%s', thisBuoy)) = cmemsData;
        clear cmemsData
    end
else
    SBcmems = struct();
    for iSmartBuoy = 1:length(smartBuoyNames)
        load(fullfile('.','data','processed',strcat('cmemsSmartBuoy',thisBuoy,'.mat')),'cmemsData')
        SBcmems.(sprintf('%s',thisBuoy)) = cmemsData;
        clear cmemsData
    end
end


if isPlotComparisonCefasVsCmemsSmartBuoy
    
end





