% function [AScarbon] = calculateAndPlotCarbonStocksFromOceanColourProducts(...
%     ASbicep,AScmems,bicepTimeSeriesDatasetNames,pathAreaStudyShapefile,...
%     pathEnduranceShapefile,filenameCarbonStocksAreaStudy,isCheckRegridding,...
%     isCheckSceneAreaCalc,isCheckIntegrationCalc,isPlotHistogramsAverageStocks,...
%     isPlotScenesAverageStocks)

% CALCULATEANDPLOTCARBONSTOCKSFROMOCEANCOLOURPRODUCTS Calculate carbon
% stocks products for the Endurance site. There are options to turn on/off 
% the possibility of plotting various properties.
%
%   INPUT:
%       ASbicep                       - table with time-series data from BICEP for our area of study
%       AScmems                       - table with time-series data from CMEMS for our area of study
%       bicepTimeSeriesDatasetNames   - string containing the names of the BICEP datasets
%       pathAreaStudyShapefile        - shapefile with our area of study
%       pathEnduranceShapefile        - shapefile with the Endurance GCS site
%       filenameCarbonStocksAreaStudy - .mat file containing AScarbon
%       isCheckRegridding             - choice for plotting
%       isCheckSceneAreaCalc          - choice for plotting
%       isCheckIntegrationCalc        - choice for plotting
%       isPlotHistogramsAverageStocks - choice for plotting
%       isPlotScenesAverageStocks     - choice for plotting
%
%   OUTPUT:
%       AScarbon                      - Matlab table with the blue carbon data processed
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 4 May 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

fprintf("\nReading carbon products products...")

%% Presets

% Conversion factors
% 1 Gg C = 1000 t C
% 1 Tg C = 1 Mt C
MILIGRAM_TO_TERAGRAM = 1e-3*1e-12; % constant to convert from mg C to Tg C 
MILIGRAM_TO_GIGAGRAM = 1e-3*1e-9; % constant to convert from mg C to Gg C
MILIGRAM_TO_GRAM = 1e-3; % constant to convert from mg C to g C
INTEGRATION_TIME_NPP = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]; % number of days in each month

% Load Endurance area
areaStudy = m_shaperead(pathAreaStudyShapefile); 
enduranceSite = m_shaperead(pathEnduranceShapefile);

% Initialise structure
AScarbon = struct('ID', {},...
                  'carbonDataArray', {},...
                  'carbonDataArrayRegriddedToAreaStudy', {},...
                  'carbonDataArrayInt', {},...
                  'latVectorCarbon', {},...
                  'lonVectorCarbon', {},...
                  'latVectorAreaStudy', {},... 
                  'lonVectorAreaStudy', {},...            
                  'timeVectorCarbon', {},...
                  'depthIntegrationDataArray', {},...
                  'depthIntegrationDataArrayRegriddedToAreaStudy', {},...
                  'depthIntegrationDataArrayInt', {},...
                  'latVectorDepthIntegrationArray', {},...
                  'lonVectorDepthIntegrationArray', {},...
                  'depthVectorCarbon', {},... 
                  'nResolvedDepthLevels', {},...
                  'depthResolutionCarbon', {},...
                  'areaPixel', {},...
                  'sceneArea', {},...
                  'pixelConcentrationPerUnitArea', {},...     
                  'pixelStock', {},...                    
                  'pixelWaterColumnVolume', {},...           
                  'sceneConcentrationPerUnitArea_lit', {},...  
                  'sceneConcentrationPerUnitArea_mine', {},... 
                  'sceneStock_lit', {},...                     
                  'sceneStock_mine', {},...                    
                  'sceneWaterVolume', {});                  

%% Load data and regrid to area of study

% The loaded carbon data have units of mg C, data have been regridded to our 
% area of study and the interpolation method used to fill in the data gaps 
% does not affect depths that should have no data.

[AScarbon] = processCmemsCphytoDataset(AScmems,AScarbon,areaStudy);
[AScarbon] = processBicepDatasets(ASbicep,bicepTimeSeriesDatasetNames,AScarbon,areaStudy);

% Retrieve ID field names
populatedIDs = ~cellfun(@isempty, {AScarbon.ID});
nCarbonDatasets = sum(populatedIDs);
carbonDatasetNames = cell(nCarbonDatasets,1);
for i = 1:nCarbonDatasets
    carbonDatasetNames{i} = AScarbon(i).ID;
end

%% Checks

% After interpolating, it is noticeable that there are still persistent 
% missing data points.

if isCheckRegridding
    % Plot before and after regridding & interpolating
    testDate = '2019-04-15';
    plotComparisonBeforeAndAfterRegridding(AScarbon,carbonDatasetNames,testDate,areaStudy,enduranceSite)
    % Plot the percentage missing data after regriding and interpolating
    plotScenePercentageMissingDataByYear(AScarbon,carbonDatasetNames,areaStudy,enduranceSite)
    plotScenePercentageMissingDataByMonth(AScarbon,carbonDatasetNames,areaStudy,enduranceSite)
    % Plot the products used to integrate carbon stocks in the water column
    plotComparisonDepthIntegrationProducts(AScmems,AScarbon)
end

%% Calculate the area of each pixel in the new grid

% Showing two approaches: Matlab's package from the Climate Data Toolbox vs
% Matlab's package from the Mapping Toolbox

for iDataset = 1:nCarbonDatasets
    
    prodName = carbonDatasetNames{iDataset};
    idxProd = find(strcmp({AScarbon.ID}, prodName));
    fprintf('\nCalculating area of a pixel for product %s...',prodName)
    
    % Get regridded latitude and longitude vectors
    latVectorRegridded = AScarbon(idxProd).latVectorAreaStudy;
    lonVectorRegridded = AScarbon(idxProd).lonVectorAreaStudy;
    
    % Calculate pixel area using cdtarea from the Climate Data Toolbox
    [latGrid,lonGrid] = ndgrid(latVectorRegridded,lonVectorRegridded);
    areaPixel = cdtarea(latGrid,lonGrid); 
    sceneArea = sum(reshape(areaPixel(:,:),[],1),'omitnan'); % m2
    fprintf('\n...our region of study covers an area of %.0f km2',sceneArea*1e-6)

    % ... or using wgs84Ellipsoid from the Mapping Toolbox 
    % wgs84 = wgs84Ellipsoid; % requires mapping toolbox (this accounts for the Earth's curvature)
    % areaPixel = zeros(numel(latVectorRegridded),numel(lonVectorRegridded));
    % for iLat = 1:numel(latVectorRegridded)
    % %     thisLat = latVectorRegridded(iLat) - 0.5*latResolution;
    % %     nextLat = latVectorRegridded(iLat) + 0.5*latResolution;
    %     thisLat = latVectorRegridded(iLat);
    %     if (iLat == numel(latVectorRegridded))
    %         nextLat = latVectorRegridded(iLat) + latResolution;
    %     else
    %         nextLat = latVectorRegridded(iLat+1);
    %     end
    %     for iLon = 1:numel(lonVectorRegridded)
    % %         thisLon = lonVectorRegridded(iLat) - 0.5*lonResolution;
    % %         nextLon = lonVectorRegridded(iLat) + 0.5*lonResolution;
    %         thisLon = lonVectorRegridded(iLon);
    %         if (iLon == numel(lonVectorRegridded))
    %             nextLon = lonVectorRegridded(iLon) + lonResolution;
    %         else
    %             nextLon = lonVectorRegridded(iLon+1);
    %         end
    %         [xNorth_m,yEast_m,~] = geodetic2ned(nextLat,nextLon,1,...
    %             thisLat,thisLon,1,wgs84,'degrees'); % m
    %         areaPixel(iLat,iLon) = yEast_m*xNorth_m; % m2
    %     end
    % end 
    % sceneArea = sum(reshape(areaPixel(:,:),[],1),'omitnan'); % m2
    % disp(['Our region of study covers an area of ',num2str(sceneArea*1e-6,'%.0f'),' km2'])

    % Save into output array
    AScarbon(idxProd).areaPixel = areaPixel; % m2
    AScarbon(idxProd).sceneArea = sceneArea; % m2
    
    % Check that the area of a pixel looks as it should
    if isCheckSceneAreaCalc
        plotAreaOfPixel(areaPixel*1e-6,latVectorRegridded,lonVectorRegridded,...
            areaStudy,enduranceSite,sprintf('Scene area: %.0f km^{2}',sceneArea*1e-6),prodName)
    end
    
end % iDataset

% Notice how the areas of the scenes vary as a function of the latitude and
% longitude vectors used

%% Integrate quantities over depth and, for NPP, over time too

% Month timer for NPP calculations
iMonth = 1;

for iDataset = 1:numel(carbonDatasetNames)
    
    prodName = carbonDatasetNames{iDataset};
    idxProd = find(strcmp({AScarbon.ID}, prodName));

    nLatPixels = numel(AScarbon(idxProd).latVectorAreaStudy);
    nLonPixels = numel(AScarbon(idxProd).lonVectorAreaStudy);
    nTimeSteps = numel(AScarbon(idxProd).timeVectorCarbon);
    
    pixelConcentrationPerUnitArea      = NaN(nLatPixels,nLonPixels,nTimeSteps); % mg m-2 (depth-integrated)
    pixelStock                         = NaN(nLatPixels,nLonPixels,nTimeSteps); % mg (depth-integrated)
    pixelWaterColumnVolume             = NaN(nLatPixels,nLonPixels,nTimeSteps); % m3
    sceneConcentrationPerUnitArea_lit  = NaN(nTimeSteps,1); % mg m-2 (approach presented in Allison et al. 2010, which considers missing, invalid pixels)
    sceneConcentrationPerUnitArea_mine = NaN(nTimeSteps,1); % mg m-2 (my approach)
    sceneStock_lit                     = NaN(nTimeSteps,1); % mg (approach presented in Allison et al. 2010, which considers missing, invalid pixels)
    sceneStock_mine                    = NaN(nTimeSteps,1); % mg (my approach)
    sceneWaterVolume                   = NaN(nTimeSteps,1); % m3 

    % NOTE: this time step, for NPP, is monthly (tye first day of each
    % month). Thus, we can multiply each value by the number of days in the 
    % respective month to get the total monthly production.
    for iTimeStep = 1:nTimeSteps
        for iLon = 1:nLonPixels
            for iLat = 1:nLatPixels
                
                thisPixelArea = AScarbon(idxProd).areaPixel(iLat,iLon); % m2
            
                if (strcmp(prodName,'cmems_cphyto')) 
                    
                    iLastResolvedDepth = AScarbon(idxProd).nResolvedDepthLevels(iLat,iLon,iTimeStep);
                    waterColumnCarbon = squeeze(AScarbon(idxProd).carbonDataArrayInt(iLat,iLon,iTimeStep,1:iLastResolvedDepth)); % mg m-3

                    % Only operate on valid water columns (i.e., those without
                    % a single NaN pixel)
                    if (sum(isnan(waterColumnCarbon)) == 0)
                        depthResolution = AScarbon(idxProd).depthResolutionCarbon(1:iLastResolvedDepth);
                        pixelConcentrationPerUnitArea(iLat,iLon,iTimeStep) =... 
                            sum(waterColumnCarbon.*depthResolution); % mg m-3 --> mg m-2
                        pixelWaterColumnVolume(iLat,iLon,iTimeStep) =... 
                            sum(depthResolution).*thisPixelArea; % m3
                    end
                    
                elseif (strcmp(prodName,'bicep_poc_4km') || strcmp(prodName,'bicep_cphyto_9km'))

                    pixelConcentrationPerUnitArea(iLat,iLon,iTimeStep) =... 
                        AScarbon(idxProd).carbonDataArrayInt(iLat,iLon,iTimeStep)...
                        .*AScarbon(idxProd).depthIntegrationDataArrayInt(iLat,iLon,iTimeStep); % mg m-3 --> mg m-2
                    pixelWaterColumnVolume(iLat,iLon,iTimeStep) =... 
                        AScarbon(idxProd).depthIntegrationDataArrayInt(iLat,iLon,iTimeStep).*thisPixelArea; % m3

                elseif (strcmp(prodName,'bicep_npp_9km'))    
                    
                    pixelConcentrationPerUnitArea(iLat,iLon,iTimeStep) =... 
                        AScarbon(idxProd).carbonDataArrayInt(iLat,iLon,iTimeStep)...
                        .*INTEGRATION_TIME_NPP(iMonth); % mg m-2 d-1 --> mg m-2 month-1
                    pixelWaterColumnVolume(iLat,iLon,iTimeStep) =... 
                        AScarbon(idxProd).depthIntegrationDataArrayInt(iLat,iLon,iTimeStep).*thisPixelArea; % m3
                    
                end
            
                pixelStock(iLat,iLon,iTimeStep) =...
                    pixelConcentrationPerUnitArea(iLat,iLon,iTimeStep).*thisPixelArea; % mg m-2 (or mg m-2 month-1) --> mg (or mg month-1)
            
            end % iLat
        end % iLon
    
        % This water volume only uses 'valid' (non-NaN pixels). The water
        % volume could be bigger if all pixels were valid. Thus, this is not
        % the actual scene water volume, but the potential
        if ~all(isnan(reshape(pixelWaterColumnVolume(:,:,iTimeStep),[],1))) % only sum when there is at least one value that is not NaN (otherwise, leave the sum as NaN)
            sceneWaterVolume(iTimeStep) =... 
                sum(squeeze(pixelWaterColumnVolume(:,:,iTimeStep)),'all','omitnan'); % m3
        end

        % As recommended in Allison et al. (2010). Consider the fact that
        % there will be pixels without values so it would not be precise
        % calculating scene stocks without those
        sceneConcentrationPerUnitArea_lit(iTimeStep) =... 
            mean(pixelConcentrationPerUnitArea(:,:,iTimeStep),'all','omitnan'); % mg m-2
        sceneStock_lit(iTimeStep) =...
            sceneConcentrationPerUnitArea_lit(iTimeStep).*AScarbon(idxProd).sceneArea; % mg

        % My approach (works well when there are no missing pixels)
        if ~all(isnan(reshape(pixelStock(:,:,iTimeStep),[],1)))
            sceneStock_mine(iTimeStep) =...
                sum(squeeze(pixelStock(:,:,iTimeStep)),'all','omitnan'); % mg
        end

        % Update month timer
        iMonth = iMonth + 1;
        % Reset month timer
        if (iMonth > 12)
            iMonth = 1;
        end
    
    end % iTimeStep
    
    % My approach
    sceneConcentrationPerUnitArea_mine(:) = sceneStock_mine(:)./AScarbon(idxProd).sceneArea; % mg m-2
    
    % Save into output array
    AScarbon(idxProd).pixelConcentrationPerUnitArea = pixelConcentrationPerUnitArea;
    AScarbon(idxProd).pixelStock = pixelStock;                      
    AScarbon(idxProd).pixelWaterColumnVolume = pixelWaterColumnVolume;            
    AScarbon(idxProd).sceneConcentrationPerUnitArea_lit = sceneConcentrationPerUnitArea_lit;
    AScarbon(idxProd).sceneConcentrationPerUnitArea_mine = sceneConcentrationPerUnitArea_mine;
    AScarbon(idxProd).sceneStock_lit = sceneStock_lit;                   
    AScarbon(idxProd).sceneStock_mine = sceneStock_mine;                   
    AScarbon(idxProd).sceneWaterVolume = sceneWaterVolume;

end % iDataset

% % Check
% daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
% obsNpp_mgCm2d = ASbicep(2).dataset(:,:,13:24); % mg C m-2 d-1
% obsNpp_gCm2d = obsNpp_mgCm2d.*1e-3; % g C m-2 d-1
% 
% 
% obsNpp_mgCm2_per_month = NaN(size(obsNpp_mgCm2d));
% for iMonth = 1:12
%     obsNpp_mgCm2_per_month(:,:,iMonth) = obsNpp_mgCm2d(:,:,iMonth).*daysInMonth(iMonth); % mg C m-2 month-1
% end
% obsNpp_mgCm2_per_year = sum(obsNpp_mgCm2_per_month,3,'omitnan'); % mg C m-2 year-1
% obsNpp_mgC_per_year = obsNpp_mgCm2_per_year.*areaPixel; % mg C year-1
% obsNpp_MtC_year = sum(obsNpp_mgC_per_year,'all','omitnan')*1e-3*1e-6*1e-6; % Mt C year-1


%% Checks

if isCheckIntegrationCalc
    % Plot scene water volume where carbon stocks have been integrated
%     plotSceneWaterVolumeCarbonStock(AScarbon,carbonDatasetNames)
    % Plot comparison of stock integration approaches
%     plotComparisonDepthIntegrationApproaches(AScarbon,carbonDatasetNames,MILIGRAM_TO_GIGAGRAM,MILIGRAM_TO_GRAM)
    % Plot scene carbon stocks
    testDate = '2019-04-15';
    plotSceneCarbonStockProperties(AScarbon,carbonDatasetNames,testDate,MILIGRAM_TO_GIGAGRAM,MILIGRAM_TO_GRAM,areaStudy,enduranceSite)
    % Plot time-series properties
%     plotTimeSeriesCarbonStockProperties(AScarbon,carbonDatasetNames,MILIGRAM_TO_GIGAGRAM,MILIGRAM_TO_GRAM)
end

%% Group data by month to calculate monthly quantities

for iDataset = 1:numel(carbonDatasetNames)
    disp(iDataset)

    prodName = carbonDatasetNames{iDataset};
    idxProd = find(strcmp({AScarbon.ID}, prodName));

    latVector = AScarbon(idxProd).latVectorAreaStudy;
    lonVector = AScarbon(idxProd).lonVectorAreaStudy;
    time = AScarbon(idxProd).timeVectorCarbon;
    nLatPixels = numel(AScarbon(idxProd).latVectorAreaStudy);
    nLonPixels = numel(AScarbon(idxProd).lonVectorAreaStudy);

    pixelStock = MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).pixelStock; 
    sceneStock = MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).sceneStock_lit; 
    sceneConcentrationPerUnitArea = MILIGRAM_TO_GRAM.*AScarbon(idxProd).sceneConcentrationPerUnitArea_lit;

    % Use a table to arrange the scene data
    TS = table(sceneStock,sceneConcentrationPerUnitArea,time,...
        'VariableNames',{'stock','concperarea','date'});
    TS.month = month(TS.date);
    TS.year = year(TS.date);
    yearVectorData = min(TS.year):1:max(TS.year);
    
    % Calculate monthly mean variables
    sceneStockDistribByMonthAndYear_monthlyMean                    = NaN(12,numel(yearVectorData)); % Gg C
    sceneConcentrationPerUnitAreaDistribByMonthAndYear_monthlyMean = NaN(12,numel(yearVectorData)); % g C m-2
    pixelStock_monthlyMean                                         = NaN(nLatPixels,nLonPixels,12,numel(yearVectorData)); % Gg C

    for iYear = 1:numel(yearVectorData)
        for iMonth = 1:12

            % For tabulated variables (scene), extract the corresponding monthly values
            thisMonthAndYearStock =...
                TS.stock(TS.month(:) == iMonth & TS.year == yearVectorData(iYear));
            thisMonthAndYearConcentrationPerUnitArea =...
                TS.concperarea(TS.month(:) == iMonth & TS.year == yearVectorData(iYear));
            sceneStockDistribByMonthAndYear_monthlyMean(iMonth,iYear) =... 
                mean(thisMonthAndYearStock,'omitnan'); 
            sceneConcentrationPerUnitAreaDistribByMonthAndYear_monthlyMean(iMonth,iYear) =... 
                mean(thisMonthAndYearConcentrationPerUnitArea,'omitnan');

            % For non-tabulated variables (pixels in scene), find indices 
            % corresponding to the current month and year
            monthStartDate = datetime(yearVectorData(iYear),iMonth,1);
            monthEndDate = dateshift(monthStartDate,'end','month');
            thisMonthAndYearIdxs = time >= monthStartDate & time <= monthEndDate;
            % Extract data for the current month
            thisMonthAndYearPixelData = squeeze(pixelStock(:,:,thisMonthAndYearIdxs)); 
            % Calculate monthly mean
            for iLon = 1:nLonPixels
                for iLat = 1:nLatPixels
                    pixelStock_monthlyMean(iLat,iLon,iMonth,iYear) =... 
                        mean(squeeze(thisMonthAndYearPixelData(iLat,iLon,:)),'omitnan');
                end
            end

        end
    end
    
    % Save into output array
    AScarbon(idxProd).sceneStockDistribByMonthAndYear_monthlyMean = sceneStockDistribByMonthAndYear_monthlyMean;
    AScarbon(idxProd).sceneConcentrationPerUnitAreaDistribByMonthAndYear_monthlyMean = sceneConcentrationPerUnitAreaDistribByMonthAndYear_monthlyMean;                      
    AScarbon(idxProd).pixelStock_monthlyMean = pixelStock_monthlyMean;            

    % Plot histograms 
    if (isPlotHistogramsAverageStocks)
        plotHistogramsOfStocks(sceneStockDistribByMonthAndYear_monthlyMean,...
            sceneConcentrationPerUnitAreaDistribByMonthAndYear_monthlyMean,...
            yearVectorData,prodName)
    end
    
end % iDataset

%% Videos

for iDataset = 1:numel(carbonDatasetNames)
    disp(iDataset)

    prodName = carbonDatasetNames{iDataset};
    idxProd = find(strcmp({AScarbon.ID}, prodName));

        latVector = AScarbon(idxProd).latVectorAreaStudy;
    lonVector = AScarbon(idxProd).lonVectorAreaStudy;
    time = AScarbon(idxProd).timeVectorCarbon;
    nLatPixels = numel(AScarbon(idxProd).latVectorAreaStudy);
    nLonPixels = numel(AScarbon(idxProd).lonVectorAreaStudy);

    pixelStock = MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).pixelStock; 
    sceneStock = MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).sceneStock_lit; 
    sceneConcentrationPerUnitArea = MILIGRAM_TO_GRAM.*AScarbon(idxProd).sceneConcentrationPerUnitArea_lit;

    % Use a table to arrange the scene data
    TS = table(sceneStock,sceneConcentrationPerUnitArea,time,...
        'VariableNames',{'stock','concperarea','date'});
    TS.month = month(TS.date);
    TS.year = year(TS.date);
    yearVectorData = min(TS.year):1:max(TS.year);
    
    % We will create a video
    videoFileName = strcat('video_scene_stocks_',prodName,'.mp4');
    videoObject = VideoWriter(fullfile('.','data','processed',videoFileName), 'MPEG-4');
    videoObject.FrameRate = 2; % frames per second
    open(videoObject); % open the video file for writing

    if (isPlotScenesAverageStocks)
        plotScenesOfStocks(AScarbon(idxProd).pixelStock_monthlyMean,...
            latVector,lonVector,yearVectorData,areaStudy,enduranceSite,prodName,videoObject)
    end
    
end % iDataset
    
%% Save

save(fullfile('.','data','processed',filenameCarbonStocksAreaStudy),'AScarbon','-v7.3')
fprintf("\n...finished reading carbon products.\n")

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function plotScenePercentageMissingDataByYear(AScarbon,carbonDatasetNames,...
    areaStudy,enduranceSite)

nRows = 5;
nCols = 5;
maxSubplots = nRows*nCols;
    
for iDataset = 1:numel(carbonDatasetNames)
    
    prodName = carbonDatasetNames{iDataset};
    
    switch prodName
        case 'cmems_cphyto'
            depthLevel = 1;
            titleStr = 'CMEMS Cphyto';
        case 'bicep_poc_4km'
            titleStr = 'BICEP POC';
        case 'bicep_npp_9km'
            titleStr = 'BICEP NPP';
        case 'bicep_cphyto_9km'
            titleStr = 'BICEP Cphyto';
    end
        
    idxProd = find(strcmp({AScarbon.ID}, prodName));
    
    if (strcmp(carbonDatasetNames{iDataset},'cmems_cphyto'))
        dataArray = squeeze(AScarbon(idxProd).carbonDataArrayInt(:,:,:,depthLevel));
        cbStr = '% missing data at 1 m depth';
    else
        dataArray = AScarbon(idxProd).carbonDataArrayInt;
        cbStr = '% missing data';
    end
    timeVector = AScarbon(idxProd).timeVectorCarbon;
    latVector = AScarbon(idxProd).latVectorAreaStudy;
    lonVector = AScarbon(idxProd).lonVectorAreaStudy;

    yearVector = (min(year(timeVector)):1:max(year(timeVector)))';
    yearFracMisingData = zeros(numel(latVector),numel(lonVector),numel(yearVector));
    for iYear = 1:numel(yearVector)
        thisYearIdxs = find(year(timeVector) == yearVector(iYear));
        nImagesInThisYear = numel(thisYearIdxs);
        for iCol = 1:numel(lonVector)
            for iRow = 1:numel(latVector)
                yearFracMisingData(iRow,iCol,iYear) =...
                    100*(sum(isnan(squeeze(dataArray(iRow,iCol,thisYearIdxs))))/nImagesInThisYear);
            end
        end
    end
    
    figure()               
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.75],'Color','w') 
    haxis = zeros(nRows,nCols);
        
    for iYear = 1:numel(yearVector)
        
        % If the iteration exceeds 25, stop plotting
        if iYear > maxSubplots
            break;
        end

        haxis(iYear) = subaxis(nRows,nCols,iYear,'Spacing',0.015,'Padding',0.010,'Margin',0.04);
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
        m_pcolor(lonVector,latVector,yearFracMisingData(:,:,iYear)) 
        shading flat
        colormap(brewermap(10,'*Blues')) % same as maximum of colourbar axis
        caxis([0 10])
        set(gca,'color','w'); % or whatever color you want for the background

        % Show axis labels for border subplots
        if (iYear == 1 || iYear == 6 || iYear == 11 || iYear == 16)
            m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[]);
        elseif (iYear == 22 || iYear == 23 || iYear == 24 || iYear == 25)
            m_grid('linewi',1,'tickdir','out','FontSize',8,'yticklabel',[]);
        elseif (iYear == 21)
            m_grid('linewi',1,'tickdir','out','FontSize',8)
        else
            m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
        end

        m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
        title(num2str(yearVector(iYear)),'FontSize',11)
        box on
        hold on

        if (iYear == 5)
            cb = colorbar('Location','eastoutside');
            cb.Position(1) = cb.Position(1) + 0.06;
            cb.Position(2) = cb.Position(2) - 0.67;
            cb.Position(3) = 0.020; % WIDTH
            cb.Position(4) = 0.70; % LENGTH
            cb.Label.String = cbStr; 
            cb.FontSize = 10;
        end

    end % iYear
    
    % Title over group of subplots
    a = axes;
    t = title(titleStr,'FontSize',12);
    a.Visible = 'off';
    t.Visible = 'on';
    a.Position(1) = 0.50 - a.Position(3)/2;
    a.Position(2) = a.Position(2) + 0.05; 

    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_check_carbon_scene_percentage_missing_data_by_year_',prodName,'.png')),'Resolution',600)
   
end % iDataset

end % plotScenePercentageMissingDataByYear

% *************************************************************************

function plotScenePercentageMissingDataByMonth(AScarbon,carbonDatasetNames,...
    areaStudy,enduranceSite)

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for iDataset = 1:numel(carbonDatasetNames)
    
    prodName = carbonDatasetNames{iDataset};
    
    switch prodName
        case 'cmems_cphyto'
            depthLevel = 1;
            titleStr = 'CMEMS Cphyto';
        case 'bicep_poc_4km'
            titleStr = 'BICEP POC';
        case 'bicep_npp_9km'
            titleStr = 'BICEP NPP';
        case 'bicep_cphyto_9km'
            titleStr = 'BICEP Cphyto';
    end
        
    idxProd = find(strcmp({AScarbon.ID}, prodName));
    
    if (strcmp(carbonDatasetNames{iDataset},'cmems_cphyto'))
        dataArray = squeeze(AScarbon(idxProd).carbonDataArrayInt(:,:,:,depthLevel));
        cbStr = '% missing data at 1 m depth in the period 1998-2022';
    else
        dataArray = AScarbon(idxProd).carbonDataArrayInt;
        cbStr = '% missing data in the period 1998-2020';
    end
    timeVector = AScarbon(idxProd).timeVectorCarbon;
    latVector = AScarbon(idxProd).latVectorAreaStudy;
    lonVector = AScarbon(idxProd).lonVectorAreaStudy;

    monthFracMisingData = zeros(numel(latVector),numel(lonVector),12);
    for iMonth = 1:12
        thisMonthIdxs = find(month(timeVector) == iMonth);
        nImagesInMonthForPeriod = numel(thisMonthIdxs);
        for iCol = 1:numel(lonVector)
            for iRow = 1:numel(latVector)
                monthFracMisingData(iRow,iCol,iMonth) =...
                    100*(sum(isnan(squeeze(dataArray(iRow,iCol,thisMonthIdxs))))/nImagesInMonthForPeriod);
            end
        end
    end

    figure()               
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.35 0.48],'Color','w') 
    haxis = zeros(3,4);

    for iMonth = 1:12

        haxis(iMonth) = subaxis(3,4,iMonth,'Spacing',0.015,'Padding',0.010,'Margin',0.04);
        ax(iMonth).pos = get(haxis(iMonth),'Position');

        % Move subplots
        if (iMonth == 1 || iMonth == 5 || iMonth == 9)
            ax(iMonth).pos(1) = ax(iMonth).pos(1) + 0.025;
        elseif (iMonth == 2 || iMonth == 6 || iMonth == 10)
            ax(iMonth).pos(1) = ax(iMonth).pos(1) + 0;
        elseif (iMonth == 3 || iMonth == 7 || iMonth == 11)
            ax(iMonth).pos(1) = ax(iMonth).pos(1) - 0.025;
        elseif (iMonth == 4 || iMonth == 8 || iMonth == 12)
            ax(iMonth).pos(1) = ax(iMonth).pos(1) - 0.050;
        end 
        ax(iMonth).pos(2) = ax(iMonth).pos(2) - 0.01;
        set(haxis(iMonth),'Position',ax(iMonth).pos) 

        m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
        m_pcolor(lonVector,latVector,monthFracMisingData(:,:,iMonth)) 
        shading flat
        colormap(brewermap(10,'*Blues')) % same as maximum of colourbar axis
        caxis([0 10])
        set(gca,'color','w'); % or whatever color you want for the background

        % Show axis labels for border subplots
        if (iMonth == 1 || iMonth == 5)
            m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[]);
        elseif (iMonth == 10 || iMonth == 11 || iMonth == 12)
            m_grid('linewi',1,'tickdir','out','FontSize',8,'yticklabel',[]);
        elseif (iMonth == 9)
            m_grid('linewi',1,'tickdir','out','FontSize',8)
        else
            m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
        end

        m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
        title(months{iMonth},'FontSize',11)
        box on
        hold on

    end % iYear

    % Colour bar settings
    cb = colorbar('Location','eastoutside');
    cb.Position(1) = cb.Position(1) + 0.075;
    cb.Position(2) = cb.Position(2) + 0.045;
    cb.Position(3) = 0.020; % WIDTH
    cb.Position(4) = 0.70; % LENGTH
    cb.Label.String = cbStr; 
    cb.FontSize = 10;

    % Title over group of subplots
    a = axes;
    t = title(titleStr,'FontSize',12);
    a.Visible = 'off';
    t.Visible = 'on';
    a.Position(1) = 0.50 - a.Position(3)/2;
    a.Position(2) = a.Position(2) + 0.03; 
    
    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_check_carbon_scene_percentage_missing_data_by_month_',prodName,'.png')),'Resolution',600)
 
end % iDataset

end % plotScenePercentageMissingDataByMonth

% *************************************************************************

function plotComparisonBeforeAndAfterRegridding(AScarbon,carbonDatasetNames,...
    testDate,areaStudy,enduranceSite)
  
cmap = brewermap(1000,'*RdYlBu');

for iDataset = 1:numel(carbonDatasetNames)
    
    prodName = carbonDatasetNames{iDataset};

    switch prodName
        case 'cmems_cphyto'
            depthLevel = 1;
            cbStr = 'CMEMS Cphyto stock (mg m^{-3}) at 1st m depth';
            cbMin = 0;
            cbMax = 100;
        case 'bicep_poc_4km'
            cbStr = 'BICEP POC stock (mg m^{-3}) water-column integrated';
            cbMin = 0;
            cbMax = 400;
        case 'bicep_npp_9km'
            cbStr = 'BICEP NPP (mg m^{-2} d^{-1}) water-column integrated';
            cbMin = 100;
            cbMax = 1000;
        case 'bicep_cphyto_9km'
            cbStr = 'BICEP Cphyto stock (mg m^{-3}) water-column integrated';
            cbMin = 10;
            cbMax = 40;
    end
    
    idxProd = find(strcmp({AScarbon.ID}, prodName));
    
    time = AScarbon(idxProd).timeVectorCarbon;
    [~, idxClosestTime] = min(abs(time - testDate));
    titleStr1 = datestr(time(idxClosestTime), 'yyyy-mm-dd');
 
    if (strcmp(carbonDatasetNames{iDataset},'cmems_cphyto'))
        arrayOriginal = squeeze(AScarbon(idxProd).carbonDataArray(:,:,:,depthLevel));
        arrayInterpolated = squeeze(AScarbon(idxProd).carbonDataArrayInt(:,:,:,depthLevel));
    else
        arrayOriginal = AScarbon(idxProd).carbonDataArray;
        arrayInterpolated = AScarbon(idxProd).carbonDataArrayInt;
    end

    latOriginal = AScarbon(idxProd).latVectorCarbon;
    lonOriginal = AScarbon(idxProd).lonVectorCarbon;
    latRegridded = AScarbon(idxProd).latVectorAreaStudy;
    lonRegridded = AScarbon(idxProd).lonVectorAreaStudy;
    
    figure()               
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.35],'Color','w') 
    haxis = zeros(1,2);

    for iSubplot = 1:2

        if (iSubplot == 1)
            data = squeeze(arrayOriginal(:,:,idxClosestTime));
            lat = latOriginal;
            lon = lonOriginal;
            titleStr2 = 'Original data product';
        else
            data = squeeze(arrayInterpolated(:,:,idxClosestTime));
            lat = latRegridded;
            lon = lonRegridded;
            titleStr2 = 'Regridded data product';
        end

        haxis(iSubplot) = subaxis(1,2,iSubplot,'Spacing',0.020,'Padding',0.020,'Margin',0.08);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.035;
        set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

        m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
        m_pcolor(lon,lat,data)
        shading flat
        colormap(cmap)
        caxis([cbMin cbMax])
        set(gca, 'color', 'w'); % or whatever color you want for the background

        m_grid('linewi',1,'tickdir','out','FontSize',8);
        m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
        title([titleStr1, newline, titleStr2],'FontSize',11)
        box on
        hold on
        
    end % iSubplot

    % Colour bar settings
    cb = colorbar('Location','southoutside');
    cb.Position(1) = cb.Position(1) - 0.35;
    cb.Position(2) = cb.Position(2) - 0.07;
    cb.Position(3) = 0.60; 
    cb.Position(4) = 0.015;
    cb.Label.String = cbStr; 
    cb.FontSize = 11;

    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_check_carbon_regridding_before_and_after_',titleStr1,'_',prodName,'.png')),'Resolution',600)
    
end % iDataset

end % plotComparisonBeforeAndAfterRegridding

% *************************************************************************

function plotComparisonDepthIntegrationProducts(AScmems,AScarbon)

% Define legend strings and other plotting parameters
nDepthIntegrationProducts = 5;
lgStr = {'MLD from CMEMS physical reanalysis (daily)',...
         'MLD emergent from CMEMS-Cphyto integration (daily)',...
         'MLD from MIMOC, used in BICEP-POC (monthly climatology)',...
         'MLD that comes with BICEP-Cphyto (monthly climatology)',...
         'Zeu calculated from CMEMS biogechemical reanalysis (daily)'};
colourProduct = flipud(brewermap(nDepthIntegrationProducts,'*Paired'));

startDate = datetime('1998-01-01');
endDate = datetime('2020-12-31');

% Retrieve data arrays and time vectors for each product

% MLD product from CMEMS (regrid and interpolate)
idxMldCmems = find(strcmp({AScmems.ID}, 'mod_phy_reg_mld'));
cmemsMld = AScmems(idxMldCmems).dataset;
cmemsMldTime = AScmems(idxMldCmems).time;

% Zeu calculated from CMEMS's Kd (regrid and interpolate)
idxKdCmems = find(strcmp({AScmems.ID}, 'mod_bgc_reg_kd'));
[cmemsZeu] = calculateZeuFromCmemsKd(AScmems);
cmemsZeuTime = AScmems(idxKdCmems).time;

% Emergent MLD from CMEMS's Cphyto product (counting depth layers
% occupied by values)
idxCphytoCmems = find(strcmp({AScarbon.ID}, 'cmems_cphyto'));
cmemsCphytoEmergentMldTime = AScarbon(idxCphytoCmems).timeVectorCarbon;
cmemsCphytoDepths = AScarbon(idxCphytoCmems).depthVectorCarbon;
cmemsCphytoNoResolvedDepths = AScarbon(idxCphytoCmems).nResolvedDepthLevels;
cmemsCphytoEmergentMld = NaN(size(cmemsCphytoNoResolvedDepths));
nonZeroDepthIndices = cmemsCphytoNoResolvedDepths > 0;
cmemsCphytoEmergentMld(nonZeroDepthIndices) =... 
    cmemsCphytoDepths(cmemsCphytoNoResolvedDepths(nonZeroDepthIndices));
   
% MLD used in BICEP's POC (MIMOC)
idxBicepPoc = find(strcmp({AScarbon.ID}, 'bicep_poc_4km'));
mimocMldTime = AScarbon(idxBicepPoc).timeVectorCarbon;
mimocMld = AScarbon(idxBicepPoc).depthIntegrationDataArrayInt;

% MLD used in BICEP's Cphyto product
idxCphytoBicep = find(strcmp({AScarbon.ID}, 'bicep_cphyto_9km'));
bicepCphytoMldTime = AScarbon(idxCphytoBicep).timeVectorCarbon;
bicepCphytoMld = AScarbon(idxCphytoBicep).depthIntegrationDataArrayInt;

% Plot

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.57],'Color','w') 
haxis = zeros(2,1);

for iSubplot = 1:2
    
    haxis(iSubplot) = subaxis(2,1,iSubplot,'Spacing',0.020,'Padding',0.020,'Margin',0.08);
    ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.04;
    ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.015;
    set(haxis(iSubplot),'Position',ax(iSubplot).pos) 
    
    for iProduct = 1:nDepthIntegrationProducts

        if (iProduct == 1)
            depthArray = cmemsMld;
            timeVector = cmemsMldTime;
        elseif (iProduct == 2)
            depthArray = cmemsCphytoEmergentMld;
            timeVector = cmemsCphytoEmergentMldTime;
        elseif (iProduct == 3) 
            depthArray = mimocMld;
            timeVector = mimocMldTime;
        elseif (iProduct == 4)
            depthArray = bicepCphytoMld;
            timeVector = bicepCphytoMldTime;
        elseif (iProduct == 5)
            depthArray = cmemsZeu;
            timeVector = cmemsZeuTime;
        end

        depthArraySceneAveraged = NaN(1,length(timeVector)); 
        for iTimeStep = 1:length(timeVector)
            thisScene = reshape(depthArray(:,:,iTimeStep),[],1);
            depthArraySceneAveraged(iTimeStep) = mean(thisScene,'omitnan');
        end

        % For daily products, we will use scatter
        if (iProduct == 1 || iProduct == 2 || iProduct == 5)
            scatter(haxis(iSubplot),timeVector,depthArraySceneAveraged,...
                6,colourProduct(iProduct,:),'filled'); hold on; 
        % For climatologies, we will use lines
        else
            plot(haxis(iSubplot),timeVector,depthArraySceneAveraged,...
                'Color',colourProduct(iProduct,:),'LineWidth',1); hold on; % longer time steps, a line provides better visualisation than scatter dots
        end
        if (iSubplot == 1)
            xlim([startDate endDate])
        elseif (iSubplot == 2) % zoom in
            xlim([datetime('1998-01-01') datetime('2003-12-31')])
            xlabel('Time')
        end

    end
    hold off

    set(gca,'Ydir','reverse')
    ylim([0 70])
    ylabel('Depth (m)')
    grid on
    box on

    if (iSubplot == 1)
        lg = legend(lgStr,'NumColumns',1);
        lg.Position(1) = 0.08; lg.Position(2) = 0.88; 
        lg.Orientation = 'horizontal';
        lg.ItemTokenSize = [10,20];
        lg.FontSize = 10;
        set(lg,'Box','off','Color','w')
    end
    
end

exportgraphics(gcf,fullfile('.','figures','depth_carbon_products.png'),'Resolution',600)
    
end % plotComparisonDepthIntegrationProducts

% *************************************************************************

function plotAreaOfPixel(areaPixel,latVectorAreaStudy,lonVectorAreaStudy,...
    areaStudy,enduranceSite,titleStr,prodName)

figure()               
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.20 0.35],'Color','w') 
m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
m_pcolor(lonVectorAreaStudy,latVectorAreaStudy,areaPixel)
shading flat
cmap = brewermap(100,'*RdYlBu');
colormap(cmap)
caxis([min(areaPixel,[],'all') max(areaPixel,[],'all')])
set(gca,'color','w');
m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
box on
title(titleStr)
cb = colorbar('Location','southoutside');
cb.Label.String = 'Area of a pixel (km^{2})'; 
cb.FontSize = 11;

exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_check_carbon_area_pixel_',prodName,'.png')),'Resolution',600)
    
end % plotAreaOfPixel

% *************************************************************************

function plotSceneWaterVolumeCarbonStock(AScarbon,carbonDatasetNames)

startDate = datetime('1998-01-01');
endDate = datetime('2020-12-31');

nCarbonDatasets = numel(carbonDatasetNames);
colourDataset = brewermap(nCarbonDatasets,'*Spectral');

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.27],'Color','w') 
axh = axes('Position', [0.10 0.15 0.70 0.75]);

for iDataset = 1:nCarbonDatasets
    prodName = carbonDatasetNames{iDataset};
    idxProd = find(strcmp({AScarbon.ID}, prodName));
    plot(axh,AScarbon(idxProd).timeVectorCarbon,AScarbon(idxProd).sceneWaterVolume,...
        'Color',colourDataset(iDataset,:),'LineWidth',1); hold on;
    xlim([startDate endDate])
    if (iDataset == nCarbonDatasets)
        lg = legend('CMEMS Cphyto','BICEP POC','BICEP NPP','BICEP Cphyto','NumColumns',1);
        lg.Position(1) = 0.80; lg.Position(2) = 0.68; 
        lg.Orientation = 'vertical';
        lg.ItemTokenSize = [10,20];
        lg.FontSize = 10;
        set(lg,'Box','on','Color','w')
    end
end
hold off    
xlabel('Time')
ylabel('Scene water volume (m^{3})')
grid on
box on

exportgraphics(gcf,fullfile('.','figures','oc_check_carbon_scene_volume.png'),'Resolution',600)
    
end % plotSceneWaterVolumeCarbonStock

% *************************************************************************

function plotComparisonDepthIntegrationApproaches(AScarbon,carbonDatasetNames,...
    MILIGRAM_TO_GIGAGRAM,MILIGRAM_TO_GRAM)

startDate = datetime('1998-01-01');
endDate = datetime('2020-12-31');

nDatasetsLoaded = numel(carbonDatasetNames);
colourDataset = brewermap(nDatasetsLoaded*2,'*Paired');

lgStr = {'CMEMS Cphyto - A2010',...
         'CMEMS Cphyto - mine',...
         'BICEP POC - A2010',...
         'BICEP POC - mine',...
         'BICEP Cphyto - A2010',...
         'BICEP Cphyto - mine'};

lgStrNpp = {'BICEP NPP - A2010',...
            'BICEP NPP - mine'};
     
for iFigure = 1:2
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.47],'Color','w')
    haxis = zeros(2,1);

    % 1st sublot = all years, 2nd subplot = zooms in the first few years for betetr visualisation
    for iSubplot = 1:2 
        
        haxis(iSubplot) = subaxis(2,1,iSubplot,'Spacing',0.020,'Padding',0.020,'Margin',0.08);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.015;
        set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

        iScatter = 1;

        for iDataset = 1:nDatasetsLoaded % we won't plot the NPP one as it's been integrated over 30 days and, thus, the values will be very high
            
            prodName = carbonDatasetNames{iDataset};
            if (strcmp(prodName,'bicep_npp_9km'))
                continue % skip plotting 'bicep_npp_9km'
            end
            idxProd = find(strcmp({AScarbon.ID}, prodName));

            if (iFigure == 1)
                scatter(haxis(iSubplot),AScarbon(idxProd).timeVectorCarbon,...
                    MILIGRAM_TO_GRAM.*AScarbon(idxProd).sceneConcentrationPerUnitArea_lit,...
                    10,colourDataset(iScatter,:),'filled'); hold on;
                iScatter = iScatter + 1;
                scatter(haxis(iSubplot),AScarbon(idxProd).timeVectorCarbon,...
                    MILIGRAM_TO_GRAM.*AScarbon(idxProd).sceneConcentrationPerUnitArea_mine,...
                    4,colourDataset(iScatter,:),'filled'); hold on;
                iScatter = iScatter + 1;
                ylabelStr = sprintf('Scene average concentration \nper unit area (g C m^{-2})');
            elseif (iFigure == 2)
                scatter(haxis(iSubplot),AScarbon(idxProd).timeVectorCarbon,...
                    MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).sceneStock_lit,...
                    10,colourDataset(iScatter,:),'filled'); hold on; 
                iScatter = iScatter + 1;
                scatter(haxis(iSubplot),AScarbon(idxProd).timeVectorCarbon,...
                    MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).sceneStock_mine,...
                    4,colourDataset(iScatter,:),'filled'); hold on;
                iScatter = iScatter + 1;
                ylabelStr = 'Scene stock (Gg C)';
            end

        end % iDataset
        
        if (iSubplot == 1)
            xlim([startDate endDate])
        elseif (iSubplot == 2)
            xlim([datetime('1998-01-01') datetime('2001-12-31')])
            xlabel('Time')
        end

        ylabel(ylabelStr)
        grid on
        box on

        if (iSubplot == 1 && iFigure < 3)
            lg = legend(lgStr,'NumColumns',2);
            lg.Position(1) = 0.25;
            lg.Position(2) = 0.90; 
            lg.Orientation = 'horizontal';
            lg.ItemTokenSize = [10,20];
            set(lg,'Box','on','Color','w')
        end
        
    end % iSubplot
    
    if (iFigure == 1)
        suffixStr = 'concarea';
    elseif (iFigure == 2)
        suffixStr = 'stock';
    end

    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_check_carbon_depth_integration_approaches_',suffixStr,'.png')),'Resolution',600)
 
end % iFigure

% Plot NPP separately, y axis very different

for iFigure = 1:2
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.47],'Color','w')
    haxis = zeros(2,1);

    % 1st sublot = all years, 2nd subplot = zooms in the first few years for betetr visualisation
    for iSubplot = 1:2 

        haxis(iSubplot) = subaxis(2,1,iSubplot,'Spacing',0.020,'Padding',0.020,'Margin',0.08);
        ax(iSubplot).pos = get(haxis(iSubplot),'Position');
        ax(iSubplot).pos(2) = ax(iSubplot).pos(2) - 0.015;
        set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

        idxProd = find(strcmp({AScarbon.ID}, 'bicep_npp_9km'));

        if (iFigure == 1)
            iScatterNpp1 = iScatter;
            scatter(haxis(iSubplot),AScarbon(idxProd).timeVectorCarbon,...
                MILIGRAM_TO_GRAM.*AScarbon(idxProd).sceneConcentrationPerUnitArea_lit,...
                10,colourDataset(iScatterNpp1,:),'filled'); hold on; 
            iScatterNpp2 = iScatterNpp1 + 1;
            scatter(haxis(iSubplot),AScarbon(idxProd).timeVectorCarbon,...
                MILIGRAM_TO_GRAM.*AScarbon(idxProd).sceneConcentrationPerUnitArea_mine,...
                4,colourDataset(iScatterNpp2,:),'filled'); hold on;
            ylabelStr = sprintf('Monthly-integrated concentration \nper unit area (g C m^{-2})');
        elseif (iFigure == 2)
            scatter(haxis(iSubplot),AScarbon(idxProd).timeVectorCarbon,...
                MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).sceneStock_lit,...
                10,colourDataset(iScatterNpp1,:),'filled'); hold on; 
            scatter(haxis(iSubplot),AScarbon(idxProd).timeVectorCarbon,...
                MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).sceneStock_mine,...
                4,colourDataset(iScatterNpp2,:),'filled'); hold on;
            ylabelStr = sprintf('Monthly-integrated scene \nstock (Gg C)');   
        end

        if (iSubplot == 1)
            xlim([startDate endDate])
        elseif (iSubplot == 2)
            xlim([datetime('1998-01-01') datetime('2001-12-31')])
            xlabel('Time')
        end

        ylabel(ylabelStr)
        grid on
        box on

        if (iSubplot == 1)
            lg = legend(lgStrNpp,'NumColumns',2);
            lg.Position(1) = 0.25;
            lg.Position(2) = 0.90; 
            lg.Orientation = 'horizontal';
            lg.ItemTokenSize = [10,20];
            set(lg,'Box','on','Color','w')
        end

    end % iSubplot

    if (iFigure == 1)
        suffixStr = 'npp_concarea';
    elseif (iFigure == 2)
        suffixStr = 'npp_stock';
    end
 
    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_check_carbon_depth_integration_approaches_',suffixStr,'.png')),'Resolution',600)

end % iFigure

end % plotComparisonDepthIntegrationApproaches

% *************************************************************************

function plotSceneCarbonStockProperties(AScarbon,carbonDatasetNames,...
    testDate,MILIGRAM_TO_GIGAGRAM,MILIGRAM_TO_GRAM,areaStudy,enduranceSite)

cmap = brewermap(1000,'*RdYlBu');

for iDataset = 1:numel(carbonDatasetNames) 
    
    prodName = carbonDatasetNames{iDataset};
    
    switch prodName
        case 'cmems_cphyto'
            titleStr1 = 'CMEMS Cphyto';
            cbStr1 = 'Concentration per unit area (g C m-2)';
            cbStr2 = 'Depth-integrated stock (Gg C)';
            cbMin1 = 0.5;
            cbMax1 = 3.0;
            cbMin2 = 0.01;
            cbMax2 = 0.1;
        case 'bicep_poc_4km'
            titleStr1 = 'BICEP POC';
            cbStr1 = 'Concentration per unit area (g C m-2)';
            cbStr2 = 'Depth-integrated stock (Gg C)';
            cbMin1 = 3.0;
            cbMax1 = 6.0;
            cbMin2 = 0.01;
            cbMax2 = 0.1;
        case 'bicep_npp_9km'
            titleStr1 = 'BICEP NPP';
            cbStr1 = 'Monthly integrated (g C m-2 month-1)';
            cbStr2 = 'Monthly integrated (Gg C month-1)';
            cbMin1 = 10;
            cbMax1 = 30;
            cbMin2 = 0.10;
            cbMax2 = 1.0;
        case 'bicep_cphyto_9km'
            titleStr1 = 'BICEP Cphyto';
            cbStr1 = 'Concentration of per unit area (g C m-2)';
            cbStr2 = 'Depth-integrated stock (Gg C)';
            cbMin1 = 0.5;
            cbMax1 = 3.0;
            cbMin2 = 0.01;
            cbMax2 = 0.1;
    end
    idxProd = find(strcmp({AScarbon.ID}, prodName));

    if (strcmp(prodName,'cmems_cphyto') || strcmp(prodName,'bicep_poc_4km') || strcmp(prodName,'bicep_cphyto_9km'))
        
        disp(prodName)
        disp(idxProd)
        
        [~, idxClosestTime] = min(abs(AScarbon(idxProd).timeVectorCarbon - testDate));
        titleStr2 = datestr(AScarbon(idxProd).timeVectorCarbon(idxClosestTime), 'yyyy-mm-dd');

        figure()               
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.38 0.35],'Color','w') 
        haxis = zeros(2,1);
    
        for iSubplot = 1:2
            % Get data for pixel concentration per unit area
            if (iSubplot == 1)
                data = MILIGRAM_TO_GRAM.*squeeze(AScarbon(idxProd).pixelConcentrationPerUnitArea(:,:,idxClosestTime));
            % Get data for depth-integrated stock
            elseif (iSubplot == 2)
                data = MILIGRAM_TO_GIGAGRAM.*squeeze(AScarbon(idxProd).pixelStock(:,:,idxClosestTime)); 
            end

            haxis(iSubplot) = subaxis(1,2,iSubplot,'Spacing',0.015,'Padding',0.010,'Margin',0.08);
            ax(iSubplot).pos = get(haxis(iSubplot),'Position');
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.035;
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

            m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
            m_pcolor(AScarbon(idxProd).lonVectorAreaStudy,AScarbon(idxProd).latVectorAreaStudy,data)
            shading flat
            if (iSubplot == 1)
                caxis([cbMin1 cbMax1])
            else
                caxis([cbMin2 cbMax2])
            end
            colormap(cmap)
            set(gca,'color','w'); 
            m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
            m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
            title([titleStr1, newline, titleStr2],'FontSize',11)
            box on
            hold on

            cb = colorbar('Location','southoutside');
            cb.Position(1) = cb.Position(1);
            cb.Position(2) = cb.Position(2) - 0.05;
            cb.Position(3) = 0.35;
            cb.Position(4) = 0.015; 
            if (iSubplot == 1)
                cb.Label.String = cbStr1; 
            elseif (iSubplot == 2)
                cb.Label.String = cbStr2; 
            end
            cb.FontSize = 11;

        end % iSubplot
    
    elseif (strcmp(prodName,'bicep_npp_9km'))
        
        [~, idxClosestTime] = min(abs(AScarbon(idxProd).timeVectorCarbon - testDate));

        idxsMonthsTargetYear = find(year(AScarbon(idxProd).timeVectorCarbon) == year(datetime(testDate)));

        figure()               
        set(gcf,'Units','Normalized','Position',[0.01 0.05 0.55 0.35],'Color','w') 
        haxis = zeros(3,1);
        
        for iSubplot = 1:3

            if (iSubplot == 1)
                data = AScarbon(idxProd).carbonDataArrayInt(:,:,idxClosestTime); % mg C m-2 d-1
                
                % Convert the date vector to a serial date number
                serialDateNumber = datenum(testDate);
                monthName = datestr(serialDateNumber, 'mmmm');
                titleStr2 = sprintf('Daily NPP for an average day in %s %d',monthName,year(datetime(testDate)));
                
            % Get data for depth-integrated stock for the month
            elseif (iSubplot == 2)
%                 data = MILIGRAM_TO_GIGAGRAM.*squeeze(AScarbon(idxProd).pixelStock(:,:,idxClosestTime)); % Gg C month-1
                data = squeeze(AScarbon(idxProd).pixelConcentrationPerUnitArea(:,:,idxClosestTime)).*1e-3; % mg m-2 month-1 --> g m-2 month-1
                
                % Convert the date vector to a serial date number
                serialDateNumber = datenum(testDate);
                monthName = datestr(serialDateNumber, 'mmmm');
%                 titleStr2 = sprintf('%s stock per pixel',monthName);
                titleStr2 = sprintf('%s monthly productivity',monthName);
                
            % Get data for depth-integrated stock for the year
            elseif (iSubplot == 3)
%                 data = MILIGRAM_TO_GIGAGRAM.*sum(AScarbon(idxProd).pixelStock(:,:,idxsMonthsTargetYear),3); % Gg C year-1
                data = sum(AScarbon(idxProd).pixelConcentrationPerUnitArea(:,:,idxsMonthsTargetYear),3).*1e-3; % mg m-2 month-1 --> g m-2 year-1
                
%                 sceneStock = sum(data,'all','omitnan').*1e-3; % Gg C year
%                 --> Mt C year-1

                sceneStock = sum(data.*AScarbon(idxProd).areaPixel(:,:),'all','omitnan').*1e-12; % g m-2 year-1 --> Mt C year-1
                fprintf('\nAnnually integrated carbon fixation of %4.1f mega tonnes of carbon', sceneStock)
%                 titleStr2 = 'Annual stock per pixel';
                titleStr2 = 'Annual productivity';
            end

            haxis(iSubplot) = subaxis(1,3,iSubplot,'Spacing',0.015,'Padding',0.010,'Margin',0.08);
            ax(iSubplot).pos = get(haxis(iSubplot),'Position');
            ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.035;
            set(haxis(iSubplot),'Position',ax(iSubplot).pos) 

            m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
            m_pcolor(AScarbon(idxProd).lonVectorAreaStudy,AScarbon(idxProd).latVectorAreaStudy,data)
            shading flat
            if (iSubplot == 1)
                caxis([200 1000])
            elseif (iSubplot == 2)
%                 caxis([0.2 1.0])
                caxis([10 30])
            elseif (iSubplot == 3)
%                 caxis([2 10])
                caxis([100 200])
            end
            colormap(cmap)
            set(gca,'color','w'); 
            m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
            m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
            title([titleStr1, newline, titleStr2],'FontSize',11)
            box on
            hold on

            cb = colorbar('Location','southoutside');
%             cb.Position(1) = cb.Position(1);
%             cb.Position(2) = cb.Position(2) - 0.05;
%             cb.Position(3) = 0.35;
            cb.Position(4) = 0.015; 
            if (iSubplot == 1)
                cb.Label.String = 'mg C m^{-2} d ^{-1}'; 
            elseif (iSubplot == 2)
%                 cb.Label.String = '10^{-3} Gt C month^{-1}'; 
                cb.Label.String = 'g C m^{-2} month^{-1}';
            elseif (iSubplot == 3)
%                 cb.Label.String = '10^{-3} Gt C year^{-1}'; 
                cb.Label.String = 'g C m^{-2} year^{-1}'; 
            end
            cb.FontSize = 11;

        end % iSubplot

    end
    
    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_check_carbon_scene_stocks_',prodName,'_2','.png')),'Resolution',600)

end % iDataset

end % plotSceneCarbonStockProperties

% *************************************************************************

function plotTimeSeriesCarbonStockProperties(AScarbon,carbonDatasetNames,...
    MILIGRAM_TO_GIGAGRAM,MILIGRAM_TO_GRAM)

for iDataset = 1:numel(carbonDatasetNames)
    
    prodName = carbonDatasetNames{iDataset};

    switch prodName
        case 'cmems_cphyto'
            titleStr = 'CMEMS Cphyto';
            ylabelStr1 = 'Stock (Gg C)';
            ylabelStr2 = 'Concentration per unit area (g C m-2)';
            yMin1 = 0;
            yMax1 = 100;
            yMin2 = 0;
            yMax2 = 4;
        case 'bicep_poc_4km'
            titleStr = 'BICEP POC';
            ylabelStr1 = 'Stock (Gg C)';
            ylabelStr2 = 'Concentration per unit area (g C m-2)';
            yMin1 = 0;
            yMax1 = 300;
            yMin2 = 0;
            yMax2 = 15;
        case 'bicep_npp_9km'
            titleStr = 'BICEP NPP';
            ylabelStr1 = 'Monthly integrated stock (Gg C)';
            ylabelStr2 = sprintf('Monthly integrated concentration \nper unit area (g C m-2)');
            yMin1 = 0;
            yMax1 = 500;
            yMin2 = 0;
            yMax2 = 40;
        case 'bicep_cphyto_9km'
            titleStr = 'BICEP Cphyto';
            ylabelStr1 = 'Stock (Gg C)';
            ylabelStr2 = 'Concentration per unit area (g C m-2)';
            yMin1 = 0;
            yMax1 = 100;
            yMin2 = 0;
            yMax2 = 4;
    end
    
    idxProd = find(strcmp({AScarbon.ID}, prodName));

    sceneConcentrationPerUnitArea = MILIGRAM_TO_GRAM.*AScarbon(idxProd).sceneConcentrationPerUnitArea_lit;
    sceneStock = MILIGRAM_TO_GIGAGRAM.*AScarbon(idxProd).sceneStock_lit;                   
    time = AScarbon(idxProd).timeVectorCarbon;

    % Perform linear regression to estimate trends
    nonNanIndices = ~isnan(sceneStock);
    [pstock,~] = polyfit(datenum(time(nonNanIndices)),sceneStock(nonNanIndices),1);
    trendLineStock = polyval(pstock,datenum(time));
    nonNanIndices = ~isnan(sceneConcentrationPerUnitArea);
    [pconc,~] = polyfit(datenum(time(nonNanIndices)),sceneConcentrationPerUnitArea(nonNanIndices),1);
    trendLineConcentrationPerUnitArea = polyval(pconc,datenum(time));

    % Choose start and end date
    startDate = '1998-01-01';
    endDate = max(time);
    
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.47],'Color','w') 
    haxis = zeros(2,1);
    
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
            scatter(haxis(iSubplot),time,sceneStock,6,'r','filled'); hold on;
            plot(haxis(iSubplot),time,trendLineStock,'Color','b','LineWidth',1.5); hold on;
        else
            scatter(haxis(iSubplot),time,sceneConcentrationPerUnitArea,6,'r','filled'); hold on;
            plot(haxis(iSubplot),time,trendLineConcentrationPerUnitArea,'Color','b','LineWidth',1.5); hold on;
        end
        
        if (iSubplot == 1)
            ylim([yMin1 yMax1])
        else
            ylim([yMin2 yMax2])
        end
            
        xlim([startDate endDate])
%         currentYLim = ylim;
%         if (iSubplot == 2 && currentYLim(2) < 10)
%             ytickformat('%.1f');
%         end

        if (iSubplot == 1)
            ylabel(ylabelStr1)
            title(titleStr)
        else
            ylabel(ylabelStr2)
        end
        if (iSubplot == 2)
            xlabel('Time')
        end
        grid on
        box on

    end % iSubplot

    lg = legend('Data','Trend line');
    lg.Position(1) = 0.87; lg.Position(2) = 0.835; 
    lg.Orientation = 'vertical';
    lg.ItemTokenSize = [10,20];
    set(lg,'Box','on','Color','w')

    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_carbon_timeseries_trends_',prodName,'.png')),'Resolution',600)
    
end % iDataset

end % plotTimeSeriesCarbonStockProperties

% *************************************************************************

function plotHistogramsOfStocks(...
    sceneStockDistribByMonthAndYear_monthlyMean,...
    sceneConcentrationPerUnitAreaByMonthAndYear_monthlyMean,...
    yearVectorData,prodName)

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
barColors = flipud(brewermap(numel(yearVectorData),'*Oranges'));

switch prodName
    case 'cmems_cphyto'
        titleStr1 = sprintf('Scene monthly mean CMEMS Cphyto stock (Gg C)');
        titleStr2 = sprintf('Scene monthly mean CMEMS Cphyto concentration per unit area (g C m-2)');
        yMaxSceneConcArea = 4;
        yMaxSceneStock = 40;
    case 'bicep_poc_4km'
        titleStr1 = sprintf('Scene monthly mean BICEP POC stock (Gg C)');
        titleStr2 = sprintf('Scene monthly mean BICEP POC concentration per unit area (g C m-2)');
        yMaxSceneConcArea = 15;
        yMaxSceneStock = 150;
    case 'bicep_npp_9km'
        titleStr1 = sprintf('Scene monthly integrated BICEP NPP stock (Gg C)');
        titleStr2 = sprintf('Scene monthly integrated BICEP NPP concentration per unit area (g C m-2)');
        yMaxSceneConcArea = 30;
        yMaxSceneStock = 350;
    case 'bicep_cphyto_9km'
        titleStr1 = sprintf('Scene monthly mean BICEP Cphyto stock (Gg C)');
        titleStr2 = sprintf('Scene monthly mean BICEP Cphyto concentration per unit area (g C m-2)');
        yMaxSceneConcArea = 4;
        yMaxSceneStock = 40;
end
    
for iFigure = 1:2

    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.65],'Color','w')

    if (iFigure == 1)
        data = sceneStockDistribByMonthAndYear_monthlyMean;
        ylabelStr = 'Gg C';
        titleStr = titleStr1;
        yMax = yMaxSceneStock;
    else
        data = sceneConcentrationPerUnitAreaByMonthAndYear_monthlyMean;
        ylabelStr = 'g C m^{-2}';
        titleStr = titleStr2;
        yMax = yMaxSceneConcArea;
    end

    for iMonth = 1:12

        haxis(iMonth) = subaxis(6,2,iMonth,'Spacing',0.01,'Padding',0.02,'Margin', 0.04);
        ax(iMonth).pos = get(haxis(iMonth),'Position');
        ax(iMonth).pos(1) = ax(iMonth).pos(1) + 0.025;
        set(haxis(iMonth),'Position',ax(iMonth).pos)

        hbar = bar(haxis(iMonth),data(iMonth,:),'BarWidth',0.80,'FaceColor','flat');
        hbar.CData(:,:) = barColors;

        ylim([0 yMax])
        nYTicks = 3;
        yticks(linspace(0, yMax, nYTicks))
        yticklabels(linspace(0, yMax, nYTicks));
        ytickformat('%.0f');


        xlim([0.5 numel(yearVectorData)+0.5])
        xticks(1:numel(yearVectorData));
        if (iMonth == 11 || iMonth == 12)
            xlabel('Year');
            xticklabels(yearVectorData);
            xtickangle(90); 
        else
            xticklabels([]);
        end
        axh = gca;
        axh.XAxis.FontSize = 8; 

        title(months(iMonth),'FontSize',12);

    end

    % Create a wider title spanning both subplots
    annotation('textbox',[0 0.5 1 0.5],... 
        'String',titleStr,...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',14,'FontWeight','bold');

    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_carbon_distrib_month_and_year_',prodName,'_',iFigure,'.png')),'Resolution',600)

end % iFigure
    
end % plotHistogramsOfStocks

% *************************************************************************

function plotScenesOfStocks(...
    pixelDepthIntegratedStock_monthlyMean,latVector,lonVector,yearVectorData,...
    areaStudy,enduranceSite,prodName,videoObject)

cmap = brewermap(1000,'*RdYlBu');
months = {'January','February','March','April','May','June','July','August',...
    'September','October','November','December'};
yearLabels = num2cell(yearVectorData(1):yearVectorData(end));

nRows = 5;
nCols = 5;
maxSubplots = nRows*nCols;

switch prodName
    case 'cmems_cphyto'
        cbStr = sprintf('Monthly mean CMEMS Cphyto stock (Gg C)');
        yMinPixelStock = 0.010;
        yMaxPixelStock = 0.10;
    case 'bicep_poc_4km'
        cbStr = sprintf('Monthly mean BICEP POC stock (Gg C)');
        yMinPixelStock = 0.010;
        yMaxPixelStock = 0.10;
    case 'bicep_npp_9km'
        cbStr = sprintf('Monthly integrated BICEP NPP stock (Gg C)');
        yMinPixelStock = 0.1;
        yMaxPixelStock = 1.0;
    case 'bicep_cphyto_9km'
        cbStr = sprintf('Monthly mean BICEP Cphyto stock (Gg C)');
        yMinPixelStock = 0.010;
        yMaxPixelStock = 0.10;
end
    
for iMonth = 1:12

    figure()               
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.75],'Color','w') 
    haxis = zeros(nRows,nCols);

    for iYear = 1:numel(yearLabels) 

        % If the iteration exceeds 25, stop plotting
        if iYear > maxSubplots
            break;
        end

        haxis(iYear) = subaxis(nRows,nCols,iYear,'Spacing',0.015,'Padding',0.010,'Margin',0.04);
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

        thisMonthAndYear = squeeze(pixelDepthIntegratedStock_monthlyMean(:,:,iMonth,iYear));
        thisMonthAndYear(thisMonthAndYear == 0) = NaN;
        m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy);
        m_pcolor(lonVector,latVector,thisMonthAndYear) 
        shading flat
        colormap(cmap)
        caxis([yMinPixelStock yMaxPixelStock])
        set(gca,'color','w'); 

        % Show axis labels for border subplots
        if (iYear == 1 || iYear == 6 || iYear == 11 || iYear == 16)
            m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[]);
        elseif (iYear == 22 || iYear == 23 || iYear == 24 || iYear == 25)
            m_grid('linewi',1,'tickdir','out','FontSize',8,'yticklabel',[]);
        elseif (iYear == 21)
            m_grid('linewi',1,'tickdir','out','FontSize',8)
        else
            m_grid('linewi',1,'tickdir','out','FontSize',8,'xticklabel',[],'yticklabel',[]);
        end

        if (iYear == 5)
            cb = colorbar('Location','eastoutside');
            cb.Position(1) = cb.Position(1) + 0.06;
            cb.Position(2) = cb.Position(2) - 0.67;
            cb.Position(3) = 0.0150; % WIDTH
            cb.Position(4) = 0.70; % LENGTH
            cb.Label.String = cbStr;
            cb.FontSize = 9;
%             cb.Ruler.TickLabelFormat = '%.1f';
        end

        m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
        title(yearLabels{iYear},'FontSize',9)
        box on
        hold on

    end
    
    hold off

    % Add a common title for the entire group of subplots
    annotation('textbox',[0.1 0.905 0.8 0.1],'String',months{iMonth},...
        'EdgeColor','none','HorizontalAlignment','center','FontSize',14);

    % Save to file
    exportgraphics(gcf,fullfile('.','figures',...
        strcat('oc_carbon_total_month_',num2str(iMonth),'_',prodName,'.png')),'Resolution',600)
 
    % Write the frame to the video
    frame = getframe(gcf);
    writeVideo(videoObject, frame);
    
    % Close the current figure
    close(gcf);
        
end % iMonth

% Close the video file
close(videoObject);

end % plotScenesOfStocks

% *************************************************************************

function [AScarbon] = processCmemsCphytoDataset(AScmems,AScarbon,areaStudy)

% Load Cphyto data
idxProd = find(strcmp({AScmems.ID}, 'mod_bgc_reg_phy'));
cphytoLat = AScmems(idxProd).lat;
cphytoLon = AScmems(idxProd).lon;
cphytoTime = AScmems(idxProd).time;
cphytoDepth = AScmems(idxProd).depth; 
cphyto = AScmems(idxProd).dataset(:,:,:,:).*12.011; % units conversion to match BICEP's units: mmol C m-3 --> mg C m-3

% According to the description of the product
% (https://catalogue.marine.copernicus.eu/documents/PUM/CMEMS-NWS-PUM-004-009-011.pdf),
% the surface level is not interpolated, it is the first model level. 
% It has a thickness of 1 m where the bathymetry is deeper than 50 m. 
% Its thickness is less than 1 m where the bathymetry is shallower than 
% 50 m. Since the bathymetry at the Endurance is around 50 m 
% (https://emodnet.ec.europa.eu/geoviewer/), let's replace the first 
% vertical depth by 1 m and add a 0 m element
cphytoDepth(1) = 1; % m
depthResolutionCarbon = [1;diff(cphytoDepth)]; % add 1 m to account for the difference 0-1 m depth

% Regrid data for our study area (before interpolating!)
areaStudyLat = linspace(areaStudy.MBRy(1),areaStudy.MBRy(2),numel(cphytoLat))';
areaStudyLon = linspace(areaStudy.MBRx(1),areaStudy.MBRx(2),numel(cphytoLon))';
[cphyto_regridded, ~] = regridAndInterpolateData(cphyto,cphytoLat,cphytoLon,...
    datenum(cphytoTime),cphytoDepth,areaStudyLat,areaStudyLon);

% Interpolate missing values within depth intervals.
% To accurately integrate carbon stocks over depth, we must avoid water
% columns with missing (NaN) values within depth intervals that contain data.
% Therefore, we interpolate through these missing values.
% For those water columns where all pixel values are NaN, we will
% not interpolate latitude x longitude x time as those are 
% missing/invalid water columns. Notice we cannot apply Racault
% interpolation as there is the danger to interpolate to those depths 
% where there should not be data.
cphyto_interp = cphyto_regridded;
nResolvedDepthLevels = zeros(numel(cphytoLat),numel(cphytoLon),numel(cphytoTime));
for iTimeStep = 1:numel(cphytoTime)
    for iLon = 1:numel(cphytoLon)
        for iLat = 1:numel(cphytoLat)
            % Extract carbon
            waterColumnCphyto = squeeze(cphyto_interp(iLat,iLon,iTimeStep,:));
            if (sum(~isnan(waterColumnCphyto)) > 0)
                % Compute depth levels
                nonNanIndices = find(~isnan(waterColumnCphyto));
                lastNonNanPosition = nonNanIndices(end);
                nResolvedDepthLevels(iLat,iLon,iTimeStep) = lastNonNanPosition;
                % Check if any NaN values are in between non-NaN values to
                % interpolate through those
                nanIndices = find(isnan(waterColumnCphyto));
                hasNanBetweenNonNan = any(diff(nanIndices) > 1);
                if (hasNanBetweenNonNan) 
                    % Interpolate using a window of 3 pixels and calculating
                    % the mean (same rationale as in Racault et al. 2014)
                    values = waterColumnCphyto(1:lastNonNanPosition);
                    interpolatedValues = values;
                    for i = 2:(lastNonNanPosition - 1)
                        window = values(i - 1:i + 1);
                        if isnan(interpolatedValues(i))
                            if any(~isnan(window))
                                interpolatedValues(i) = mean(window, 'omitnan');
                            end
                        end
                    end
                    cphyto_interp(iLat,iLon,iTimeStep,1:lastNonNanPosition) = interpolatedValues;
                end
            end
        end % iLat
    end % iLon
end % iTimeStep

% Find latest position occupied in the structure AScarbon
populatedIDs = ~cellfun(@isempty, {AScarbon.ID});
iLatestPositionOccupied = sum(populatedIDs);
iPositionFree = iLatestPositionOccupied + 1;

% Save information into output array
AScarbon(iPositionFree).ID = 'cmems_cphyto';
AScarbon(iPositionFree).carbonDataArray = cphyto;
AScarbon(iPositionFree).carbonDataArrayRegriddedToAreaStudy = cphyto_regridded;
AScarbon(iPositionFree).carbonDataArrayInt = cphyto_interp;
AScarbon(iPositionFree).latVectorCarbon = cphytoLat;
AScarbon(iPositionFree).lonVectorCarbon = cphytoLon;
AScarbon(iPositionFree).latVectorAreaStudy = areaStudyLat;
AScarbon(iPositionFree).lonVectorAreaStudy = areaStudyLon;
AScarbon(iPositionFree).timeVectorCarbon = cphytoTime;
AScarbon(iPositionFree).depthVectorCarbon = cphytoDepth;
AScarbon(iPositionFree).nResolvedDepthLevels = nResolvedDepthLevels;
AScarbon(iPositionFree).depthResolutionCarbon = depthResolutionCarbon;

end % processCmemsCphytoDataset

% *************************************************************************

function [AScarbon] = processBicepDatasets(ASbicep,bicepTimeSeriesDatasetNames,...
    AScarbon,areaStudy)

% Find latest position occupied in the structure AScarbon
populatedIDs = ~cellfun(@isempty, {AScarbon.ID});
iLatestPositionOccupied = sum(populatedIDs);
iPositionFree = iLatestPositionOccupied + 1;

% Find idx to depth integration products
idxMldMimoc = find(strcmp({ASbicep.ID}, 'mld_mimoc_poc'));
idxZeuCmems = find(strcmp({ASbicep.ID}, 'cmems_zeu'));

% Load BICEP's carbon data
for iDataset = 1:size(bicepTimeSeriesDatasetNames, 1)
    
    datasetName = bicepTimeSeriesDatasetNames{iDataset,2};
    idxCarbonProd = find(strcmp({ASbicep.ID}, datasetName));
    
    carbonArrayLat = ASbicep(idxCarbonProd).lat;
    carbonArrayLon = ASbicep(idxCarbonProd).lon;
    carbonArrayTime = ASbicep(idxCarbonProd).time;
    nYears = numel(unique(year(carbonArrayTime)));
    
    % Perform dataset-specific processing
    if (strcmp(datasetName,'bicep_poc_4km'))

        carbonArray = ASbicep(idxCarbonProd).dataset;
        
        depthIntegrationArray_tmp = ASbicep(idxMldMimoc).dataset;
        depthIntegrationArrayLat = ASbicep(idxMldMimoc).lat;
        depthIntegrationArrayLon = ASbicep(idxMldMimoc).lon;
        
        % Time arrangements. The MIMOC array is a monthly climatology, so 
        % we want to repeat the 12 months's data for nYears so that it
        % matches the time frequency of the POC data array.
        depthIntegrationArray = repmat(depthIntegrationArray_tmp,[1,1,nYears]);

    elseif (strcmp(datasetName,'bicep_npp_9km'))
        
        carbonArray = ASbicep(idxCarbonProd).dataset;

        depthIntegrationArray_tmp = ASbicep(idxZeuCmems).dataset;
        depthIntegrationArrayLat = ASbicep(idxZeuCmems).lat;
        depthIntegrationArrayLon = ASbicep(idxZeuCmems).lon;
        depthIntegrationArrayTime = ASbicep(idxZeuCmems).time;
        
        % Time arrangements. The zeu array uses a different time frequency 
        % compared to the NPP data array. Let's match the frequency of the
        % zeu array to that of the NPP array.
        startDate = carbonArrayTime(1);
        endDate = carbonArrayTime(end);
        inRange = depthIntegrationArrayTime >= startDate & depthIntegrationArrayTime <= endDate;
        % 1 - Extract data for the range of dates we have carbon data for
        depthIntegrationArray_inrange = depthIntegrationArray_tmp(:,:,inRange);
        depthIntegrationArrayTime_inrange = depthIntegrationArrayTime(inRange);
        % 2 - Group data by years and months
        [groups,groupIds] = findgroups(year(depthIntegrationArrayTime_inrange),...
            month(depthIntegrationArrayTime_inrange));
        % 3 - Calculate the mean zeu value for each group (year-month combination)
        depthIntegrationArray = NaN(numel(depthIntegrationArrayLat),numel(depthIntegrationArrayLon),12*nYears);
        for iLat = 1:numel(depthIntegrationArrayLat)
            for iLon = 1:numel(depthIntegrationArrayLon)
                localArray = squeeze(depthIntegrationArray_inrange(iLat,iLon,:));
                localArray_monthlymean = splitapply(@mean, localArray, groups);
                depthIntegrationArray(iLat,iLon,:) = localArray_monthlymean;
            end
        end

    elseif (strcmp(datasetName,'bicep_cphyto_9km'))
            
        idxCphytoVars = find(contains(ASbicep(idxCarbonProd).varNames, 'C_'));
        carbonArray = ASbicep(idxCarbonProd).dataset(:,:,:,idxCphytoVars);

        % The MLD product is included in the bicep_cphyto_9km dataset, it
        % has the same time arrangement
        idxMldVar = find(contains(ASbicep(idxCarbonProd).varNames, 'mld'));
        depthIntegrationArray = ASbicep(idxCarbonProd).dataset(:,:,:,idxMldVar);
        depthIntegrationArrayLat = carbonArrayLat;
        depthIntegrationArrayLon = carbonArrayLon;
            
    end

    % Regrid data carbonArray and depthArray for our study area 
    % (before interpolating!). To properly calculate carbon stocks, we cannot 
    % have missing data. Thus, apply interpolation in latitude, longitude and time to fill in gaps.
    
    areaStudyLat = linspace(areaStudy.MBRy(1),areaStudy.MBRy(2),numel(carbonArrayLat))';
    areaStudyLon = linspace(areaStudy.MBRx(1),areaStudy.MBRx(2),numel(carbonArrayLon))';

    [carbonArray_regridded, carbonArray_interp] = regridAndInterpolateData(...
        carbonArray,carbonArrayLat,carbonArrayLon,(1:12*nYears)',[],areaStudyLat,areaStudyLon);
    [depthIntegrationArray_regridded, depthIntegrationArray_interp] =...
        regridAndInterpolateData(depthIntegrationArray,depthIntegrationArrayLat,...
        depthIntegrationArrayLon,(1:12*nYears)',[],areaStudyLat,areaStudyLon);

    % Save information into output array
    AScarbon(iPositionFree).ID = char(datasetName);
    AScarbon(iPositionFree).carbonDataArray = carbonArray;
    AScarbon(iPositionFree).carbonDataArrayRegriddedToAreaStudy = carbonArray_regridded;
    AScarbon(iPositionFree).carbonDataArrayInt = carbonArray_interp;
    AScarbon(iPositionFree).latVectorCarbon = carbonArrayLat;
    AScarbon(iPositionFree).lonVectorCarbon = carbonArrayLon;
    AScarbon(iPositionFree).latVectorAreaStudy = areaStudyLat;
    AScarbon(iPositionFree).lonVectorAreaStudy = areaStudyLon;
    AScarbon(iPositionFree).timeVectorCarbon = carbonArrayTime;
    AScarbon(iPositionFree).depthIntegrationDataArray = depthIntegrationArray;
    AScarbon(iPositionFree).depthIntegrationDataArrayRegriddedToAreaStudy = depthIntegrationArray_regridded;
    AScarbon(iPositionFree).depthIntegrationDataArrayInt = depthIntegrationArray_interp;
    AScarbon(iPositionFree).latVectorDepthIntegrationArray = depthIntegrationArrayLat;
    AScarbon(iPositionFree).lonVectorDepthIntegrationArray = depthIntegrationArrayLon;    
    iPositionFree = iPositionFree + 1;
    
end % iDataset

end % processBicepDatasets

% *************************************************************************

% end % calculateAndPlotCarbonStocksFromOceanColourProducts
