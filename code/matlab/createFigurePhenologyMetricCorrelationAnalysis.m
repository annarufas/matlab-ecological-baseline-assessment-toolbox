
clc
clear all

fullPathMainDir         = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/';
fullPathDataCMEMSdir    = strcat(fullPathMainDir,'CMEMS_data/data_timeseries_new/');
fullPathPlotsDir        = strcat(fullPathMainDir,'plots/phenology/');
fullPathAreaStudyCoords = strcat(fullPathMainDir,'coords_boundarybox/bbox_50km');
fullPathEnduranceCoords = strcat(fullPathMainDir,'coords_endurance/endurance_2023');
addpath(genpath(strcat(fullPathMainDir)))

areaStudy = m_shaperead(fullPathAreaStudyCoords); % lat/lon coordinates
enduranceSite = m_shaperead(fullPathEnduranceCoords); % lat/lon coordinates

% Box of 100 x 100 km around the Endurance site
minLatAreaStudy = areaStudy.MBRy(1); % 53.758
maxLatAreaStudy = areaStudy.MBRy(2); % 54.656
minLonAreaStudy = areaStudy.MBRx(1); % 0.265
maxLonAreaStudy = areaStudy.MBRx(2); % 1.801

% First, average scenes. Second, regrid in time. In this order please!!!

%% Load datasets

% Clotophyll a
load(strcat(fullPathMainDir,'continuousSatelliteChlaFiveDay4km.mat'),...
    'continuousSatelliteChlaFiveDay4km','timeVectorCompleteFiveDay4km',...
    'latSatellite','lonSatellite')

% POC, Cphyto and NPP
load(strcat(fullPathMainDir,'carbonStocks.mat'))

% CMEMS biogeochemical and physical model reanalysis products
load(strcat(fullPathDataCMEMSdir,'OC_cmems_bbox50km.mat'))

% Suspended particulate matter
load(strcat(fullPathMainDir,'spmSilva.mat'))

% Phenology metrics
load(strcat(fullPathMainDir,'phenologyMetrics.mat'))

% Climate indices (from https://psl.noaa.gov/data/climateindices/)
filePathNaoIndex = strcat(fullPathMainDir,'NAO_index.txt');
data = importdata(filePathNaoIndex);
valuesNaoIndex = data(:,2);

filePathMeiIndex = strcat(fullPathMainDir,'MEI_index.txt');
data = importdata(filePathMeiIndex);
valuesMeiIndex = data(:,2);

filePathAmoIndex = strcat(fullPathMainDir,'AMO_index.txt');
data = importdata(filePathAmoIndex);
valuesAmoIndex = data(:,2);

filePathGloTemp = strcat(fullPathMainDir,'globalMeanTemperature.txt');
data = importdata(filePathGloTemp);
valuesGloTempIndex = data(:,2);

%% Process climate indexes

timeClimateIndexes = datetime('1998-01-15'):calmonths(1):datetime('2023-12-15');

% Use a table to arrange the scene data and average per year
TS = table(valuesNaoIndex,valuesMeiIndex,valuesAmoIndex,valuesGloTempIndex,...
    timeClimateIndexes','VariableNames',{'nao','mei','amo','temp','date'});
TS.year = year(TS.date);
yearVectorData = min(TS.year):1:max(TS.year);
naoAnnualMean = NaN(numel(yearVectorData),1); 
meiAnnualMean = NaN(numel(yearVectorData),1); 
amoAnnualMean = NaN(numel(yearVectorData),1); 
tempAnnualMean = NaN(numel(yearVectorData),1);

for iYear = 1:numel(yearVectorData)
    thisYearNao = TS.nao(TS.year == yearVectorData(iYear));
    thisYearMei = TS.mei(TS.year == yearVectorData(iYear));
    thisYearAmo = TS.amo(TS.year == yearVectorData(iYear));
    thisYearTemp = TS.temp(TS.year == yearVectorData(iYear));
    naoAnnualMean(iYear) = mean(thisYearNao,'omitnan'); 
    meiAnnualMean(iYear) = mean(thisYearMei,'omitnan'); 
    amoAnnualMean(iYear) = mean(thisYearAmo,'omitnan'); 
    tempAnnualMean(iYear) = mean(thisYearTemp,'omitnan'); 
end

% Adjust years
iYearEnd = find(yearVectorData == 2020);
valuesNaoAdj = naoAnnualMean(1:iYearEnd);
valuesMeiAdj = meiAnnualMean(1:iYearEnd);
valuesAmoAdj = amoAnnualMean(1:iYearEnd);
valuesTempAdj = tempAnnualMean(1:iYearEnd);

%% Process CMEMS products

% Time query points for interpolation
qTimeCmems = datenum(datetime('1993-01-15'):calmonths(1):datetime('2022-12-15'));

listIdCmemsProd = {'mod_bgc_reg_kd',...
                   'mod_bgc_reg_no3',...
                   'mod_bgc_reg_po4',...
                   'mod_bgc_reg_o2',...
                   'mod_bgc_reg_ph',...
                   'mod_bgc_reg_pco2',...
                   'mod_phy_reg_mld',...
                   'mod_phy_reg_sal',...
                   'mod_phy_reg_temp',...
                   'mod_phy_reg_ssh',...
                   'mod_phy_reg_velo'};
               
nCmemsProdFinalList = numel(listIdCmemsProd) + 1; % to open up velocity
             
iProcessedProduct = 1;

% Regrid (lat, lon, time), fill in gaps and average per scene
for iProd = 1:numel(listIdCmemsProd)
    
    idProd = listIdCmemsProd{iProd};
    posIdProd = find(strcmp({OC_cmems.ID}, idProd));
    cmemsTime = OC_cmems(posIdProd).time;
    cmemsLat = OC_cmems(posIdProd).lat;
    cmemsLon = OC_cmems(posIdProd).lon;
    nTimeSteps = length(cmemsTime);
    nLatPixelsCmems = length(cmemsLat);
    nLonPixelsCmems = length(cmemsLon);
    qLatVectorAreaStudy = linspace(minLatAreaStudy,maxLatAreaStudy,nLatPixelsCmems)';
    qLonVectorAreaStudy = linspace(minLonAreaStudy,maxLonAreaStudy,nLonPixelsCmems)';
    
    cmemsProd = OC_cmems(posIdProd).dataset;

    % Initialise array
    if (iProd == 1)
        cmemsProcessedProd = NaN(nLatPixelsCmems,nLonPixelsCmems,nTimeSteps,nCmemsProdFinalList);
        cmemsFilledProduct = NaN(nLatPixelsCmems,nLonPixelsCmems,nTimeSteps,nCmemsProdFinalList);
    end
    
    % We want to calculate Zeu from kd
    if (strcmp(OC_cmems(posIdProd).ID,'mod_bgc_reg_kd'))
        
        cmemsKd = cmemsProd;
        cmemsZeu = NaN(size(cmemsProd,1:3));
        for iTimeStep = 1:nTimeSteps
            for iLon = 1:nLonPixelsCmems
                for iLat = 1:nLatPixelsCmems
                    % Notice we use "1" to take only first value provided
                    if (cmemsKd(iLat,iLon,iTimeStep,1) > 0 && ~isnan(cmemsKd(iLat,iLon,iTimeStep,1)))
                        cmemsZeu(iLat,iLon,iTimeStep) = 4.6./cmemsKd(iLat,iLon,iTimeStep,1); % standard formula to calculate zeu
                    end
                end
            end
        end
        cmemsProcessedProd(:,:,:,iProcessedProduct) = cmemsZeu;
        iProcessedProduct = iProcessedProduct + 1;
    
    % These products don't need any processing (unique depth) 
    elseif (strcmp(OC_cmems(posIdProd).ID,'mod_phy_reg_mld') ||...
            strcmp(OC_cmems(posIdProd).ID,'mod_phy_reg_ssh') ||...
            strcmp(OC_cmems(posIdProd).ID,'mod_bgc_reg_pco2')) 
    
        cmemsProcessedProd(:,:,:,iProcessedProduct) = cmemsProd;
        iProcessedProduct = iProcessedProduct + 1;
        
    % Products with depth levels    
    elseif (strcmp(OC_cmems(posIdProd).ID,'mod_bgc_reg_no3') ||...
            strcmp(OC_cmems(posIdProd).ID,'mod_bgc_reg_po4') ||...
            strcmp(OC_cmems(posIdProd).ID,'mod_bgc_reg_o2') ||...
            strcmp(OC_cmems(posIdProd).ID,'mod_bgc_reg_ph') ||...
            strcmp(OC_cmems(posIdProd).ID,'mod_phy_reg_sal') ||...
            strcmp(OC_cmems(posIdProd).ID,'mod_phy_reg_temp'))
        
        depthAveragedProd = NaN(size(cmemsProd,1:3));
        for iTimeStep = 1:nTimeSteps
            for iLon = 1:nLonPixelsCmems
                for iLat = 1:nLatPixelsCmems
                    depthAveragedProd(iLat,iLon,iTimeStep) = mean(cmemsProd(iLat,iLon,iTimeStep,:),'omitnan');
                end
            end
        end
        cmemsProcessedProd(:,:,:,iProcessedProduct) = depthAveragedProd;
        iProcessedProduct = iProcessedProduct + 1;
 
    % Split velocities    
    elseif (strcmp(OC_cmems(posIdProd).ID,'mod_phy_reg_velo'))
        
        for iVeloProd = 1:2
            depthAveragedProd = NaN(size(cmemsProd,1:3));
            for iTimeStep = 1:nTimeSteps
                for iLon = 1:nLonPixelsCmems
                    for iLat = 1:nLatPixelsCmems
                        depthAveragedProd(iLat,iLon,iTimeStep) =... 
                            mean(squeeze(cmemsProd(iLat,iLon,iTimeStep,:,iVeloProd)),'omitnan');
                    end
                end
            end
            cmemsProcessedProd(:,:,:,iProcessedProduct) = depthAveragedProd;
            iProcessedProduct = iProcessedProduct + 1;
        end
        
    end
    
    % New grid (query points for interpolation)
    [qX,qY,qT] = ndgrid(qLatVectorAreaStudy,qLonVectorAreaStudy,datenum(cmemsTime));
    % Actual grid
    [X,Y,T] = ndgrid(cmemsLat,cmemsLon,datenum(cmemsTime));
    
    % Regrid and fill
    if (iProd < numel(listIdCmemsProd))
        % Interpolant
        F = griddedInterpolant(X,Y,T,cmemsProcessedProd(:,:,:,iProd));
        % Regridded product
        cmemsRegriddedProductToAreaStudy = F(qX,qY,qT); 
        cmemsFilledProduct(:,:,:,iProd) = Racault2014interpolation(cmemsRegriddedProductToAreaStudy);
    elseif (iProd == numel(listIdCmemsProd)) % to account for the two velocity products
        F = griddedInterpolant(X,Y,T,cmemsProcessedProd(:,:,:,iProd));
        cmemsRegriddedProductToAreaStudy = F(qX,qY,qT); 
        cmemsFilledProduct(:,:,:,iProd) = Racault2014interpolation(cmemsRegriddedProductToAreaStudy);
        F = griddedInterpolant(X,Y,T,cmemsProcessedProd(:,:,:,iProd+1));
        cmemsRegriddedProductToAreaStudy = F(qX,qY,qT); 
        cmemsFilledProduct(:,:,:,iProd+1) = Racault2014interpolation(cmemsRegriddedProductToAreaStudy);
    end    
     
end

% Average per scene
cmemsSceneAverage = NaN(numel(cmemsTime),nCmemsProdFinalList);
for iProd = 1:nCmemsProdFinalList
    for iTimeStep = 1:numel(cmemsTime)
        thisScene = reshape(squeeze(cmemsFilledProduct(:,:,iTimeStep,iProd)),[],1);
        cmemsSceneAverage(iTimeStep,iProd) = mean(thisScene,'omitnan');
    end
end

% Use a table to arrange the scene data and average per year
for iProd = 1:nCmemsProdFinalList
    thisProd = cmemsSceneAverage(:,iProd);
    TS = table(thisProd,cmemsTime,'VariableNames',{'data','date'});
    TS.year = year(TS.date);
    yearVectorData = min(TS.year):1:max(TS.year);
    if (iProd == 1) % initialise
        prodAnnualMean = NaN(numel(yearVectorData),nCmemsProdFinalList); 
    end
    for iYear = 1:numel(yearVectorData)
        thisYearData = TS.data(TS.year == yearVectorData(iYear));
        prodAnnualMean(iYear,iProd) = mean(thisYearData,'omitnan'); 
    end
end

% Adjust years
iYearEnd = find(yearVectorData == 2020);
iYearStart = find(yearVectorData == 1998);
valuesCmemsAdj = prodAnnualMean(iYearStart:iYearEnd,:);

%% Process C stock products

% These products are already lat x lon regridded to my area and with the
% right time vector

for iVar = 1:3
    if (iVar == 1)
        posId = find(strcmp({carbonStocks.ID}, 'bicep_poc_4km'));
    elseif (iVar == 2)    
        posId = find(strcmp({carbonStocks.ID}, 'bicep_cphyto_9km'));
    elseif (iVar == 3)
        posId = find(strcmp({carbonStocks.ID}, 'bicep_npp_9km'));
    end
    
    if (iVar == 1) % initialise
        timeBicep = carbonStocks(posId).timeVectorCarbon; % same for all
        carbonStockSceneAverage = NaN(numel(timeBicep),1);
        carbonConcentrationPerUnitAreaSceneAverage = NaN(numel(timeBicep),1);
    end

    carbonStockSceneAverage(:,iVar) = 1e-12.*carbonStocks(posId).sceneStock_lit; % mg C --> Mg C
    carbonConcentrationPerUnitAreaSceneAverage(:,iVar) = 1e-3.*carbonStocks(posId).sceneConcentrationPerUnitArea_lit; % mg C m-2 --> g C m-2

end

% Use a table to arrange the scene data and average per year
TS = table(carbonConcentrationPerUnitAreaSceneAverage(:,1),...
    carbonConcentrationPerUnitAreaSceneAverage(:,2),...
    carbonConcentrationPerUnitAreaSceneAverage(:,3),...
    timeBicep,'VariableNames',{'poc','cphyto','npp','date'});
TS.year = year(TS.date);
yearVectorData = min(TS.year):1:max(TS.year);
stockPocAnnualMean = NaN(numel(yearVectorData),1); 
stockCphytoAnnualMean = NaN(numel(yearVectorData),1); 
stockNppAnnualMean = NaN(numel(yearVectorData),1); 
for iYear = 1:numel(yearVectorData)
    thisYearPocData = TS.poc(TS.year == yearVectorData(iYear));
    thisYearCphytoData = TS.cphyto(TS.year == yearVectorData(iYear));
    thisYearNppData = TS.npp(TS.year == yearVectorData(iYear));
    stockPocAnnualMean(iYear) = mean(thisYearPocData,'omitnan'); 
    stockCphytoAnnualMean(iYear) = mean(thisYearCphytoData,'omitnan'); 
    stockNppAnnualMean(iYear) = mean(thisYearNppData,'omitnan'); 
end

% Adjust years
iYearEnd = find(yearVectorData == 2020);
valuesPocStockAdj = stockPocAnnualMean(1:iYearEnd);
valuesCphytoStockAdj = stockCphytoAnnualMean(1:iYearEnd);
valuesNppStockAdj = stockNppAnnualMean(1:iYearEnd);

%% Process chla product

% Chla goes 1 Jan 1998 to 27 Dec 2023

% Average per scene
chlaSceneAverage = NaN(numel(timeVectorCompleteFiveDay4km),1);
for iTimeStep = 1:length(timeVectorCompleteFiveDay4km)
    thisScene = reshape(continuousSatelliteChlaFiveDay4km(:,:,iTimeStep),[],1);
    chlaSceneAverage(iTimeStep) = mean(thisScene,'omitnan');
end

% Use a table to arrange the scene data and average per year
TS = table(chlaSceneAverage,timeVectorCompleteFiveDay4km,'VariableNames',{'data','date'});
TS.year = year(TS.date);
yearVectorData = min(TS.year):1:max(TS.year);
chlaAnnualMean = NaN(numel(yearVectorData),1); 
for iYear = 1:numel(yearVectorData)
    thisYearData = TS.data(TS.year == yearVectorData(iYear));
    chlaAnnualMean(iYear) = mean(thisYearData,'omitnan'); 
end

% Adjust years
iYearEnd = find(yearVectorData == 2020);
valuesChlaAdj = chlaAnnualMean(1:iYearEnd);

%% Process SPM product

[~,~,nYears] = size(spmAnnualClimatologyMean);
%yearsVector = 1998:1:2015;
        
% Create scene average    
spmSceneAveraged = NaN(nYears,1); 
for iYear = 1:nYears
    thisScene = reshape(spmAnnualClimatologyMean(:,:,iYear),[],1);
    spmSceneAveraged(iYear) = mean(thisScene,'omitnan');
end   

% Append 5 more years data
runningMean = mean(spmSceneAveraged,'omitnan');
valuesSpmAdj = [spmSceneAveraged(:);runningMean.*ones(5,1)];

%% Process phenology metrics

[~,~,nMetrics,nYears] = size(seasonalityDataFirstBloom);

% Scene average
for iProd = 1:2 % 1st bloom, 2nd bloom
    if (iProd == 1)
        dataArray = seasonalityDataFirstBloom;
    elseif (iProd == 2)
        dataArray = seasonalityDataSecondBloom;
    end
    sceneAveraged = NaN(nMetrics,nYears);
    for iYear = 1:nYears
        for iMetric = 1:nMetrics
            thisScene = reshape(dataArray(:,:,iMetric,iYear),[],1);
            sceneAveraged(iMetric,iYear) = mean(thisScene,'omitnan');
        end
    end
    if (iProd == 1)
        firstBloomSceneAveraged = sceneAveraged;
    elseif (iProd == 2)
        secondBloomSceneAveraged = sceneAveraged;
    end
end 

% Adjust years
yearsVector = unique(year(phenologyTime));
iYearEnd = find(yearsVector == 2020);
valuesPhenologyFirstBloomAdj = firstBloomSceneAveraged(:,1:iYearEnd);
valuesPhenologySecondBloomAdj = secondBloomSceneAveraged(:,1:iYearEnd);

% Swap metric and time
valuesPhenologyFirstBloomAdj = permute(valuesPhenologyFirstBloomAdj, [2 1]);
valuesPhenologySecondBloomAdj = permute(valuesPhenologySecondBloomAdj, [2 1]);

%% Spearman rank

predictorsData = [valuesChlaAdj,...         % Chla
                  valuesPocStockAdj,...     % POC
                  valuesCphytoStockAdj,...  % Cphyto
                  valuesNppStockAdj,...     % NPP
                  valuesCmemsAdj(:,1),...   % Zeu
                  valuesSpmAdj,...          % SPM
                  valuesCmemsAdj(:,2),...   % NO3
                  valuesCmemsAdj(:,4:5),... % O2, pH
                  valuesCmemsAdj(:,7:end),... % MLD, Sal, T, SSH, u0, v0
                  valuesNaoAdj,...          % NAO
                  valuesMeiAdj,...          % MEI
                  valuesAmoAdj,...          % AMO
                  valuesTempAdj];           % global mean T

responseDataPhenoMainBloom = [valuesPhenologyFirstBloomAdj(:,1),... % tstart
                              valuesPhenologyFirstBloomAdj(:,5),... % tmid
                              valuesPhenologyFirstBloomAdj(:,2),... % tend
                              valuesPhenologyFirstBloomAdj(:,3),... % delta_t
                              valuesPhenologyFirstBloomAdj(:,4),... % chla base
                              valuesPhenologyFirstBloomAdj(:,6),... % chla peak
                              valuesPhenologyFirstBloomAdj(:,7),... % delta_chla
                              valuesPhenologyFirstBloomAdj(:,8),... % rate increase
                              valuesPhenologyFirstBloomAdj(:,9)];   % rate decrease

responseDataPhenoSecondaryBloom = [valuesPhenologySecondBloomAdj(:,1),...
                                   valuesPhenologySecondBloomAdj(:,5),...
                                   valuesPhenologySecondBloomAdj(:,2),...
                                   valuesPhenologySecondBloomAdj(:,3),...
                                   valuesPhenologySecondBloomAdj(:,4),...
                                   valuesPhenologySecondBloomAdj(:,6),...
                                   valuesPhenologySecondBloomAdj(:,7),...
                                   valuesPhenologySecondBloomAdj(:,8),...
                                   valuesPhenologySecondBloomAdj(:,9)];

% Calculate Spearman rank correlation coefficients
[rhoFirstBloom, pvalFirstBloom] = corr(responseDataPhenoMainBloom, predictorsData, 'Type', 'Spearman');
[rhoSecondBloom, pvalSecondBloom] = corr(responseDataPhenoSecondaryBloom, predictorsData, 'Type', 'Spearman');

[nYears,nPredictors] = size(predictorsData);
[nRows,nCols] = size(rhoFirstBloom);

% Plot
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.78 0.30],'Color','w')
haxis = zeros(1,2);

for iSubplot = 1:2
    
    if (iSubplot == 1)
        correlationMatrix = rhoFirstBloom;
        pvalueMatrix = pvalFirstBloom;
    elseif (iSubplot == 2)   
        correlationMatrix = rhoSecondBloom;
        pvalueMatrix = pvalSecondBloom;
    end
    
    haxis(iSubplot) = subaxis(1,2,iSubplot,'Spacing',0.010,'Padding',0.02,'Margin',0.03);
    ax(iSubplot).pos = get(haxis(iSubplot),'Position');
    if (iSubplot == 1)
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1);
    else
        ax(iSubplot).pos(1) = ax(iSubplot).pos(1) - 0.03;
    end
    ax(iSubplot).pos(2) = ax(iSubplot).pos(2) + 0.045;
    set(haxis(iSubplot),'Position',ax(iSubplot).pos)

    imagesc(correlationMatrix(:,:)); hold on;
    colormap(brewermap(20,'*RdBu'));
    caxis([-1 1])

    % Display correlation values in the middle of each cell
    for iRow = 1:nRows
        for iCol = 1:nCols
            thisPvalue = pvalueMatrix(iRow, iCol);
            thisCorrCoeff = correlationMatrix(iRow, iCol);
            % Determine if the correlation value should appear
            isBold = thisPvalue < 0.05;
            if (isBold)
                % Display the correlation value in the middle of the cell
                text(iCol, iRow, num2str(abs(thisCorrCoeff),'%.2f'), ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle',...
                    'FontWeight', 'bold','FontSize',8);
            end
        end
    end
    
    axh = gca;
    set(axh, 'TickLength', [0 0]);

    yticklabels([]);
    xticklabels([]);
    
    box on;
    axh.LineWidth = 2;  % Set the line width according to your preference

end

% Colour bar settings
cb = colorbar('Location','southoutside');
cb.Position(1) = cb.Position(1);
cb.Position(2) = cb.Position(2) - 0.12;
cb.Position(3) = 0.425; % LENGTH
cb.Position(4) = 0.035; % WIDTH

exportgraphics(gcf,fullfile(fullPathPlotsDir,'plot_corr_matrix.png'),'Resolution',600)

