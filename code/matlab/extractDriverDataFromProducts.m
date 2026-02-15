function [S,cmemsSelectedDatasets,bicepSelectedDatasets,chlaSelectedDatasets] =... 
    extractDriverDataFromProducts(AScmems,ASbicep,ASnasa,ASoccci)

% EXTRACTDRIVERDATAFROMPRODUCTS Extract driver data from (i) CMEMS
% biogeochemical and physical reanalysis products, (ii) BICEP carbon products
% products and (iii) chlorophyll a products and average data by scene and
% depth, and organise the data into a structure for easy access and analysis.
%
%   INPUT: 
%       AScmems               - table with time-series data from CMEMS for our area of study
%       ASbicep               - table with time-series data from BICEP for our area of study
%       ASnasa
%       ASoccci
%
%   OUTPUT:
%       S                     - data structure with data extracted from selected datasets
%       cmemsSelectedDatasets - string containing selected product names
%       bicepSelectedDatasets - string containing selected product names
%       chlaSelectedDatasets
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

% Define dataset labels and names for selected products

cmemsSelectedDatasets = {...
    'mod_bgc_reg_diat', {'Diatoms'};
    'mod_bgc_reg_dino', {'Dinoflagellates'};
    'mod_bgc_reg_nano', {'Nanophytoplankton'};
    'mod_bgc_reg_pico', {'Picophytoplankton'};
    'mod_bgc_reg_kd',   {'Zeu'};
    'mod_bgc_reg_no3',  {'Nitrate'};
    'mod_bgc_reg_po4',  {'Phosphate'};
    'mod_bgc_reg_o2',   {'Oxygen'};
    'mod_bgc_reg_ph',   {'pH'};
    'mod_bgc_reg_pco2', {'pCO2'};
    'mod_phy_reg_mld',  {'MLD'};
    'mod_phy_reg_sal',  {'Salinity'};
    'mod_phy_reg_temp', {'Temperature'};
    'mod_phy_reg_ssh',  {'SSH'};
    'mod_phy_reg_velo', {'Eastward velocity component','Northward velocity component'}
};
                                 
bicepSelectedDatasets = {...
    'bicep_npp_9km',    {'NPP'};
    'bicep_poc_4km',    {'POC'};
    'bicep_cphyto_9km', {'Cphyto','Microphytoplankton','Nanophytoplankton','Picophytoplankton'};
};

% Chlorophyll datasets
[chlaStruct,chlaProdLabel,chlaProdName,~] = ...
    extractChlorophyllFromOceanColourProducts(AScmems,ASnasa,ASoccci);
chlaSelectedDatasets = [chlaProdName',chlaProdLabel'];

% Initialise a structure to hold the extracted data
S = struct();  
iDatasetLoc = 1;

%% Get CMEMS datasets and average by scene and depth
           
for iDataset = 1:size(cmemsSelectedDatasets, 1)
    
    datasetName = cmemsSelectedDatasets{iDataset, 1};
    idxDataset = find(strcmp({AScmems.ID}, datasetName));
    
    % Common dataset properties
    if (iDataset == 1)
        cmemsTime = AScmems(idxDataset).time;
        nTimeSteps = length(cmemsTime);
    end

    % We want to calculate Zeu from kd
    if (strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_kd'))
        
        S(iDatasetLoc).ID = char(cmemsSelectedDatasets{iDataset, 2});
        
        [cmemsZeu] = calculateZeuFromCmemsKd(AScmems);
        [cmemsSceneAverage,~] = calcSceneAverage(cmemsZeu,cmemsTime);
        [cmemsSceneAverageYearlyMean,cmemsSceneAverageMonthlyMean] =... 
            calcSceneTimeAveragedData(cmemsSceneAverage,cmemsTime);
        
        [S,iDatasetLoc] = addSceneAverageToDataStruct(S,cmemsTime,...
            cmemsSceneAverage,cmemsSceneAverageYearlyMean,...
            cmemsSceneAverageMonthlyMean,iDatasetLoc);

    % These products don't need any processing (unique depth) 
    elseif (strcmp(AScmems(idxDataset).ID,'mod_phy_reg_mld') ||...
            strcmp(AScmems(idxDataset).ID,'mod_phy_reg_ssh') ||...
            strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_pco2')) 

        S(iDatasetLoc).ID = char(cmemsSelectedDatasets{iDataset, 2});
        
        cmemsDataset = AScmems(idxDataset).dataset; 
        [cmemsSceneAverage,~] = calcSceneAverage(cmemsDataset,cmemsTime);
        [cmemsSceneAverageYearlyMean,cmemsSceneAverageMonthlyMean] =... 
            calcSceneTimeAveragedData(cmemsSceneAverage,cmemsTime);
        [S,iDatasetLoc] = addSceneAverageToDataStruct(S,cmemsTime,...
            cmemsSceneAverage,cmemsSceneAverageYearlyMean,...
            cmemsSceneAverageMonthlyMean,iDatasetLoc);
        
    % Products with depth levels, pondered depth average    
    elseif (strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_diat') ||...
            strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_dino') ||...
            strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_nano') ||...
            strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_pico') ||...
            strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_no3') ||...
            strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_po4') ||...
            strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_o2') ||...
            strcmp(AScmems(idxDataset).ID,'mod_bgc_reg_ph') ||...
            strcmp(AScmems(idxDataset).ID,'mod_phy_reg_sal') ||...
            strcmp(AScmems(idxDataset).ID,'mod_phy_reg_temp'))
        
        S(iDatasetLoc).ID = char(cmemsSelectedDatasets{iDataset, 2});
        
        nDepthLevels = numel(AScmems(idxDataset).depth);
        cmemsDataset = AScmems(idxDataset).dataset; 
        cmemsSceneAverage = calcScenePonderedDepthAverage(cmemsDataset,...
            cmemsTime,nDepthLevels,nTimeSteps);
        [cmemsSceneAverageYearlyMean,cmemsSceneAverageMonthlyMean] =... 
            calcSceneTimeAveragedData(cmemsSceneAverage,cmemsTime);
        [S,iDatasetLoc] = addSceneAverageToDataStruct(S,cmemsTime,...
            cmemsSceneAverage,cmemsSceneAverageYearlyMean,...
            cmemsSceneAverageMonthlyMean,iDatasetLoc);

    % Velocity has two components  
    elseif (strcmp(AScmems(idxDataset).ID,'mod_phy_reg_velo'))
        for i = 1:2
            
            % Access the nested cell array in the second column of the row
            nestedCellArray = cmemsSelectedDatasets{iDataset, 2};
            S(iDatasetLoc).ID = char(nestedCellArray{i});
            
            cmemsDataset = squeeze(AScmems(idxDataset).dataset(:,:,:,:,i));
            nDepthLevels = numel(AScmems(idxDataset).depth);
            cmemsSceneAverage = calcScenePonderedDepthAverage(cmemsDataset,...
                cmemsTime,nDepthLevels,nTimeSteps);
            [cmemsSceneAverageYearlyMean,cmemsSceneAverageMonthlyMean] =... 
                calcSceneTimeAveragedData(cmemsSceneAverage,cmemsTime);
            [S,iDatasetLoc] = addSceneAverageToDataStruct(S,cmemsTime,...
                cmemsSceneAverage,cmemsSceneAverageYearlyMean,...
                cmemsSceneAverageMonthlyMean,iDatasetLoc);

        end
    end
 
end % iDataset

%% Get BICEP's datasets and average by scene
           
for iDataset = 1:size(bicepSelectedDatasets, 1)
    
    datasetName = bicepSelectedDatasets{iDataset, 1};
    idxDataset = find(strcmp({ASbicep.ID}, datasetName));
    
    % Common dataset properties
    if (iDataset == 1)
        bicepTime = ASbicep(idxDataset).time;
        nTimeSteps = length(bicepTime);
    end

    % These products do not need any processing
    if (strcmp(ASbicep(idxDataset).ID,'bicep_poc_4km') ||...
        strcmp(ASbicep(idxDataset).ID,'bicep_npp_9km'))
    
        S(iDatasetLoc).ID = char(bicepSelectedDatasets{iDataset, 2});
        
        bicepDataset = ASbicep(idxDataset).dataset; 
        [bicepSceneAverage,~] = calcSceneAverage(bicepDataset,bicepTime);
        [bicepSceneAverageYearlyMean,bicepSceneAverageMonthlyMean] = calcSceneTimeAveragedData(bicepSceneAverage,bicepTime);
        [S,iDatasetLoc] = addSceneAverageToDataStruct(S,bicepTime,bicepSceneAverage,bicepSceneAverageYearlyMean,bicepSceneAverageMonthlyMean,iDatasetLoc);
        
    elseif (strcmp(ASbicep(idxDataset).ID,'bicep_cphyto_9km'))
        
        % Find indices of plankton variable names from the dataset
        bicepPlkVarNames = {'C_phyto','C_microphyto','C_nanophyto','C_picophyto'};
        prodVarNames = ASbicep(idxDataset).varNames;
        idxsPlkVars = [];
        for i = 1:numel(prodVarNames)
            for j = 1:numel(bicepPlkVarNames)
                % Match indices of plankton variable names using strcmp (case-sensitive)
                if strcmp(prodVarNames{i}, bicepPlkVarNames{j})
                    idxsPlkVars(end + 1) = i; % add the index to the array
                end
            end
        end
        
        for i = 1:numel(bicepPlkVarNames)
            
            % Access the nested cell array in the second column of the row
            nestedCellArray = bicepSelectedDatasets{iDataset, 2};
            S(iDatasetLoc).ID = char(nestedCellArray{i});
   
            bicepDataset = squeeze(ASbicep(idxDataset).dataset(:,:,:,idxsPlkVars(i)));
            bicepSceneAverage = calcSceneAverage(bicepDataset,bicepTime);
            [bicepSceneAverageYearlyMean,bicepSceneAverageMonthlyMean] = calcSceneTimeAveragedData(bicepSceneAverage,bicepTime);
            [S,iDatasetLoc] = addSceneAverageToDataStruct(S,bicepTime,bicepSceneAverage,bicepSceneAverageYearlyMean,bicepSceneAverageMonthlyMean,iDatasetLoc);
 
        end
        
    end
    
end % iDataset

%% Get chla datasets and average by scene
           
for iDataset = 1:size(chlaSelectedDatasets, 1)
    
    thisDatasetName = chlaSelectedDatasets{iDataset, 1};
    chlaTime = chlaStruct.(thisDatasetName).time;

    S(iDatasetLoc).ID = char(chlaSelectedDatasets{iDataset, 2});

    chlaDataset = chlaStruct.(thisDatasetName).dataset;
    [chlaSceneAverage,~] = calcSceneAverage(chlaDataset,chlaTime);
    [chlaSceneAverageYearlyMean,chlaSceneAverageMonthlyMean] = calcSceneTimeAveragedData(chlaSceneAverage,chlaTime);
    [S,iDatasetLoc] = addSceneAverageToDataStruct(S,chlaTime,chlaSceneAverage,chlaSceneAverageYearlyMean,chlaSceneAverageMonthlyMean,iDatasetLoc);
    
end 

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function [S,iDatasetLoc] = addSceneAverageToDataStruct(S,time,sceneAverage,...
    sceneAverageYearlyMean,sceneAverageMonthlyMean,iDatasetLoc)
    
S(iDatasetLoc).sceneAverage = sceneAverage;
S(iDatasetLoc).time = time;
S(iDatasetLoc).yearlyMean = sceneAverageYearlyMean;
S(iDatasetLoc).monthlyMean = sceneAverageMonthlyMean;
iDatasetLoc = iDatasetLoc + 1;

end % addSceneAverageToDataStruct

%%
% for iDataset = 1:numel(cmemsDatasetNames)+1
%     TS = table(S(iDataset).sceneAverage,S(iDataset).time,...
%         'VariableNames',{'sceneaverage','date'});
%     TS.year = year(TS.date);
%     TS.month = month(TS.date);
%     yearVectorData = min(TS.year):1:max(TS.year);
%     S(iDataset).yearVector = yearVectorData;
%     
%     % Calculate mean by year
%     datasetMeanByYear = NaN(numel(yearVectorData),1); 
%     for iYear = 1:numel(yearVectorData)
%         thisYearData = TS.sceneaverage(TS.year == yearVectorData(iYear));
%         datasetMeanByYear(iYear) = mean(thisYearData,'omitnan'); 
%     end
%     S(iDataset).meanByYear = datasetMeanByYear;
%     
%     % Calculate monthly mean by year
%     datasetMonthlyMeanByYear = NaN(12,numel(yearVectorData)); 
%     for iYear = 1:numel(yearVectorData)
%         for iMonth = 1:12
%             thisMonthAndYearSceneAverage = TS.sceneaverage(TS.month(:) == iMonth & TS.year == yearVectorData(iYear));
%             datasetMonthlyMeanByYear(iMonth,iYear) = mean(thisMonthAndYearSceneAverage,'omitnan'); 
%         end
%     end
%     S(iDataset).monthlyMeanByYear = datasetMonthlyMeanByYear;
% 
% end
% 

end % extractDriverDataFromProducts