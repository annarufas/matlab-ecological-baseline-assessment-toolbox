function plotPftsFromOceanColourProducts(AScmems,ASbicep)

% PLOTPFTSFROMOCEANCOLOURPRODUCTS Plot the mean relative distribution of
% PFTs in CMEMS ocean colour products.
%
%   INPUT: 
%       AScmems  - table with time-series data from CMEMS for our area of study
%       AScbicep - table with time-series data from BICEP for our area of study
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

% Products selected
cmemsProdLabel = {'1 km, regional, CMEMS multi',...
                  '7 km, regional, CMEMS reanalysis'};
cmemsProdName =  {'regional_cmemsmulti_r1km',...
                  'regional_cmemsreanalysis_r7km'};

bicepProdLabel = '9 km, global, BICEP';  
bicepProdName  = 'bicep_cphyto_r9km';  

% Initialise structure to hold data
cmemsPlkStruct = struct();
bicepPlkStruct = struct();

% Define variable names to extract
cmemsPlkVarNames = {'DIATO','DINO','NANO','PICO','CHL'};
bicepPlkVarNames = {'C_microphyto','C_nanophyto','C_picophyto','chl_a'};

%% Read in time-series data from CMEMS multi 1 km resolution

idxProd = find(strcmp({AScmems.ID}, 'obs_satell_reg_cmems_multi_1km_plk'));

% Find indices of plankton variable names from the dataset
prodVarNames = AScmems(idxProd).varNames;
idxsPlanktonVars = [];
for i = 1:numel(prodVarNames)
    for j = 1:numel(cmemsPlkVarNames)
        % Match indices of plankton variable names using strcmp (case-insensitive)
        if strcmpi(prodVarNames{i}, cmemsPlkVarNames{j})
            idxsPlanktonVars(end + 1) = i; % add the index to the array
        end
    end
end

% Extract data for each plankton variable
for i = 1:numel(cmemsPlkVarNames)
    cmemsPlkStruct.regional_cmemsmulti_r1km.dataset.(cmemsPlkVarNames{i}) =...
        AScmems(idxProd).dataset(:,:,:,idxsPlanktonVars(i));
end
cmemsPlkStruct.regional_cmemsmulti_r1km.time = AScmems(idxProd).time;
cmemsPlkStruct.regional_cmemsmulti_r1km.lat = AScmems(idxProd).lat;
cmemsPlkStruct.regional_cmemsmulti_r1km.lon = AScmems(idxProd).lon;

%% Read in time-series data from CMEMS biogeochemical reanalysis 
 
idxProd = find(strcmp({AScmems.ID}, 'mod_bgc_reg_diat'));
cmemsPlkStruct.regional_cmemsreanalysis_r7km.dataset.(cmemsPlkVarNames{1}) = AScmems(idxProd).dataset(:,:,:,:);
cmemsPlkStruct.regional_cmemsreanalysis_r7km.time = AScmems(idxProd).time;
cmemsPlkStruct.regional_cmemsreanalysis_r7km.lat = AScmems(idxProd).lat;
cmemsPlkStruct.regional_cmemsreanalysis_r7km.lon = AScmems(idxProd).lon;
cmemsPlkStruct.regional_cmemsreanalysis_r7km.depth = AScmems(idxProd).depth;

idxProd = find(strcmp({AScmems.ID}, 'mod_bgc_reg_dino'));
cmemsPlkStruct.regional_cmemsreanalysis_r7km.dataset.(cmemsPlkVarNames{2}) = AScmems(idxProd).dataset(:,:,:,:);
idxProd = find(strcmp({AScmems.ID}, 'mod_bgc_reg_nano'));
cmemsPlkStruct.regional_cmemsreanalysis_r7km.dataset.(cmemsPlkVarNames{3}) = AScmems(idxProd).dataset(:,:,:,:);
idxProd = find(strcmp({AScmems.ID}, 'mod_bgc_reg_pico'));
cmemsPlkStruct.regional_cmemsreanalysis_r7km.dataset.(cmemsPlkVarNames{4}) = AScmems(idxProd).dataset(:,:,:,:);
idxProd = find(strcmp({AScmems.ID}, 'mod_bgc_reg_chl'));
cmemsPlkStruct.regional_cmemsreanalysis_r7km.dataset.(cmemsPlkVarNames{5}) = AScmems(idxProd).dataset(:,:,:,:);

%% Read in time-series data from BICEP's Cphyto product

idxProd = find(strcmp({ASbicep.ID}, 'bicep_cphyto_9km'));

% Find indices of plankton variable names from the dataset
prodVarNames = ASbicep(idxProd).varNames;
idxsPlanktonVars = [];
for i = 1:numel(prodVarNames)
    for j = 1:numel(bicepPlkVarNames)
        % Match indices of plankton variable names using strcmp (case-sensitive)
        if strcmp(prodVarNames{i}, bicepPlkVarNames{j})
            idxsPlanktonVars(end + 1) = i; % add the index to the array
        end
    end
end

% Extract data for each plankton variable
for i = 1:numel(bicepPlkVarNames)
    bicepPlkStruct.bicep_cphyto_r9km.dataset.(bicepPlkVarNames{i}) =...
        ASbicep(idxProd).dataset(:,:,:,idxsPlanktonVars(i));
end
bicepPlkStruct.bicep_cphyto_r9km.time = ASbicep(idxProd).time;
bicepPlkStruct.bicep_cphyto_r9km.lat = ASbicep(idxProd).lat;
bicepPlkStruct.bicep_cphyto_r9km.lon = ASbicep(idxProd).lon;

%% Compute average data for our area of study

% CMEMS
for iProduct = 1:numel(cmemsProdName)
    thisProductName = cmemsProdName{iProduct};
    time = cmemsPlkStruct.(thisProductName).time;
    nTimeSteps = length(time);
    [cmemsPlkStruct] = calculateAndAddAverageDataToTheStructure(thisProductName,...
        cmemsPlkVarNames,cmemsPlkStruct,time,nTimeSteps); 
end 

% BICEP
thisProductName = bicepProdName;
time = bicepPlkStruct.(thisProductName).time;
nTimeSteps = length(time);
[bicepPlkStruct] = calculateAndAddAverageDataToTheStructure(thisProductName,...
    bicepPlkVarNames,bicepPlkStruct,time,nTimeSteps);

%% Plot the relative abundance of PFTs for CMEMS products

for iProduct = 1:numel(cmemsProdName)

    thisProductName = cmemsProdName{iProduct};
    
    % Extract variables
    plkDataset = cmemsPlkStruct.(thisProductName).dataset;
    chla = plkDataset.CHL_monthlyMean.sceneAverage;
    time = plkDataset.CHL_monthlyMean.time;
    diat = plkDataset.DIATO_monthlyMean.sceneAverage; % micro
    dino = plkDataset.DINO_monthlyMean.sceneAverage; % micro
    nano = plkDataset.NANO_monthlyMean.sceneAverage; % haptophytes + green algae
    pico = plkDataset.PICO_monthlyMean.sceneAverage; % prochlorophytes + prokaryotes
    % hapto = plkDataset.HAPTO_monthlyMean; % coccolithophores
    % green = plkDataset.GREEN_monthlyMean; % green algae 
    % prokar = plkDataset.PROKAR_monthlyMean; % prokaryotes
    % prochlo = plkDataset.PROCHLO_monthlyMean; % prochlorophytes

    % Compute total and relative abundances
    total = diat + dino + nano + pico;
    diatRelat = 100 * (diat./total);
    dinoRelat = 100 * (dino./total);
    nanoRelat = 100 * (nano./total);
    picoRelat = 100 * (pico./total);
    pftArray = [diatRelat, dinoRelat, nanoRelat, picoRelat];
    
    % Plotting
    figure()
    set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.35],'Color','w') 
    ax = axes('Position', [0.10 0.20 0.77 0.70]); % position in within the figure axes (DO NOT CHANGE)
    lgLabel = {'Diatoms','Dinoflagellates','Nanophytoplankton','Picophytoplankton'};
    setPftPlotFeatures(ax,time,pftArray,chla,cmemsProdLabel{iProduct},lgLabel)
    exportgraphics(gcf,fullfile('.','figures',strcat('PFT_CMEMS_',num2str(iProduct),'.png')),'Resolution',600)
    clear ax
    
end % iProduct

%% Plot the relative abundance of PFTs for BICEP's Cphyto

thisProductName = bicepProdName;

% Extract variables
plkDataset = bicepPlkStruct.(thisProductName).dataset;
chla = plkDataset.chl_a_monthlyMean.sceneAverage;
time = plkDataset.chl_a_monthlyMean.time;
micro = plkDataset.C_microphyto_monthlyMean.sceneAverage; 
nano = plkDataset.C_nanophyto_monthlyMean.sceneAverage; 
pico = plkDataset.C_picophyto_monthlyMean.sceneAverage; 

% Compute total and relative abundances
total = micro + nano + pico;
microRelat = 100 * (micro./total);
nanoRelat = 100 * (nano./total);
picoRelat = 100 * (pico./total);
pftArray = [microRelat, nanoRelat, picoRelat];

% Plotting
figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.35],'Color','w') 
ax = axes('Position', [0.10 0.20 0.77 0.70]); % position in within the figure axes (DO NOT CHANGE)
lgLabel = {'Microphytoplankton','Nanophytoplankton','Picophytoplankton'};
setPftPlotFeatures(ax,time,pftArray,chla,bicepProdLabel,lgLabel)
exportgraphics(gcf,fullfile('.','figures','PFT_BICEP.png'),'Resolution',600)
clear ax

% =========================================================================
%%
% -------------------------------------------------------------------------
% LOCAL FUNCTIONS
% -------------------------------------------------------------------------

% *************************************************************************

function [plkStruct] = calculateAndAddAverageDataToTheStructure(thisProductName,...
    plkVarNames,plkStruct,time,nTimeSteps)

for iPlkVar = 1:numel(plkVarNames)

    thisPlkVarName = plkVarNames{iPlkVar};
    plkVar = plkStruct.(thisProductName).dataset.(thisPlkVarName);

    if (isfield(plkStruct.(thisProductName),'depth') == 0) % if the dataset has no depth dependency
        [plkVarSceneAverage,~] = calcSceneAverage(plkVar,time);
    else % if the dataset has depth-dependency, compute the pondered depth average 
        nDepthLevels = numel(plkStruct.(thisProductName).depth);
        plkVarSceneAverage = calcScenePonderedDepthAverage(plkVar,time,nDepthLevels,nTimeSteps);
    end
    [plkVarSceneAverageYearlyMean,plkVarSceneAverageMonthlyMean] =... 
        calcSceneTimeAveragedData(plkVarSceneAverage,time);

    % Add fields to the structure array
    plkStruct.(thisProductName).dataset = setfield(plkStruct.(thisProductName).dataset,...
        strcat(thisPlkVarName,'_dailyMean'), plkVarSceneAverage);
    plkStruct.(thisProductName).dataset = setfield(plkStruct.(thisProductName).dataset,...
        strcat(thisPlkVarName,'_yearlyMean'), plkVarSceneAverageYearlyMean);
    plkStruct.(thisProductName).dataset = setfield(plkStruct.(thisProductName).dataset,...
        strcat(thisPlkVarName,'_monthlyMean'), plkVarSceneAverageMonthlyMean);   

end % iPlkVar

end % calculateAndAddAverageDataToTheStructure
 
% *************************************************************************

function setPftPlotFeatures(ax,time,pftArray,chla,prodLabel,lgLabel)

yyaxis left
h = bar(ax,time,pftArray,'stacked');
ylim([0 100])
yticks([0 20 40 60 80 100])
yticklabels({'0','20','40','60','80','100'})
ylabel('Relative abundance (%)','FontSize',12)
set(gca,'YColor','k')
hold on

yyaxis right
plot(ax,time,chla,'Color','k','LineWidth',1);
ylim([0 4])
ylabel('Chlorophyll a (mg m^{-3})','FontSize',12)
set(gca,'YColor','k')
hold off

% Determine start month for data trimming
iStartMonth = find(~isnan(pftArray(:,1)),1);
dmin = time(iStartMonth);
dmax = time(end);

xlim([dmin dmax])
startYear = year(dmin);
endYear = year(dmax);
xtickTimes = datetime(startYear, 1, 1):calyears(2):datetime(endYear, 12, 31);
xtickTimes.Format = 'yyyy'; % use 'yyyy' to display only the year
xtickLabels = string(xtickTimes);
xticks(xtickTimes);
xticklabels(xtickLabels);
xtickangle(45)
set(gca,'TickLength',[0.02, 0.02])

% Change the colors of each bar segment
nPfts = size(pftArray,2); % num. columns
coloursPft = colormap(brewermap(nPfts,'*set3'));
coloursPft = repelem(coloursPft,size(h,1),1); 
coloursPft = mat2cell(coloursPft,ones(size(coloursPft,1),1),3);
set(h,{'FaceColor'},coloursPft)

title(prodLabel,'FontSize',12)

lg = legend(lgLabel);
lg.Position(1) = 0.37; lg.Position(2) = -0.04;
lg.Orientation = 'horizontal';
lg.ItemTokenSize = [15,15];
lg.FontSize = 12; 
set(lg,'Box','off') 

set(ax,'FontSize',12)
box on
    
end % setPftPlotFeatures

% *************************************************************************

end % plotPftsFromOceanColourProducts