
clc
clear all

fullPathMainDir         = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/';
fullPathDataCMEMSdir    = strcat(fullPathMainDir,'CMEMS_data/data_timeseries_new/');
fullPathPlotsDir        = strcat(fullPathMainDir,'plots/ecosystem/');
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

%% Load datasets

% Clorophyll a
load(strcat(fullPathMainDir,'continuousSatelliteChlaFiveDay4km.mat'),...
    'continuousSatelliteChlaFiveDay4km','timeVectorCompleteFiveDay4km',...
    'latSatellite','lonSatellite')

% POC, Cphyto and NPP
load(strcat(fullPathMainDir,'carbonStocks.mat'))

% CMEMS biogeochemical and physical reanalysis products
load(strcat(fullPathDataCMEMSdir,'OC_cmems_bbox50km.mat'))

% Suspended aprticulate matter
load(strcat(fullPathMainDir,'spmSilva.mat'))

%% Create period averages for CMEMS products

listIdCmemsProd = {'mod_bgc_reg_kd',...
                   'mod_bgc_reg_no3',...
                   'mod_bgc_reg_o2',...
                   'mod_bgc_reg_ph',...
                   'mod_phy_reg_mld',...
                   'mod_phy_reg_sal',...
                   'mod_phy_reg_temp',...
                   'mod_phy_reg_ssh',...
                   'mod_phy_reg_velo'};
               
nCmemsProdFinalList = numel(listIdCmemsProd) + 1;

% Do some processing first

iProcessedProduct = 1;
             
for iProd = 1:numel(listIdCmemsProd)
    
    idProd = listIdCmemsProd{iProd};
    posIdProd = find(strcmp({OC_cmems.ID}, idProd));
    cmemsTime = OC_cmems(posIdProd).time;
    cmemsLat = OC_cmems(posIdProd).lat;
    cmemsLon = OC_cmems(posIdProd).lon;
    nTimeSteps = length(cmemsTime);
    nLatPixelsCmems = length(cmemsLat);
    nLonPixelsCmems = length(cmemsLon);
    latVectorAreaStudy = linspace(minLatAreaStudy,maxLatAreaStudy,nLatPixelsCmems)';
    lonVectorAreaStudy = linspace(minLonAreaStudy,maxLonAreaStudy,nLonPixelsCmems)';
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
    
    % These products don't need any processing as are only provided for one depth   
    elseif (strcmp(OC_cmems(posIdProd).ID,'mod_phy_reg_mld') ||...
            strcmp(OC_cmems(posIdProd).ID,'mod_phy_reg_ssh')) 
    
        cmemsProcessedProd(:,:,:,iProcessedProduct) = cmemsProd;
        iProcessedProduct = iProcessedProduct + 1;
        
    % Products with depth levels    
    elseif (strcmp(OC_cmems(posIdProd).ID,'mod_bgc_reg_no3') ||...
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
 
    % Split into eastward and northward velocities    
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
    [qXm,qYm,qTm] = ndgrid(latVectorAreaStudy,lonVectorAreaStudy,datenum(cmemsTime));
    % Actual grid
    [Xm,Ym,Tm] = ndgrid(cmemsLat,cmemsLon,datenum(cmemsTime));
    
    if (iProd < numel(listIdCmemsProd))
        % Interpolant
        Fm = griddedInterpolant(Xm,Ym,Tm,cmemsProcessedProd(:,:,:,iProd));
        % Regridded product
        cmemsRegriddedProductToAreaStudy = Fm(qXm,qYm,qTm); 
        cmemsFilledProduct(:,:,:,iProd) = Racault2014interpolation(cmemsRegriddedProductToAreaStudy);
    elseif (iProd == numel(listIdCmemsProd)) % to account for the two velocity products
        Fm = griddedInterpolant(Xm,Ym,Tm,cmemsProcessedProd(:,:,:,iProd));
        cmemsRegriddedProductToAreaStudy = Fm(qXm,qYm,qTm); 
        cmemsFilledProduct(:,:,:,iProd) = Racault2014interpolation(cmemsRegriddedProductToAreaStudy);
        Fm = griddedInterpolant(Xm,Ym,Tm,cmemsProcessedProd(:,:,:,iProd+1));
        cmemsRegriddedProductToAreaStudy = Fm(qXm,qYm,qTm); 
        cmemsFilledProduct(:,:,:,iProd+1) = Racault2014interpolation(cmemsRegriddedProductToAreaStudy);
    end    
     
end

% Create period averages 

periodSceneCmems = NaN(12,nCmemsProdFinalList);

for iProd = 1:nCmemsProdFinalList
    
    dataArray = cmemsFilledProduct(:,:,:,iProd);
    timeVector = cmemsTime;

    % Average per scene
    dataArraySceneAveraged = NaN(length(timeVector),1); 
    for iTimeStep = 1:length(timeVector)
        thisScene = reshape(dataArray(:,:,iTimeStep),[],1);
        dataArraySceneAveraged(iTimeStep) = mean(thisScene,'omitnan');
    end

    % Use a table to arrange the scene data and average per month and year
    TS = table(dataArraySceneAveraged,timeVector,'VariableNames',{'data','date'});
    TS.month = month(TS.date);
    TS.year = year(TS.date);
    yearVectorData = min(TS.year):1:max(TS.year);

    dataMonthlyMean = NaN(12,numel(yearVectorData)); 
    for iYear = 1:numel(yearVectorData)
        for iMonth = 1:12
            thisMonthAndYear = TS.data(TS.month(:) == iMonth & TS.year == yearVectorData(iYear));
            dataMonthlyMean(iMonth,iYear) = mean(thisMonthAndYear,'omitnan'); 
        end
    end

    % Average per month in all years
    periodSceneCmems(:,iProd) = mean(dataMonthlyMean, 2, 'omitnan');

end

%% Create period averages for C stock products

% BICEP-Cphyto, BICEP-NPP and BICEP-POC goes 15 Jan 1998 to 15 Dec 2020
% CMEMS-Cphyto goes 1 Jan 1998 to 31 Dec 2022

periodSceneStock = NaN(12,4); % Mg
periodSceneConcentrationPerUnitArea = NaN(12,4); % g C m-2

for iVar = 1:4
    if (iVar == 1)
        posId = find(strcmp({carbonStocks.ID}, 'bicep_poc_4km'));
    elseif (iVar == 2)    
        posId = find(strcmp({carbonStocks.ID}, 'bicep_cphyto_9km'));
    elseif (iVar == 3)
        posId = find(strcmp({carbonStocks.ID}, 'cmems_cphyto'));
    elseif (iVar == 4)
        posId = find(strcmp({carbonStocks.ID}, 'bicep_npp_9km'));
    end

    % These arrays have size 12 x nYears
    sceneStock = carbonStocks(posId).sceneStockDistribByMonthAndYear_monthlyMean; % Mg
    sceneConcentrationPerUnitArea = carbonStocks(posId).sceneConcentrationPerUnitAreaDistribByMonthAndYear_monthlyMean; % g C m-2

    % Average per month in all years
    periodSceneStock(:,iVar) = mean(sceneStock, 2, 'omitnan');
    periodSceneConcentrationPerUnitArea(:,iVar) = mean(sceneConcentrationPerUnitArea, 2, 'omitnan'); 
end

%% Create period averages for chla product

% Chla goes 1 Jan 1998 to 27 Dec 2023

periodSceneChla = NaN(12,1);
dataArray = continuousSatelliteChlaFiveDay4km;
timeVector = timeVectorCompleteFiveDay4km;

% Average per scene
dataArraySceneAveraged = NaN(length(timeVector),1); 
for iTimeStep = 1:length(timeVector)
    thisScene = reshape(dataArray(:,:,iTimeStep),[],1);
    dataArraySceneAveraged(iTimeStep) = mean(thisScene,'omitnan');
end

% Use a table to arrange the scene data and average per month and year
TS = table(dataArraySceneAveraged,timeVector,'VariableNames',{'data','date'});
TS.month = month(TS.date);
TS.year = year(TS.date);
yearVectorData = min(TS.year):1:max(TS.year);

dataMonthlyMean = NaN(12,numel(yearVectorData)); 
for iYear = 1:numel(yearVectorData)
    for iMonth = 1:12
        thisMonthAndYear = TS.data(TS.month(:) == iMonth & TS.year == yearVectorData(iYear));
        dataMonthlyMean(iMonth,iYear) = mean(thisMonthAndYear,'omitnan'); 
    end
end

% Average per month in all years
periodSceneChla(:) = mean(dataMonthlyMean, 2, 'omitnan');

%% Regrid SPM and create scene average SPM

qLat = linspace(minLatAreaStudy,maxLatAreaStudy,numel(latSpm))';
qLon = linspace(minLonAreaStudy,maxLonAreaStudy,numel(lonSpm))';
    
% % Regrid to lat x lon
% % New grid (query points for interpolation)
% [qX,qY,qT] = ndgrid(qLat,qLon,(1:12)');
% % Actual grid
% [X,Y,T] = ndgrid(latSpm,lonSpm,(1:12)');
% % Interpolant
% F = griddedInterpolant(X,Y,T,spmMonthlyClimatologyMean);
% % Regridded product
% spmRegriddedToAreaStudy = F(qX,qY,qT); % 1e-3 kg m-3

spmRegriddedToAreaStudy = spmMonthlyClimatologyMean;
        
% Create scene average    
spmSceneAveraged = NaN(12,1); 
for iMonth = 1:12
    thisScene = reshape(spmRegriddedToAreaStudy(:,:,iMonth),[],1);
    spmSceneAveraged(iMonth) = mean(thisScene,'omitnan');
end   

periodSceneSpm = spmSceneAveraged;

%% Plot

% Start in July for better visualisation of the MLD dip
% months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
months = {'Jul','Aug','Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun'};

colourCarbonProduct = [1 1 1; 0.7 0.7 0.7; 0 0 0];
yBreakValueCarbon = 6; % 10

colourCmemsProduct = brewermap(4,'*Paired');

shiftRight = 0.08;
heigtSubplot = 0.15;

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.50],'Color','w') 

for iSubplot = 1:4
    
    % Chla product
    if (iSubplot == 1) 
        
        ax1 = axes('Position', [shiftRight 0.825 0.80 0.12]); 
        chlaVectorShifted = [periodSceneChla(6:12);periodSceneChla(1:5)];
        plot(ax1,(1:12),chlaVectorShifted,'o-','Color','k','LineWidth',1,...
            'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor','r'); hold on;
        ylim([0 2])
        xticks(1:12);
        xlim([0.5 12.5])
        xticklabels([])
        ytickformat('%.1f');
        ax1.XGrid = 'on';
        ax1.YGrid = 'off'; 
        hold on

    elseif (iSubplot == 2) % C products
          
        for iVar = 1:3

            if (iVar == 1) % POC
                carbonVectorShifted = [periodSceneConcentrationPerUnitArea(6:12,1);...
                    periodSceneConcentrationPerUnitArea(1:5,1)];
                colourThis = colourCarbonProduct(1,:);
            elseif (iVar == 2) % Cphyto
                carbonVectorShifted = [periodSceneConcentrationPerUnitArea(6:12,2);...
                    periodSceneConcentrationPerUnitArea(1:5,2)];
                colourThis = colourCarbonProduct(2,:);
            elseif (iVar == 3) % NPP
                carbonVectorShifted = [periodSceneConcentrationPerUnitArea(6:12,4);...
                    periodSceneConcentrationPerUnitArea(1:5,4)];
                colourThis = colourCarbonProduct(3,:);
            end

            % Create axes for the bottom part of the plot
            if (iVar == 1)    
                ax2 = axes('Position',[shiftRight 0.42 0.80 0.30]); 
                ax2.XGrid = 'on';
                ax2.YGrid = 'off';
            end
            hold on

            plot(ax2,(1:12),carbonVectorShifted,'-o','Color','k','LineWidth',1,...
                'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colourThis); hold on;
            ylim([0 yBreakValueCarbon])
            if (iVar == 1)
                yticks([0:2:yBreakValueCarbon])
            end
            xticks(1:12);
            xlim([0.5 12.5])
            xticklabels([])
            ax1.TickDir = 'in';
            ax1.YColor = 'k'; % Set the color of the left y-axis to black
            box on
            hold on

            % Create axes for the upper part of the plot
            if (iVar == 1)   
                ax3 = axes('Position', [shiftRight 0.725 0.80 0.08]); 
                ax3.XAxis.Visible = 'on';
                ax3.YScale = 'log';
                ax3.XGrid = 'on';
                ax3.YGrid = 'off';
            end
            hold on

            plot(ax3,(1:12),carbonVectorShifted,'-o','Color','k','LineWidth',1,...
                'MarkerSize',8,'MarkerEdgeColor','k','MarkerFaceColor',colourThis); hold on;
            ylim([yBreakValueCarbon 29])
            if (iVar == 1)   
                yticks([10 20 30])
                yticklabels({'10','20','30'})
            end
            xticks(1:12);
            xlim([0.5 12.5])
            xticklabels([])
            box on
            hold on
            
            ax2.XGrid = 'on';
            ax2.YGrid = 'off';
            ax3.XGrid = 'on';
            ax3.YGrid = 'off'; 
            ax3.TickDir = 'in';
                
        end % iVar
        hold on
    
    % Zeu and SPM
    elseif (iSubplot == 3) 

        ax4 = axes('Position', [shiftRight 0.22 0.80 heigtSubplot]);
        cmemsZeuShifted = [periodSceneCmems(6:12,1);periodSceneCmems(1:5,1)];
        spmShifted = [periodSceneSpm(6:12);periodSceneSpm(1:5)];
        
        yyaxis left
        plot(ax4,(1:12),cmemsZeuShifted,'-','Color',colourCmemsProduct(1,:),...
            'LineWidth',2); hold on;
        ylim([25 43])
        xticks(1:12);
        xlim([0.5 12.5])
        xticklabels([])
%         ylabel('Zeu (m)')
        set(ax4, 'YDir', 'reverse');
        ax4.TickDir = 'in';
        ax4.YColor = 'k'; % Set the color of the left y-axis to black
        ax4.XGrid = 'on';
        ax4.YGrid = 'off'; 

        yyaxis right
        plot(ax4,(1:12),spmShifted,'-','Color',colourCmemsProduct(2,:),...
            'LineWidth',2); hold on;
        ylim([0.6 2.6])
        ytickformat('%.1f');
        xticks(1:12);
        xlim([0.5 12.5])
        xticklabels([])
        ax4.TickDir = 'in';
%         ylabel('SPM')
        ax4.YColor = 'k'; 
        ax4.XGrid = 'on';
        ax4.YGrid = 'off'; 

        box on
        
    % MLD and temperature    
    elseif (iSubplot == 4) 
        
        ax5 = axes('Position', [shiftRight 0.05 0.80 heigtSubplot]);
        cmemsMldShifted = [periodSceneCmems(6:12,5);periodSceneCmems(1:5,5)];
        cmemsTempShifted = [periodSceneCmems(6:12,7);periodSceneCmems(1:5,7)]; 

        yyaxis left
        plot(ax5,(1:12),cmemsMldShifted,'-','Color',colourCmemsProduct(3,:),...
            'LineWidth',2); hold on;
        ylim([15 45])
%         ylabel('MLD (m)')
        xticks(1:12);
        xlim([0.5 12.5])
        xticklabels([])
        ax5.TickDir = 'in';
        set(ax5, 'YDir', 'reverse');
        ax5.XGrid = 'on';
        ax5.YGrid = 'off'; 
        ax5.YColor = 'k'; % Set the color of the left y-axis to black
        
        yyaxis right
        plot(ax5,(1:12),cmemsTempShifted,'-','Color',colourCmemsProduct(4,:),...
            'LineWidth',2); hold on;
        ylim([5 15])
%         ylabel('Temperature (ºC)')
        xticks(1:12);
        xlim([0.5 12.5])
        xticklabels(months)
        ax5.TickDir = 'in';
        ax5.XAxis.FontSize = 10.5;
        ax5.XGrid = 'on';
        ax5.YGrid = 'off'; 
        ax5.YColor = 'k'; % Set the color of the left y-axis to black

        box on
        
end % iSubplot

exportgraphics(gcf,fullfile(fullPathPlotsDir,...
    'plot_ocean_colour_vars_and_controls.png'),'Resolution',600)
