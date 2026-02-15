function matchupAnalysisHplcData(filenameCmemsMatchupsHplc,filenameNasaMatchupsHplc,...
    filenameOccciMatchupsHplc,HPLC,AScmems,ASnasa,ASoccci,pathEnduranceShapefile,...
    pathAreaStudyShapefile)

% MATCHUPANALYSISHPLCDATA Perform matchups between chlorophyll concentration
% data from in situ HPLC measurements with L3 satellite observations from 
% various products downloaded from CMEMS, NASA and OCCCI. Regression 
% statistics are calculated and data is plotted for visual inspection, both
% regression plots and scene analysis plots.
%
%   INPUT:
%       filenameCmemsMatchupsHplc - .csv file with CMEMS chlorophyll data matchups
%       filenameNasaMatchupsHplc  - .csv file with NASA chlorophyll data matchups
%       filenameOccciMatchupsHplc - .csv file with OCCCI chlorophyll data matchups
%       HPLC                      - table with in situ HPLC data
%       AScmems                   - table with time-series data from CMEMS for our area of study
%       ASnasa                    - table with time-series data from NASA for our area of study
%       ASoccci                   - table with time-series data from OC-CCI for our area of study
%       pathEnduranceShapefile    - shapefile with the Endurance GCS site
%       pathAreaStudyShapefile    - shapefile with our area of study
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 1 May 2024 
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Read in time-series data for scenes

[chlaStruct.scene,prodLabel,prodName,nChlaProducts] = ...
    extractChlorophyllFromOceanColourProducts(AScmems,ASnasa,ASoccci);

%% Read in matchup data

cmemsMatchupTable = readtable(fullfile('.','data','raw','CMEMS_data',filenameCmemsMatchupsHplc));
chlaStruct.matchup.global_olci_r4km              = cmemsMatchupTable.obs_satell_glob_cmems_olci_4km_plk;
chlaStruct.matchup.regional_cmemsmulti_r1km      = cmemsMatchupTable.obs_satell_reg_cmems_multi_1km_plk;
chlaStruct.matchup.regional_olci_r300m           = cmemsMatchupTable.obs_satell_reg_cmems_olci_300m_plk; 
chlaStruct.matchup.regional_cmemsreanalysis_r7km = cmemsMatchupTable.mod_bgc_reg_chl;

nasaMatchupTable = readtable(fullfile('.','data','raw','NASA_data',filenameNasaMatchupsHplc));
chlaStruct.matchup.global_aquamodis_r4km = nasaMatchupTable.aquamodis_4km;
chlaStruct.matchup.global_viirssnpp_r4km = nasaMatchupTable.viirssnpp_4km;

occciMatchupTable = readtable(fullfile('.','data','raw','OCCCI_data',filenameOccciMatchupsHplc));
chlaStruct.matchup.global_occcimulti_r4km = occciMatchupTable.occci_4km_1day;
chlaStruct.matchup.global_occcimulti_r1km = occciMatchupTable.occci_1km_1day;

%% Create an empty table to store matchup statistics of regression

columnHeaders = {'N','slope','intercept','R2','r','RMSE','log_RMSE',...
    'ME','log_ME','MAE','log_MAE','MRPE','MAPE'};
variableTypes = {'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double'};
rowNames = {'4km-global-OLCI','4km-global-AquaMODIS','4km-global-VIIRSSNPP',...
    '4km-global-OCCCImulti','1km-global-OCCCImulti','1km-regional-CMEMSmulti',...
    '300m-regional-OLCI','7km-regional-CMEMSreanalysis'};
T = table('Size',[numel(rowNames),numel(columnHeaders)],...
    'VariableNames',columnHeaders,'VariableTypes',variableTypes,'RowNames',rowNames);

%% Calculate regression statistics 

for iProduct = 1:nChlaProducts
    
    thisProductName = prodName{iProduct};
    y = chlaStruct.matchup.(thisProductName);
    x = HPLC.TChlA_ug_L(:);
    
    % Remove concurrent entries where at least one of the datasets has NaN
    insituChlNonNaN = x(~isnan(x) & ~isnan(y));
    prodChlNonNaN = y(~isnan(x) & ~isnan(y));
    
    % Calculate regression statistics
    [T{iProduct,:}] = calculateMatchupRegressionStatistics(insituChlNonNaN,prodChlNonNaN);
    
end % iProduct

%% Save matchup statistics

writetable(T,fullfile('.','data','processed','statistics_matchups_hplc_l3.csv'),... 
    'Delimiter', ',', 'WriteRowNames', true);

%% Plot matchup subplots

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.60 0.40],'Color','w')
haxis = zeros(nChlaProducts,1);

for iProduct = 1:nChlaProducts

    haxis(iProduct) = subaxis(2,4,iProduct,'Spacing',0.015,'Padding',0.010,'Margin',0.10);
    ax(iProduct).pos = get(haxis(iProduct),'Position');
    if (iProduct <= 4)
        ax(iProduct).pos(2) = ax(iProduct).pos(2) + 0.06;
    elseif (iProduct > 4)
        ax(iProduct).pos(2) = ax(iProduct).pos(2) - 0.02;
    end
    set(haxis(iProduct),'Position',ax(iProduct).pos) 

    thisProductName = prodName{iProduct};
    y = chlaStruct.matchup.(thisProductName);
    x = HPLC.TChlA_ug_L(:);
    
    % Remove concurrent entries where at least one of the datasets has NaN
    insituChlNonNaN = x(~isnan(x) & ~isnan(y));
    prodChlNonNaN = y(~isnan(x) & ~isnan(y));
    
    % Plot 1:1 reference line
    xref = 10.^(-2:0.1:2); 
    yref = xref;
    plot(xref,yref,'k-','LineWidth',1.5); hold on;
    
    % Plot 2:1 reference line
    xref = 10.^(-2:0.1:2); 
    yref = 0.5.*xref;
    plot(xref,yref,'k:','LineWidth',1); hold on;
    
    % Plot 1:2 reference line
    yref = 10.^(-2:0.1:2); 
    xref = 0.5.*yref;
    plot(xref,yref,'k:','LineWidth',1); hold on;
    
    % Plot the fit line
%     xfit = linspace(min(insituChlNonNaN),max(insituChlNonNaN),100); %10.^(log(min(insituChlNonNaN)):0.1:log(max(insituChlNonNaN)));
%     yfit = polyval(regressCoeff, xfit);
%     plot(xfit,yfit,'r-','LineWidth',1.5);
 
    % Scatter plot
    c = getScatterHeatPlotColormapIndices(insituChlNonNaN,prodChlNonNaN,80);
    scatter(haxis(iProduct),insituChlNonNaN,prodChlNonNaN,50,c,'.');
    hold on
    colormap(brewermap(1000,'*YlGnBu'))
    if (iProduct == nChlaProducts)
        cb = colorbar;
        cb.Label.String = 'Density of points';
        set(cb,'FontSize',11);
        set(cb, 'Position', [0.91 0.087 0.017 0.865]);
        caxis([0 50]);
    end
    box on
    set(gca,'Xscale','log','Yscale','log') % we need log scale to better visualise separation between scatter dots

    % Fit line
    fitModel = fit(insituChlNonNaN,prodChlNonNaN,'power1'); hold on; % has to be a power so that it becomes a straight line in a log-log space
    plot(insituChlNonNaN, fitModel(insituChlNonNaN),'r-','LineWidth', 2);

    % Axis settings
    if (iProduct == 6)
        xl = xlabel("In situ chl (mg m^{-3})",'FontSize',12);
        xl.Position(1) = xl.Position(1) + 250;
        xl.Position(2) = xl.Position(2) + 0.01;
    end
    if (iProduct == 1)
        yl = ylabel("Satellite chl (mg m^{-3})",'FontSize',12);
        yl.Position(1) = yl.Position(1) + 0.005;
        yl.Position(2) = yl.Position(2) - 0.70;
    end

    xlim([0.03 110]);
    ylim([0.03 110])   
    xticks([1e-1 1 10 100]);
    yticks([1e-1 1 10 100])
    set(gca, 'XTickLabels', {'0.1','1','10','100'});
    set(gca, 'YTickLabels', {'0.1','1','10','100'});
    
    % Add statistics
    xt = min(xlim)+0.50*min(xlim); 
    yt = max(ylim)-0.02*max(ylim);
    
    text(xt, yt,... % text position relative to axis
         {[cell2mat(strcat({'R^{2} = '},num2str(T.R2(iProduct),'%.2f')))],...
         [cell2mat(strcat({'Slope = '},num2str(T.slope(iProduct),'%.2f')))],...
         [cell2mat(strcat({'{RMSE} = '},num2str(T.RMSE(iProduct),'%.2f')))],...
         [cell2mat(strcat({'\delta = '},num2str(T.ME(iProduct),'%.2f')))],...
         [cell2mat(strcat({'N = '},num2str(T.N(iProduct),'%.0f')))]},...
         'FontSize',7.5,'Horiz','left','Vert','top');

    % Title
    tl = title(prodLabel{iProduct},'FontSize',11);
    tl.Position(2) = tl.Position(2) + 20;
    
end

exportgraphics(gcf,fullfile('.','figures','HPLC_matchups.png'),'Resolution',600)
clear ax

%% Plot scenes of chl, averaged over a week

areaStudy = m_shaperead(pathAreaStudyShapefile); % lat/lon coordinates
endurance = m_shaperead(pathEnduranceShapefile);

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.60 0.40],'Color','w')
haxis = zeros(nChlaProducts,1);

for iProduct = 1:nChlaProducts

    haxis(iProduct) = subaxis(2,4,iProduct,'Spacing',0.015,'Padding',0.010,'Margin',0.1);
    ax(iProduct).pos = get(haxis(iProduct),'Position');
    if (iProduct <= 4)
        ax(iProduct).pos(2) = ax(iProduct).pos(2) + 0.06;
    elseif (iProduct > 4)
        ax(iProduct).pos(2) = ax(iProduct).pos(2) - 0.02;
    end
    set(haxis(iProduct),'Position',ax(iProduct).pos) 

    thisProductName = prodName{iProduct};
    chla = chlaStruct.scene.(thisProductName).dataset;
    lon = chlaStruct.scene.(thisProductName).lon;
    lat = chlaStruct.scene.(thisProductName).lat;
    time = chlaStruct.scene.(thisProductName).time;
    
    % Pick up 1st week in May 2020 and do a 7-days average of the data
    [~,closestIndexStart] = min(abs(time - datetime('2022-04-01 00:00:00')));
    [~,closestIndexEnd] = min(abs(time - datetime('2022-04-07 00:00:00')));
    chlaMeanWeek = mean(chla(:,:,closestIndexStart:closestIndexEnd),3,'omitnan');
    
    cmin = log10(0.5); %log10(min(chlaMeanWeek,[],'all','omitnan'));
    cmax = log10(10); %log10(max(chlaMeanWeek,[],'all','omitnan'));

    m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
    m_pcolor(lon,lat,log10(chlaMeanWeek))
    shading flat
    colormap(brewermap(100,'*YlGnBu'))
    caxis([cmin cmax])

    m_gshhs_i('color','k');
    m_grid('linewi',1,'tickdir','out','FontSize',9);
    m_line(endurance.ncst{1}(:,1),endurance.ncst{1}(:,2),'linewi',1,'color','r');

    % Axis settings
    if (iProduct == 6)
        xl = xlabel("Longitude (ºE)",'FontSize',12);
        xl.Position(1) = xl.Position(1) + 0.010;
        xl.Position(2) = xl.Position(2);
    end
    if (iProduct == 1)
        yl = ylabel("Latitude (ºN)",'FontSize',12);
        yl.Position(1) = yl.Position(1);
        yl.Position(2) = yl.Position(2) - 0.011;
    end
    
    tl = title(prodLabel{iProduct},'FontSize',11);
    tl.Position(2) = tl.Position(2) + 0.0001;

    if (iProduct == nChlaProducts)
        cb = colorbar('Location','eastoutside');
        cb.Label.String = 'Chlorophyll a (mg m^{-3})';
        set(cb,'ytick',log10([0.5 1 2 3 5 10 20 30]),'yticklabel',[0.5 1 2 3 5 10 20 30],...
            'tickdir','out','FontSize',10);
        set(cb, 'Position', [0.91 0.087 0.017 0.865]);
    end

end
    
exportgraphics(gcf,fullfile('.','figures','scenes_chlorophyll_avg7days.png'),'Resolution',600)
clear ax

end % matchupAnalysisHplcData
   