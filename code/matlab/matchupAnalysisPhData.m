function matchupAnalysisPhData(filenameCmemsMatchupsPh,PH)

% MATCHUPANALYSISPHDATA Perform matchups between in situ pH data and CMEMS
% biogeochemical reanalysis pH data. Regression statistics are calculated 
% and data is plotted for visual inspection.
%
%   INPUT:
%       filenameCmemsMatchupsPh - .csv file with CMEMS pH data matchups
%       PH                      - table with in situ pH data
%
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 29 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Read in matchup file 

cmemsMatchups = readtable(fullfile('.','data','raw','CMEMS_data',filenameCmemsMatchupsPh));

%% Create an empty table to store matchup statistics of regression

columnHeaders = {'N','slope','intercept','R2','r','RMSE','log_RMSE',...
    'ME','log_ME','MAE','log_MAE','MRPE','MAPE'};
variableTypes = {'double','double','double','double','double','double','double',...
    'double','double','double','double','double','double'};
rowNames = {'all_data','subset_data'};

T = table('Size',[numel(rowNames),numel(columnHeaders)],...
    'VariableNames',columnHeaders,'VariableTypes',variableTypes,'RowNames',rowNames);

%% Calculate regression statistics for two cases

% The data 
xvar = PH.pH;
yvar = cmemsMatchups.ph_mod_bgc_reg_ph;

% Case 1: remove concurrent entries where at least one of the datasets has NaN
insituPhAll = xvar(~isnan(xvar) & ~isnan(yvar));
cmemsPhAll = yvar(~isnan(xvar) & ~isnan(yvar));
[T{1,:}] = calculateMatchupRegressionStatistics(insituPhAll,cmemsPhAll);

% Case 2: remove concurrent entries where at least one of the datasets has NaN 
% as well as where in situ data is above 9
insituPhSubset = xvar(~isnan(xvar) & ~isnan(yvar) & xvar <= 9);
cmemsPhSubset = yvar(~isnan(xvar) & ~isnan(yvar) & xvar <= 9);
[T{2,:}] = calculateMatchupRegressionStatistics(insituPhSubset,cmemsPhSubset);

%% Save matchup statistics

writetable(T,fullfile('.','data','processed','statistics_matchups_ph.csv'),... 
    'Delimiter', ',', 'WriteRowNames', true);

%% Plot reanalysis data vs in situ data

figure()

%%%%%%%%%%%%%%%%%%%%%%% Main figure %%%%%%%%%%%%%%%%%%%%%%%
hfig = set(gcf,'Units','Normalized','Position',[0.01 0.05 0.30 0.35],'Color','w');

%%%%%%%%%%% Presets %%%%%%%%%%%
x = insituPhAll;
y = cmemsPhAll;
titleString = 'pH matchups (all)';
xmin = 7.4;
ymin = xmin;
xmax = 11.4;
ymax = xmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
% Plot 1:1 reference line
xref = xmin:0.1:xmax; 
yref = xref;
plot(xref,yref,'k:','LineWidth',1.5); hold on;

% Scatter plot
c = getScatterHeatPlotColormapIndices(x,y,80);
scatter(x,y,50,c,'.');
hold on
colormap(brewermap(1000,'*YlGnBu'))
cb = colorbar;
cb.Label.String = 'Density of points';
% set(cb,'FontSize',12);
caxis([1 40]); % density data points
box on

% Fit line
fitModel = fit(x,y,'power1'); hold on; % has to be a power so that it becomes a straight line in a log-log space
plot(x, fitModel(x),'r-','LineWidth', 2);

% Axis settings
xlabel("In situ (SSB and UKOA programs)",'FontSize',12)
ylabel("CMEMS biogeochemical reanalysis",'FontSize',12)

xlim([xmin xmax]);
ylim([ymin ymax]);   

% Add statistics
xt = max(xlim) - 0.95 * (max(xlim) - min(xlim));
yt = max(ylim) - 0.020 * (max(ylim) - min(ylim));
text(xt, yt,... % text position relative to axis
     {[cell2mat(strcat({'R^{2} = '},num2str(T.R2(1),'%.2f')))],...
     [cell2mat(strcat({'Slope = '},num2str(T.slope(1),'%.2f')))],...
     [cell2mat(strcat({'{RMSE} = '},num2str(T.RMSE(1),'%.2f')))],...
     [cell2mat(strcat({'\delta = '},num2str(T.ME(1),'%.2f')))],...
     [cell2mat(strcat({'N = '},num2str(T.N(1),'%.0f')))]},...
     'FontSize',10,'Horiz','left','Vert','top');

% Title
title(titleString,'FontSize',12);

%%%%%%%%%%%%%%%%%%%%%%% Inset figure %%%%%%%%%%%%%%%%%%%%%%%

% Create a new pair of axes inside current figure
axes('Position',[0.40 0.45 0.34 0.39])

%%%%%%%%%%% Presets %%%%%%%%%%%
x = insituPhSubset;
y = cmemsPhSubset;
titleString = 'pH matchups (subset)';
xmin = 7.4;
ymin = xmin;
xmax = 9;
ymax = xmax;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot 1:1 reference line
xref = xmin:0.1:xmax; 
yref = xref;
plot(xref,yref,'k:','LineWidth',1.5); hold on;

% Scatter plot
c = getScatterHeatPlotColormapIndices(x,y,80);
scatter(x,y,50,c,'.');
hold on
colormap(brewermap(1000,'*YlGnBu'))
% cb = colorbar;
% cb.Label.String = 'Density of points';
% set(cb,'FontSize',12);
caxis([1 40]); % density data points
box on

% Fit line
fitModel = fit(x,y,'power1'); hold on; % has to be a power so that it becomes a straight line in a log-log space
plot(x, fitModel(x),'r-','LineWidth', 2);

% Axis settings
% xlabel("In situ (SSB and UKOA programs)",'FontSize',12)
% ylabel("CMEMS biogeochemical reanalysis",'FontSize',12)

xlim([xmin xmax]);
ylim([ymin ymax]);   

% Add statistics
xt = max(xlim) - 0.95 * (max(xlim) - min(xlim));
yt = max(ylim) - 0.020 * (max(ylim) - min(ylim));
text(xt, yt,... % text position relative to axis
     {[cell2mat(strcat({'R^{2} = '},num2str(T.R2(2),'%.2f')))],...
     [cell2mat(strcat({'Slope = '},num2str(T.slope(2),'%.2f')))],...
     [cell2mat(strcat({'{RMSE} = '},num2str(T.RMSE(2),'%.2f')))],...
     [cell2mat(strcat({'\delta = '},num2str(T.ME(2),'%.2f')))],...
     [cell2mat(strcat({'N = '},num2str(T.N(2),'%.0f')))]},...
     'FontSize',8,'Horiz','left','Vert','top');

% Title
title(titleString,'FontSize',12);

exportgraphics(gcf,fullfile('.','figures','pH_matchups.png'),'Resolution',600)

end % matchupAnalysisPhData