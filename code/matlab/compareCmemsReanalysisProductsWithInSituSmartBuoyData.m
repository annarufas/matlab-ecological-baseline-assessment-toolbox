% function compareCmemsReanalaysisProductsWithInSituSmartBuoyData(SB,SBcmems)

% PLOTSMARTBUOYDATACEFASVSCMEMS Create various plots to visualise the SmartBuoy dataset.
%
%   INPUT:
%       SB - SmartBuoy dataset from CEFAS
%       cmemsSBall - 
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 23 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Preamble

% CEFAS information
cefasVars = {'CHLOR','FTU','TOXN','SILICA','O2CONC','TEMP','SAL'};
cefasVarLabels = {'Chl a (\mug/L)','Turbidity (FTU)','Nitrate (\mumol/L)',...
    'Silica (\mumol/L)','Oxygen (mg/L)','Temperature (ºC)','Salinity (PSU)'}; % 'FTU' = 'Formazin Turbidity Units'
nCefasVars = length(cefasVars);

cefasBuoyNames = categories(SB.buoy);
nCefasBuoys = length(cefasBuoyNames);

cefasDepthCategories = unique(SB.depth);
nCefasDepthCategories = numel(cefasDepthCategories);
% DOWSING: 1 m and 15 m 
% NTHDOGGER: 25 m and 31 m 
% OYSTER: 35 m
% OYSTERML: 45 m

cefasStartDate = min(SB.dateTime);
cefasEndDate = max(SB.dateTime);

% CMEMS information
idCmemsProds = {'mod_bgc_reg_chl','mod_bgc_reg_kd','mod_bgc_reg_no3',...
    'mod_bgc_reg_o2','mod_phy_reg_temp','mod_phy_reg_sal'};
cmemsVars = {'chla','kd','no3','o2','temp','sal'};
cmemsBuoyNames = fieldnames(SBcmems);

% Define the buoy depth ranges
depthRanges = struct('Dowsing', [1, 15], 'NorthDogger', [25, 31], 'OysterGround', [35, 45]);

extractCmems = struct();

% Loop over buoy names
for iBuoy = 1:numel(cmemsBuoyNames)
    
    thisBuoyName = cmemsBuoyNames{iBuoy};
    thisBuoyStruct = SBcmems.(sprintf('%s', thisBuoyName));
%     fprintf('Buoy %s\n',thisBuoyName)
    
    % Get depth range for the current buoy
    buoyDepths = depthRanges.(thisBuoyName);
    
    % Loop over depth levels
    for iDepth = 1:2
        targetDepth = buoyDepths(iDepth);
        if iDepth == 1
            depthLabel = 'shallow';
        else
            depthLabel = 'deep';
        end
%         fprintf('Depth %d\n',iDepth)
    
        % Loop over products
        for iProd = 1:numel(idCmemsProds)
            prodVar = cmemsVars{iProd};
            thisProdId = idCmemsProds{iProd};
            thisProdIndex = find(strcmp({thisBuoyStruct.ID}, thisProdId));
%             fprintf('Product %s\n',thisProdId)
            
            % Find the closest start and end time indices
            timeVector = thisBuoyStruct(thisProdIndex).time;
            [~, closestStartTimeIndex] = min(abs(timeVector - cefasStartDate));
            [~, closestEndTimeIndex] = min(abs(timeVector - cefasEndDate));
    
            % Find the closest depth index
            depthVector = thisBuoyStruct(thisProdIndex).depth;
            [~, closestDepthIndex] = min(abs(depthVector - targetDepth));
            disp('Target')
            disp(targetDepth)
            disp(closestDepthIndex)
            disp(depthVector(closestDepthIndex))

            % Extract the dataset and time
            data = squeeze(thisBuoyStruct(thisProdIndex).dataset(:,:,:,closestDepthIndex));
            timeRange = timeVector(closestStartTimeIndex:closestEndTimeIndex);
            
            % Store the extracted data and time in the extractCmems structure
            extractCmems.(thisBuoyName).(prodVar).(depthLabel).var = data;
            extractCmems.(thisBuoyName).(prodVar).(depthLabel).time = timeRange;
        end
    end
end

% We will add pH

% Uncertainty???
% Need to calculate turbidity from kd

% Well, Kd490 is not a measure of turbidity, it's a measure of diffuse 
% attenuation. Turbidity is effectively a measure of scattering. While 
% scattering is part of attenuation, it is not everything, as attenuation 
% also includes absorption. So there is a relationship between turbidity 
% and Kd490, but it is not necessarily a simple one (although it *could* 
% be, if scattering is the dominant contributor to attenuation).
% 1/Kd is generally considered the 1st optical depth (roughly the 37% light 
% level), so you might be able to use that...

% Also, add the pH data from CEFAS and SPM. Create scripts to analyse those
% data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.40 0.85],'Color','w') 
% haxis = zeros(nCefasVars,1);
% 
% for iTracer = 1:nCefasVars
%     
%     haxis(iTracer) = subaxis(nCefasVars,1,iTracer,'Spacing',0.01,'Padding',0.01,'Margin', 0.10);
%     ax(iTracer).pos = get(haxis(iTracer),'Position');
%     ax(iTracer).pos(1) = ax(iTracer).pos(1)-0.045; % move all subplots a little bit to the left
%     set(haxis(iTracer),'Position',ax(iTracer).pos) 
%     
%     for iBuoy = 1:nCefasBuoys
%     
%         thisBuoyTracer = find(SB.variable == cefasVars(iTracer) & SB.buoy == cefasBuoyNames(iBuoy));
% %     iBuoySil = find(SB.variable == "SILICA" & SB.buoy == buoyCategories(iBuoy));
% %     iBuoyChl = find(SB.variable == "CHLOR" & SB.buoy == buoyCategories(iBuoy));
% %     iBuoyTem = find(SB.variable == "TEMP" & SB.buoy == buoyCategories(iBuoy));
% %     iBuoyOxy = find(SB.variable == "O2CONC" & SB.buoy == buoyCategories(iBuoy));
% %     iBuoySal = find(SB.variable == "SAL" & SB.buoy == buoyCategories(iBuoy));
% 
%         scatter(haxis(iTracer),SB.dateTime(thisBuoyTracer),SB.value(thisBuoyTracer),...
%             6,coloursBuoys(iBuoy,:),'filled'); hold on;
%     end
%     hold off;
%     
%     xlim([startDate endDate])    
% %     xtickangle(45)
% 
%     if (iTracer == 7)
%         set(gca,'YScale','log')
%         ylim([0 1000])
%         yticks([1e-1 1 10 1e2 ])
%         yticklabels({'0.1','1','10','100'})
%     end
%     ylh = ylabel(cefasVarLabels(iTracer),'FontSize',11);
%     ylh.Position(1) = -1300;
% 
%     if (iTracer == 1)
%         title('CEFAS SmartBuoy dataset','FontSize',12)
%     end
%         
%     if (iTracer == nCefasVars)
%         xlabel('Time (date)','FontSize',11);
%         lg = legend(labelLegend(:));  
%         lg.Orientation = 'vertical';
%         lg.Position(1) = 0.84; lg.Position(2) = 0.82;
%         lg.ItemTokenSize = [15,1];
%         lg.FontSize = 11;
%         set(lg,'Box','off')
%     end        
% 
%     box on
% 
% end 
% 
% hold off
% 
% exportgraphics(gcf,fullfile('.','figures','SmartBuoy_CEFAS_timeseries_all.png'),'Resolution',600)
% 
% clear ax
% 
% %% Inspect chlorophyll closer
% 
% startDate = datetime('2006-01-01');
% endDate = datetime('2013-09-15');
% 
% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.75 0.30],'Color','w') 
% for iBuoy = 1:nCefasBuoys
%     thisBuoyChla = find(SB.variable == "CHLOR" & SB.buoy == cefasBuoyNames(iBuoy));
%     x = SB.dateTime(thisBuoyChla);
%     y = SB.value(thisBuoyChla);
%     if (isempty(y))
%         plot(NaN,NaN,'-o','Color',coloursBuoys(iBuoy,:),'LineWidth',1); hold on;
%     else
%         plot(x,y,'-o','Color',coloursBuoys(iBuoy,:),'LineWidth',1); hold on; %'LineWidth',1.5,coloursBuoys(iBuoy,:)
%     end
% %     plot(SB.dateTime(thisBuoyChla),SB.value(thisBuoyChla),...
% %         6,coloursBuoys(iBuoy,:),'filled'); hold on;
%     if (iBuoy == nCefasBuoys)
%         lg = legend(labelLegend(:));  
%         lg.Orientation = 'vertical';
%         lg.Position(1) = 0.84; lg.Position(2) = 0.72;
%         lg.ItemTokenSize = [15,1];
%         lg.FontSize = 11;
%         set(lg,'Box','off')
%     end
% end
% hold off;
% 
% xlim([startDate endDate])
% ylim([0 1])
% grid on;
% ax = gca;
% newPosition = ax.Position;
% newPosition(1) = newPosition(1) - 0.08; % Adjust the value based on your needs
% newPosition(2) = newPosition(2) + 0.02; 
% ax.Position = newPosition;
% 
% ylabel('Chla (mg m^{-3}','FontSize',11)
% xlabel('Time (date)','FontSize',11);
% 
% % lg = legend(labelLegend(:));  
% % lg.Orientation = 'vertical';
% % lg.Position(1) = 0.84; lg.Position(2) = 0.72;
% % lg.ItemTokenSize = [15,1];
% % lg.FontSize = 11;
% % set(lg,'Box','off')
% box on
% 
% exportgraphics(gcf,fullfile('.','figures','SmartBuoy_CEFAS_timeseries_chla.png'),'Resolution',600) 
% clear ax
% 
% %% Histogram with all the data
% 
% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.44 0.29],'Color','w') 
% 
% % Main plot
% histogram(SB.value(SB.variable == "CHLOR"),50,'BinLimits',[0,20],'FaceColor', [0.5, 0.5, 0.5])
% ylim([1 6e4])
% set(gca,'YScale','log');
% grid on;
% mainAxes = gca;
% mainAxes.XGrid = 'on';  
% mainAxes.YGrid = 'on';  
% mainAxes.MinorGridAlpha = 0;
% % yticks([1, 1e1, 1e2, 1e3]);
% % yticklabels({'1','10','10^{2}','10^{3}'});
% xlabel('Chlorophyll a (mg m^{-3})');
% ylabel('Frequency');
% title('Distribution of SmartBuoy observations from CEFAS');
% 
% % Inset plot
% insetAxes = axes('Position', [0.52, 0.43, 0.36, 0.45]); % define the position of the inset plot
% histogram(insetAxes,SB.value(SB.variable == "CHLOR"),19,'BinLimits',[1e-1,2],'FaceColor', [0.5, 0.5, 0.5])
% ylim([1 4e4])
% set(insetAxes,'YScale','log');
% grid on;
% insetAxes.XGrid = 'on'; 
% insetAxes.YGrid = 'on'; 
% insetAxes.MinorGridAlpha = 0;  
% insetAxes.XAxis.FontSize = 8;
% insetAxes.YAxis.FontSize = 8;
% yticks([1, 1e1, 1e2, 1e3, 1e4, 1e5]);
% yticklabels({'1','10','10^{2}','10^{3}','10^{4}','10^{5}'});
% xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
% ylabel('Frequency','FontSize',8)
% 
% exportgraphics(gcf,fullfile('.','figures','SmartBuoy_CEFAS_histogram.png'),'Resolution',600)
% 
% %% Histogram by season
% 
% seasons = {'Winter','Spring','Summer','Autumn'};
% colourScheme = brewermap(length(seasons),'*Spectral');
% 
% figure()
% set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.50],'Color','w') 
% haxis = zeros(2,2);
% 
% for iSeason = 1:length(seasons)
%     
%     haxis(iSeason) = subaxis(2,2,iSeason,'Spacing',0.018,'Padding',0.020,'Margin',0.07);
%     ax(iSeason).pos = get(haxis(iSeason),'Position');
%     if (iSeason == 1 || iSeason == 2)
%         ax(iSeason).pos(2) = ax(iSeason).pos(2) + 0.05;
%     end
%     ax(iSeason).pos(1) = ax(iSeason).pos(1) + 0.01;
%     set(haxis(iSeason),'Position',ax(iSeason).pos) 
%  
%     switch seasons{iSeason}
%         case 'Winter'
%             seasonData = SB.value(SB.season == "Winter" & SB.variable == "CHLOR"); 
%         case 'Spring'
%             seasonData = SB.value(SB.season == "Spring" & SB.variable == "CHLOR"); 
%         case 'Summer'
%             seasonData = SB.value(SB.season == "Summer" & SB.variable == "CHLOR"); 
%         case 'Autumn'
%             seasonData = SB.value(SB.season == "Autumn" & SB.variable == "CHLOR"); 
%     end
%     
%     % Main plot
%     
%     % Plot all data
%     histogram(haxis(iSeason),SB.value(SB.variable == "CHLOR"),50,'BinLimits',[0,20],'FaceColor',[0.7,0.7,0.7],'EdgeColor','none')
%     hold on
%     % Plot seasonal data
%     histogram(haxis(iSeason),seasonData,50,'BinLimits',[0,20],'FaceColor',colourScheme(iSeason,:))
%     hold on
%     ylim([1 6e4])
%     set(gca,'YScale','log');
%     grid on;
%     mainSeasonalAxes = gca;
%     mainSeasonalAxes.XGrid = 'on';
%     mainSeasonalAxes.YGrid = 'on';
%     mainSeasonalAxes.MinorGridAlpha = 0;
% 
%     if (iSeason == 1 || iSeason == 3)
%         ylabel('Frequency');
%     end
%     if (iSeason == 3 || iSeason == 4)
%         xlabel('Chlorophyll a (mg m^{-3})');
%     end
%     title(seasons{iSeason});
%     
%     % Inset plot
%     insetSeasonalAxes = axes('Position', [ax(iSeason).pos(1) + ax(iSeason).pos(3)*0.45,... 
%                                           ax(iSeason).pos(2) + ax(iSeason).pos(4)*0.47,... 
%                                           0.52*ax(iSeason).pos(3),... 
%                                           0.49*ax(iSeason).pos(4)]);
%     
%     histogram(insetSeasonalAxes,seasonData,19,'BinLimits',[1e-1,2],'FaceColor',colourScheme(iSeason,:))
%     hold on
%     ylim([1 1e4])
%     set(gca,'YScale','log');
%     grid on;
%     insetSeasonalAxes.XGrid = 'on';
%     insetSeasonalAxes.YGrid = 'on';
%     insetSeasonalAxes.MinorGridAlpha = 0;
%     insetSeasonalAxes.XAxis.FontSize = 8;
%     insetSeasonalAxes.YAxis.FontSize = 8;
%     yticks([1, 1e1, 1e2, 1e3, 1e4]);
%     yticklabels({'1','10','10^{2}','10^{3}','10^{4}'});
%     xlabel('Chlorophyll a (mg m^{-3})','FontSize',8);
%     ylabel('Frequency','FontSize',8)  
%     
% end
% hold off
% 
% exportgraphics(gcf,fullfile('.','figures','SmartBuoy_CEFAS_histogram_byseason.png'),'Resolution',600)
% clear ax
% 
% % end
% 
%  