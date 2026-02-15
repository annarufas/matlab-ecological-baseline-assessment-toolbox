function plotOpticalWaterTypes(yearOfChoice,ASoccci,pathAreaStudyShapefile,...
    pathEnduranceShapefile)

% PLOTOPTICALWATERTYPES Plot optical water types from the use OC-CCI 1 km 
% spatial resolution data product.
%
%   INPUT:
%       yearOfChoice           - year of choice to plot optical water types
%       ASoccci                - table with time-series data from OC-CCI for our area of study
%       pathAreaStudyShapefile - shapefile with our area of study
%       pathEnduranceShapefile - shapefile with the Endurance GCS site
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

%% Extract water class data from the OC-CCI daily 1 km resolution product

idxProd = find(strcmp({ASoccci.ID}, 'occci_1km_1day'));
idxsWaterClassVars = find(contains(ASoccci(idxProd).varNames, 'water_class'));
waterClass.dataset = ASoccci(idxProd).dataset(:,:,:,idxsWaterClassVars);
waterClass.time = ASoccci(idxProd).time; 
waterClass.lat = ASoccci(idxProd).lat; 
waterClass.lon = ASoccci(idxProd).lon;

waterClass.nLonPixels = length(waterClass.lon);
waterClass.nLatPixels = length(waterClass.lat);
nWaterClasses = numel(idxsWaterClassVars);

%% Calculate the per pixel dominant water class

% Generate 12 dates (one for each month of the year of choice) centred
% around the 15th of each month
dates = datetime(yearOfChoice,1,15):calmonths(1):datetime(yearOfChoice,12,15); 

% Build a 5-day composite scene around the central month's date
monthlySceneWaterClass = NaN(waterClass.nLatPixels,waterClass.nLonPixels,12);
for iMonth = 1:12
    
    centralDate = dates(iMonth);
    forwardDates = centralDate + days(1:2); % add 2 days
    backwardDates = centralDate - days(1:2); % subtract 2 days
    fiveDayString = [fliplr(backwardDates), centralDate, forwardDates];
    
    % Find multiple dates in the time array
    whereIsIt = ismember(waterClass.time,fiveDayString); % "true" if the element in "a" is a member of "b".
    idxsToDates = find(whereIsIt);
    
    % Create the water class 5-day composite
    thisDateImage = squeeze(waterClass.dataset(:,:,idxsToDates,:));
    monthlyWaterClassFiveDayMean = zeros(waterClass.nLatPixels,waterClass.nLonPixels,nWaterClasses);
    for iWaterClass = 1:nWaterClasses
        for iLon = 1:waterClass.nLonPixels
            for iLat = 1:waterClass.nLatPixels
                monthlyWaterClassFiveDayMean(iLat,iLon,iWaterClass) = mean(squeeze(thisDateImage(iLat,iLon,:,iWaterClass)),'omitnan');
            end
        end
    end
    
    % Calculate the fractional contribution of each water class
    monthlyWaterClassFrac = zeros(waterClass.nLatPixels,waterClass.nLonPixels,nWaterClasses);
    for iLon = 1:waterClass.nLonPixels
        for iLat = 1:waterClass.nLatPixels
            sumWaterClass = sum(squeeze(monthlyWaterClassFiveDayMean(iLat,iLon,:)),'omitnan'); 
            for iWaterClass = 1:nWaterClasses
                monthlyWaterClassFrac(iLat,iLon,iWaterClass) = monthlyWaterClassFiveDayMean(iLat,iLon,iWaterClass)/sumWaterClass;
            end
        end
    end
    
    % Calculate the dominant water class
    monthlyWaterClassDominant = zeros(waterClass.nLatPixels,waterClass.nLonPixels);
    for iLon = 1:waterClass.nLonPixels
        for iLat = 1:waterClass.nLatPixels
            arrayPixelWaterClasses = squeeze(monthlyWaterClassFrac(iLat,iLon,:));
            if (~isnan(sum(arrayPixelWaterClasses)))
                [~,idxMax] = max(arrayPixelWaterClasses,[],'omitnan');
                monthlyWaterClassDominant(iLat,iLon) = idxMax;
            end
        end
    end
    
    monthlySceneWaterClass(:,:,iMonth) = monthlyWaterClassDominant;

end

% Set to NaN those pixels that could not be assigned a water class
monthlySceneWaterClass(monthlySceneWaterClass == 0) = NaN;

%% Plot

areaStudy = m_shaperead(pathAreaStudyShapefile); % lat/lon coordinates
enduranceSite = m_shaperead(pathEnduranceShapefile); % lat/lon coordinates

monthsTag = month(dates,'shortname');
datesTag = cell(1,12);
for iMonth = 1:12
%     datesTag{iMonth} = [monthsTag{iMonth},' ','12-16',' ',num2str(yearOfChoice)];
    datesTag{iMonth} = [monthsTag{iMonth},' ',num2str(yearOfChoice)];
end

figure()
set(gcf,'Units','Normalized','Position',[0.01 0.05 0.45 0.60],'Color','w') 
haxis = zeros(3,4);

for iMonth = 1:12
    
    haxis(iMonth) = subaxis(3,4,iMonth,'Spacing',0.015,'Padding',0.015,'Margin',0.050);
    ax(iMonth).pos = get(haxis(iMonth),'Position');
    ax(iMonth).pos(2) = ax(iMonth).pos(2) + 0.050;
    set(haxis(iMonth),'Position',ax(iMonth).pos)
    
    m_proj('equidistant','long',areaStudy.MBRx,'lat',areaStudy.MBRy); 
    m_pcolor(waterClass.lon,waterClass.lat,monthlySceneWaterClass(:,:,iMonth)) 
    shading flat
    cmap = brewermap(nWaterClasses,'*RdYlBu');
    colormap(cmap)
    caxis([1 nWaterClasses])
    
    m_grid('linewi',1,'tickdir','out','FontSize',9);
    m_line(enduranceSite.ncst{1}(:,1),enduranceSite.ncst{1}(:,2),'linewi',1,'color','k');
    
    title(datesTag{iMonth},'FontSize',12);
    
end

cb = colorbar('Location','southoutside');

% Set tick marks to be middle of each range
dTk = diff(cb.Limits)/(2*length(cmap));
set(cb,'Ticks',[cb.Limits(1)+dTk:2*dTk:cb.Limits(2)-dTk],...
    'TickLabels',({'1','2','3','4','5','6','7','8','9','10','11','12','13','14'}))
cb.Position(1) = cb.Position(1) - 0.50;
cb.Position(2) = cb.Position(2) - 0.105;
cb.Position(3) = 0.50; % LENGTH
cb.Position(4) = 0.03; % WIDTH
cb.Label.String = 'Water class number'; 
cb.FontSize = 11;

exportgraphics(gcf,fullfile('.','figures',strcat('optical_water_types',num2str(yearOfChoice),'.png')),'Resolution',600)

end % plotOpticalWaterTypes