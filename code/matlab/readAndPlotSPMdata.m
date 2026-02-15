

%% Suspended particulate matter (SPM) data
% Monthly average non-algal Suspended Particulate Matter concentrations on the UK shelf waters
% https://data.cefas.co.uk/view/18133

clc
clear all

fullPathMainDir         = '/Users/Anna/LocalDocuments/Academic/Projects/Agile/baseline_assessment/';
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

%% By individual year

fullPathDataDir = strcat(fullPathMainDir,'CEFAS_SPM_Silva/analysed_spim_individual_years/');

% Get a list of all files in the folder
ncFiles = dir(fullfile(fullPathDataDir, '*.nc'));
years = [1998:1:2015];
nYears = length(years);

% Inspect the first file
S = ncinfo(fullfile(fullPathDataDir, ncFiles(1).name)); % short summary
latSpmTemp = double(ncread(fullfile(fullPathDataDir, ncFiles(1).name),'lat'));
lonSpmTemp = double(ncread(fullfile(fullPathDataDir, ncFiles(1).name),'lon'));

spmAnnualClimatologyMeanTemp = NaN(length(lonSpmTemp),length(latSpmTemp),nYears);
nSamplesSpm = zeros(length(lonSpmTemp),length(latSpmTemp),nYears);
for iYear = 1:nYears
    currentFile = ncFiles(iYear).name;
    fullPathFile = fullfile(fullPathDataDir, currentFile);
    spmAnnualClimatologyMeanTemp(:,:,iYear) = ncread(fullPathFile,'SMMean');
    nSamplesSpm(:,:,iYear) = ncread(fullPathFile,'SMNSam');
end

% Swap latitude and longitude
spmAnnualClimatologyMeanTemp = permute(spmAnnualClimatologyMeanTemp, [2 1 3]);

% figure()
% pcolor(spmAnnualClimatologyMeanTemp(:,:,1))
% shading interp

% Adjust to desired lat range and lon range – this reduces the size of
% the output array
[~,closestIdxMinLat] = min(abs(latSpmTemp-minLatAreaStudy));
[~,closestIdxMaxLat] = min(abs(latSpmTemp-maxLatAreaStudy));
[~,closestIdxMinLon] = min(abs(lonSpmTemp-minLonAreaStudy));
[~,closestIdxMaxLon] = min(abs(lonSpmTemp-maxLonAreaStudy));

% To be on the safe side, expand
idxMinLat = closestIdxMinLat-1;
idxMaxLat = closestIdxMaxLat+1;
idxMinLon = closestIdxMinLon-1;
idxMaxLon = closestIdxMaxLon+1;

latSpmAnnual = latSpmTemp(idxMinLat:idxMaxLat);
lonSpmAnnual = lonSpmTemp(idxMinLon:idxMaxLon);
spmAnnualClimatologyMean = spmAnnualClimatologyMeanTemp(idxMinLat:idxMaxLat,idxMinLon:idxMaxLon,:);

% figure()
% pcolor(spmAnnualClimatologyMean(:,:,1))
% shading interp

%% By individual month in each year (files corrupted!!!!)

fullPathDataDir = strcat(fullPathMainDir,'CEFAS_SPM_Silva/silva_analaysed_spim_individual_months/');

% zipFiles = dir(fullfile(fullPathDataDir, '*.gz'));
% for iFile = 1:length(zipFiles)
%     zipFilePath = fullfile(fullPathDataDir, zipFiles(iFile).name);
%     gunzip(zipFilePath, fullPathDataDir);
% end

% Get a list of all files in the folder
ncFiles = dir(fullfile(fullPathDataDir, '*.nc'));

% Inspect the first file
S = ncinfo(fullfile(fullPathDataDir, ncFiles(1).name)); % short summary
latSpmTemp = ncread(fullfile(fullPathDataDir, ncFiles(1).name),'lat');
lonSpmTemp = ncread(fullfile(fullPathDataDir, ncFiles(1).name),'lon');

spmMonthlyClimatologyMeanTemp = zeros(length(lonSpmTemp),length(latSpmTemp),12); % '1e-3 kg m-3'
nSamplesSpm = zeros(length(lonSpmTemp),length(latSpmTemp),12);
spmStd = zeros(length(lonSpmTemp),length(latSpmTemp),12);
for iMonth = 1:12
    currentFile = ncFiles(iMonth).name;
    fullPathFile = fullfile(fullPathDataDir, currentFile);
    spmMonthlyClimatologyMeanTemp(:,:,iMonth) = ncread(fullPathFile,'SMMean');
    nSamplesSpm(:,:,iMonth) = ncread(fullPathFile,'SMNSam');
    spmStd(:,:,iMonth) = ncread(fullPathFile,'SMStDev'); 
end

%% Monthly climatology

fullPathDataDir = strcat(fullPathMainDir,'CEFAS_SPM_Silva/analysed_spim_monthly_climatology_1998_2015/');

% Get a list of all files in the folder
ncFiles = dir(fullfile(fullPathDataDir, '*.nc'));

% Inspect the first file
S = ncinfo(fullfile(fullPathDataDir, ncFiles(1).name)); % short summary
latSpmTemp = double(ncread(fullfile(fullPathDataDir, ncFiles(1).name),'lat'));
lonSpmTemp = double(ncread(fullfile(fullPathDataDir, ncFiles(1).name),'lon'));

spmMonthlyClimatologyMeanTemp = zeros(length(lonSpmTemp),length(latSpmTemp),12); % '1e-3 kg m-3'
nSamplesSpm = zeros(length(lonSpmTemp),length(latSpmTemp),12);
spmStd = zeros(length(lonSpmTemp),length(latSpmTemp),12);
for iMonth = 1:12
    currentFile = ncFiles(iMonth).name;
    fullPathFile = fullfile(fullPathDataDir, currentFile);
    spmMonthlyClimatologyMeanTemp(:,:,iMonth) = ncread(fullPathFile,'SMMean');
    nSamplesSpm(:,:,iMonth) = ncread(fullPathFile,'SMNSam');
    spmStd(:,:,iMonth) = ncread(fullPathFile,'SMStDev'); 
end

% Swap latitude and longitude
spmMonthlyClimatologyMeanTemp = permute(spmMonthlyClimatologyMeanTemp, [2 1 3]);

% Adjust to desired lat range and lon range – this reduces the size of
% the output array
[~,closestIdxMinLat] = min(abs(latSpmTemp-minLatAreaStudy));
[~,closestIdxMaxLat] = min(abs(latSpmTemp-maxLatAreaStudy));
[~,closestIdxMinLon] = min(abs(lonSpmTemp-minLonAreaStudy));
[~,closestIdxMaxLon] = min(abs(lonSpmTemp-maxLonAreaStudy));

% To be on the safe side, expand
idxMinLat = closestIdxMinLat-1;
idxMaxLat = closestIdxMaxLat+1;
idxMinLon = closestIdxMinLon-1;
idxMaxLon = closestIdxMaxLon+1;

latSpm = latSpmTemp(idxMinLat:idxMaxLat);
lonSpm = lonSpmTemp(idxMinLon:idxMaxLon);
spmMonthlyClimatologyMean = spmMonthlyClimatologyMeanTemp(idxMinLat:idxMaxLat,idxMinLon:idxMaxLon,:);

%% Save

save(strcat(fullPathMainDir,'spmSilva.mat'),...
    'spmAnnualClimatologyMean','latSpmAnnual','lonSpmAnnual',...
    'spmMonthlyClimatologyMean','latSpm','lonSpm')
