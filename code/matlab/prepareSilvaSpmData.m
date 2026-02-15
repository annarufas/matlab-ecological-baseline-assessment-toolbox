function [SPMclimatol] = prepareSilvaSpmData(pathAreaStudyShapefile)

% PREPAREINSITUPHDATA Read in Suspended particulate matter (SPM) data
% Monthly average non-algal Suspended Particulate Matter concentrations on the UK shelf waters
% https://data.cefas.co.uk/view/18133
%
%   INPUT:
%       pathAreaStudyShapefile - shapefile with our area of study
%
%   OUTPUT:
%       SPMclimatol            - Matlab table with combined data from both programs and filtered pH data
%          
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 24 April 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Definition of area of study boundaries

areaStudy = m_shaperead(pathAreaStudyShapefile); % lat/lon coordinates

% Box of 100 x 100 km around the Endurance site
minLatAreaStudy = areaStudy.MBRy(1); % 53.758
maxLatAreaStudy = areaStudy.MBRy(2); % 54.656
minLonAreaStudy = areaStudy.MBRx(1); % 0.265
maxLonAreaStudy = areaStudy.MBRx(2); % 1.801

%% Read in data by individual year

pathDataIndividualYearDir = fullfile('.','data','raw','Silva_2016_SPM',...
    'analysed_spim_individual_years');

% Get a list of all files in the folder
ncFiles = dir(fullfile(pathDataIndividualYearDir, '*.nc'));
years = [1998:1:2015];
nYears = length(years);

% Inspect the first file
S = ncinfo(fullfile(pathDataIndividualYearDir, ncFiles(1).name)); % short summary
latSpm = double(ncread(fullfile(pathDataIndividualYearDir, ncFiles(1).name),'lat'));
lonSpm = double(ncread(fullfile(pathDataIndividualYearDir, ncFiles(1).name),'lon'));

spmAnnualValues = NaN(length(lonSpm),length(latSpm),nYears);
nSamplesSpm = zeros(length(lonSpm),length(latSpm),nYears);
for iYear = 1:nYears
    currentFile = ncFiles(iYear).name;
    pathFile = fullfile(pathDataIndividualYearDir, currentFile);
    spmAnnualValues(:,:,iYear) = ncread(pathFile,'SMMean');
    nSamplesSpm(:,:,iYear) = ncread(pathFile,'SMNSam');
end

% Swap latitude and longitude
spmAnnualValues = permute(spmAnnualValues, [2 1 3]);

% figure(); pcolor(spmAnnualValues(:,:,1)); shading interp

% Adjust to desired lat range and lon range – this reduces the size of
% the output array
[~,closestIdxMinLat] = min(abs(latSpm-minLatAreaStudy));
[~,closestIdxMaxLat] = min(abs(latSpm-maxLatAreaStudy));
[~,closestIdxMinLon] = min(abs(lonSpm-minLonAreaStudy));
[~,closestIdxMaxLon] = min(abs(lonSpm-maxLonAreaStudy));

% To be on the safe side, expand
idxMinLat = closestIdxMinLat-1;
idxMaxLat = closestIdxMaxLat+1;
idxMinLon = closestIdxMinLon-1;
idxMaxLon = closestIdxMaxLon+1;

latSpmAnnual = latSpm(idxMinLat:idxMaxLat);
lonSpmAnnual = lonSpm(idxMinLon:idxMaxLon);
spmAnnualValues = spmAnnualValues(idxMinLat:idxMaxLat,idxMinLon:idxMaxLon,:);

% figure(); pcolor(spmAnnualValues(:,:,1)); shading interp

%% Read in data by individual month in each year (files corrupted!)

pathDataIndividualMonthDir = fullfile('.','data','raw','Silva_2016_SPM',...
    'silva_analaysed_spim_individual_months');

% zipFiles = dir(fullfile(pathDataIndividualMonthDir, '*.gz'));
% for iFile = 1:length(zipFiles)
%     zipFilePath = fullfile(pathDataIndividualMonthDir, zipFiles(iFile).name);
%     gunzip(zipFilePath, pathDataIndividualMonthDir);
% end

% Get a list of all files in the folder
ncFiles = dir(fullfile(pathDataIndividualMonthDir, '*.nc'));

% Inspect the first file
S = ncinfo(fullfile(pathDataIndividualMonthDir, ncFiles(1).name)); % short summary
latSpm = ncread(fullfile(pathDataIndividualMonthDir, ncFiles(1).name),'lat');
lonSpm = ncread(fullfile(pathDataIndividualMonthDir, ncFiles(1).name),'lon');

spmMonthlyValues = zeros(length(lonSpm),length(latSpm),12); % '1e-3 kg m-3'
nSamplesSpm = zeros(length(lonSpm),length(latSpm),12);
spmStd = zeros(length(lonSpm),length(latSpm),12);
for iMonth = 1:12
    currentFile = ncFiles(iMonth).name;
    pathFile = fullfile(pathDataIndividualMonthDir, currentFile);
    spmMonthlyValues(:,:,iMonth) = ncread(pathFile,'SMMean');
    nSamplesSpm(:,:,iMonth) = ncread(pathFile,'SMNSam');
    spmStd(:,:,iMonth) = ncread(pathFile,'SMStDev'); 
end

%% Read in monthly climatology data

pathDataMonthlyClimatologyDir = fullfile('.','data','raw','Silva_2016_SPM',...
    'analysed_spim_monthly_climatology_1998_2015');

% Get a list of all files in the folder
ncFiles = dir(fullfile(pathDataMonthlyClimatologyDir, '*.nc'));

% Inspect the first file
S = ncinfo(fullfile(pathDataMonthlyClimatologyDir, ncFiles(1).name)); % short summary
latSpm = double(ncread(fullfile(pathDataMonthlyClimatologyDir, ncFiles(1).name),'lat'));
lonSpm = double(ncread(fullfile(pathDataMonthlyClimatologyDir, ncFiles(1).name),'lon'));

spmMonthlyClimatology = zeros(length(lonSpm),length(latSpm),12); % '1e-3 kg m-3'
nSamplesSpm = zeros(length(lonSpm),length(latSpm),12);
spmStd = zeros(length(lonSpm),length(latSpm),12);
for iMonth = 1:12
    currentFile = ncFiles(iMonth).name;
    pathFile = fullfile(pathDataMonthlyClimatologyDir, currentFile);
    spmMonthlyClimatology(:,:,iMonth) = ncread(pathFile,'SMMean');
    nSamplesSpm(:,:,iMonth) = ncread(pathFile,'SMNSam');
    spmStd(:,:,iMonth) = ncread(pathFile,'SMStDev'); 
end

% Swap latitude and longitude
spmMonthlyClimatology = permute(spmMonthlyClimatology, [2 1 3]);

% Adjust to desired lat range and lon range – this reduces the size of
% the output array
[~,closestIdxMinLat] = min(abs(latSpm-minLatAreaStudy));
[~,closestIdxMaxLat] = min(abs(latSpm-maxLatAreaStudy));
[~,closestIdxMinLon] = min(abs(lonSpm-minLonAreaStudy));
[~,closestIdxMaxLon] = min(abs(lonSpm-maxLonAreaStudy));

% To be on the safe side, expand
idxMinLat = closestIdxMinLat-1;
idxMaxLat = closestIdxMaxLat+1;
idxMinLon = closestIdxMinLon-1;
idxMaxLon = closestIdxMaxLon+1;

latSpm = latSpm(idxMinLat:idxMaxLat);
lonSpm = lonSpm(idxMinLon:idxMaxLon);
SPMclimatol = spmMonthlyClimatology(idxMinLat:idxMaxLat,idxMinLon:idxMaxLon,:);

%% We will output the monthly climatology data

save(fullfile('.','data','processed','silvaSPM.mat'),...
    'SPMclimatol','latSpm','lonSpm')

end