function [filteredTable] = addBathymetryAndFilter(dataTable,lonColName,latColName,...
    pathBathymetryFile,AVG_SHELF_DEPTH)

% ADDBATHYMETRYANDFILTER Assign bathymetry values to data points based on  
% latitude and longitude coordinates using bathymetry data. Filter out data 
% points that are not within the specified shelf sea depth range.
%
%   INPUT:
%       dataTable          - data table
%       lonColName         - name of the column containing longitude values
%       latColName         - name of the column containing latitude values
%       pathBathymetryFile - .nc file containing the bathymetric dataset
%       AVG_SHELF_DEPTH    - depth in m
% 
%   OUTPUT:
%       filteredTable      - filtered data table
% 
%   WRITTEN BY A. RUFAS, UNIVERISTY OF OXFORD
%   Anna.RufasBlanco@earth.ox.ac.uk
%
%   Version 1.0 - Completed 3 May 2024   
%
% =========================================================================
%%
% -------------------------------------------------------------------------
% PROCESSING STEPS
% -------------------------------------------------------------------------

%% Load bathymetry data

% ncdisp(pathBathymetryFile)
bathy = ncread(pathBathymetryFile, 'z');
lat = ncread(pathBathymetryFile, 'lat');
lon = ncread(pathBathymetryFile, 'lon');

%% Create a gridded interpolant for bathymetry data

[Xbathy, Ybathy] = ndgrid(lon, lat);
Fbathy = griddedInterpolant(Xbathy, Ybathy, bathy);

%% Process each data point

% Initialise logical array for filtering data points
nDataPoints = height(dataTable);
isDataPointShelfSea = false(nDataPoints, 1);

for iDataPoint = 1:nDataPoints
    qLon = dataTable.(lonColName)(iDataPoint);
    qLat = dataTable.(latColName)(iDataPoint);
    % Query points for interpolation
    qDataPointDepth = Fbathy(qLon, qLat);
    % Store bathymetry value in the data table
    dataTable.bathymetry_m(iDataPoint) = qDataPointDepth;
    % Check if data point is in the shelf sea
    if (qDataPointDepth < 0 && abs(qDataPointDepth) <= AVG_SHELF_DEPTH)
        isDataPointShelfSea(iDataPoint) = true;
    else
        fprintf('\nData point #%d, with depth %3.2f will be deleted.\n', iDataPoint, qDataPointDepth);
    end
end

% Calculate the number of data points to be deleted
nLocsToDelete = sum(~isDataPointShelfSea);
fprintf('\n%d data points that were not on the shelf seas have been deleted\n', nLocsToDelete);

% Filter data points to keep only the ones on the shelf sea
filteredTable = dataTable(isDataPointShelfSea, :);

end % addBathymetryAndFilter
